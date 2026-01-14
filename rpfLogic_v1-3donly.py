import pandas as pd
import numpy as np
import re
from Bio.PDB import PDBParser, NeighborSearch
from itertools import combinations, product
from scipy.spatial.distance import cdist
from ccpn.core.Project import Project
from collections import defaultdict
from scipy.spatial import KDTree

"""
 **Recall, Precision & F-Score; validating structures with peak lists**

 The purpose of this popup is to asses the quality of structures and structure
 generating peak lists by calculating recall, precision, and F-measure (RPF)
 scores, according to the method of Huang, Powers and Montelione (see reference
 below). Put simply, the calculation uses short distances from a structure to
 find peak locations which would be expected but which don't appear in the
 peak list; "missing peaks" and also to find peaks which do not correspond to a
 short distance in the  structure; "unexplained peaks". These results are then
 combined to form the overall RPF metrics, which gives a convenient and
 informative quality measure to say how well peak list and structure match. If
 there are good, real peaks that don't match the structure then something may
 be wrong with the structure calculation/restraints. If there are short
 distances with with no corresponding peak then the spectra should be
 investigated to find out why, which often is actually legitimate, but worth
 checking all the same. For full details read the paper!

 To perform the analysis the user simply selects the required peak lists in the
 project, the required structure and clicks *[Run RPF]*.

 **Peak Lists**

 The main peak lists table is where the user selects which peak lists in the
 current project to consider.
 Naturally these will be those peak lists that represent through-space (NOESY)
 connectivities which actually relate to the structure being analysed. Also,
 there are a some options that control how peak lists are investigated.
 Specifically, the user can exclude diagonal peaks and redundant prochrials
 (e.g. methylene protons that are too close in chemical shift to give separate
 peaks are not really missing) by setting PPM thresholds.

 **Reference**

 *Huang, Y. J.; Powers, R. & Montelione, G. T. Protein NMR recall, precision,
 and F-measure scores (RPF scores): structure quality assessment measures based
 on information retrieval statistics. J Am Chem Soc 127, 1665-74 (2005)*

 """


class RpfCalculator:
    NOE_POWER = -6  # For distance^-6 weighting

    # This dictionary says how different numbers of 13C, 15N bound
    # dimensions limits which protons would be visible in a spectrum
    # given 13C or 15N connectivity
    # Context is ('N','N),('C','N'),('C','C') for 2D NOESY i.e. any 1H to any 1H
    # ('N','N),('C','N') for 3D HSQC-NOESY, ('C','C') for 4D HCCH NOESY etc..

    CN_DIM_CONTEXT_DICT = {
        (0, 0): (('C', 'C'), ('C', 'N'), ('N', 'C'), ('N', 'N')),
        (1, 0): (('C', 'C'), ('C', 'N'), ('N', 'C')),
        (0, 1): (('C', 'N'), ('N', 'C'), ('N', 'N')),
        (1, 1): (('C', 'N'), ('N', 'C')),
        (2, 0): (('C', 'C'),),
        (0, 2): (('N', 'N'),),
        }

    def __init__(self, project: Project, pdb_path: str, peak_list_pids: list,
                 chemical_shift_pid: str, distance_threshold: float = 4.8, prochiral_tolerance=0.04):
        self.project = project
        self.pdb_path = pdb_path
        self.peak_list_pids = peak_list_pids
        self.chemical_shift_pid = chemical_shift_pid
        self.distance_threshold = distance_threshold
        self.prochiral_tolerance = prochiral_tolerance
        # Data storage

        self.proton_dists_conn = {}
        self.peaklists = []
        self.ensemble = None
        self.ensembleModels = []
        self.proton_pairs = []
        self.shift_map = None
        self.chain = None

        # dict of dicts peaklist.pid as key
        self.shift_pairs = {}
        self.matched_peaks = {}
        self.unexplained_peaks = {}
        self.missing_pairs = {}

        self.residue_scores = {}

        self.report = {}
        self.unmatched_peaks_report = {}
        self.metrics = {}
        self.syntheticPeaks = {}

    def getPeakList_objects(self):
        return [self.project.getByPid(pid) for pid in self.peak_list_pids]

    def _get_atom_context(self, atom1):
        """Get H atom context from structure"""
        model = self.ensembleModels[0]
        ns = NeighborSearch(list(model.get_atoms()))
        h1_neighbors = ns.search(atom1.coord, 1.2, level='A')
        h1_neighbors = [a for a in h1_neighbors if a != atom1]
        if not h1_neighbors :
            return None

        e1 = h1_neighbors[0].element

        return e1

    def get_ambig_noe_conn(self, shift_pairs, peak_df, tol_dict, diagonal_tol=0.3):
        """Enhanced ambiguous peak-resonance mapping with v2 features"""
        connections = defaultdict(list)
        peak_assignments = defaultdict(list)

        # Validate and format points array
        if not shift_pairs:
            print("⚠️ No shift pairs available for KDTree construction")
            return connections, peak_assignments

        try:
            # Create 2D array of shape (N, 2)
            points = np.array([[pair['shift1'], pair['shift2']] for pair in shift_pairs])
            print(f"✅ KDTree points shape: {points.shape}")  # Should be (N, 2)

            # Build KDTree with periodic boundary conditions for aliasing
            tree = KDTree(points, boxsize=[360, 360])  # Assuming ±180ppm aliasing

        except Exception as e:
            print(f"❌ KDTree creation failed: {str(e)}")
            print(f"Sample points: {points[:3] if len(points) > 3 else points}")
            return connections, peak_assignments
        for idx, peak in peak_df.iterrows():
            if pd.isna(peak['position_F1']) or pd.isna(peak['position_F2']):
                continue

            # Get dynamic tolerance from peak dimensions
            # there is an error here..
            #tol1 = tol_dict[peak['dim_F1']][0]
            #tol2 = tol_dict[peak['dim_F2']][0]
            # ... easy temporary fix....
            tol1 = tol_dict[1][0]
            tol2 = tol_dict[2][0]
            query_tol = np.sqrt(tol1**2 + tol2**2)

            # Find matches considering aliasing
            matches = tree.query_ball_point([peak['position_F1'], peak['position_F2']], query_tol)

            # Exclude diagonal peaks
            if abs(peak['position_F1'] - peak['position_F2']) < diagonal_tol:
                continue

            # Store ambiguous connections
            for i in matches:
                pair = shift_pairs[i]
                connections[(pair['atom1'], pair['atom2'])].append(peak)
                peak_assignments[peak['serial']].append(pair)

        return connections, peak_assignments

    def get_chemicalShift_of_attached_atom(self, atom):
        model = self.ensembleModels[0]
        chainCode = atom.parent.parent.id  # Chain
        resiNum = atom.parent.id[1]  # Residue number
        #atomName = atom.name  # Atom name
        chain = model[chainCode]  # Get chain A
        residue = chain[(" ", resiNum, " ")]  # (hetflag, seqnum, insertion_code)
        resiName = residue.resname

        model = self.ensembleModels[0]
        ns = NeighborSearch(list(model.get_atoms()))
        h1_neighbors = ns.search(atom.coord, 1.2, level='A')
        h1_neighbors = [a for a in h1_neighbors if a != atom]
        e1 = h1_neighbors[0].element
        atomName = h1_neighbors[0].get_name()

        NmrAtomPID = f"NA:{chainCode}.{resiNum}.{resiName.upper()}.{atomName}"

        try:
            nmrAtom = self.project.getByPid(NmrAtomPID)
            chemicalShift = nmrAtom.chemicalShifts[0]
            return chemicalShift.value, e1, h1_neighbors[0]
        except:
            return None, e1, h1_neighbors[0]

    def get_proton_dists_conn(self, contexts, dist_threshold=4.8):
        """NOE-weighted distance calculation across ensemble
        This calculation is very slow - probably revisit at some point - perhaps
        using a biopython neerest neighbour search to speed up?
        """
        noe_dict = defaultdict(float)

        # Original contexts example
        # contexts = (('N', 'C'), ('N', 'C'))

        # Flatten and get unique elements
        unique_elements = list(set([element for context in contexts for element in context]))

        # find protons attached to the atoms in the contexts within a bond distance.
        # done like this as PDB contains no connectivity info
        protons = self.get_bound_protons(search_element='H',
                                         target_elements=unique_elements,
                                         bond_cutoff=1.2)

        for model in self.ensembleModels:  # ACTUALLY ITERATE THROUGH MODELS
            # Get coordinates FOR THIS MODEL
            coords = []
            for atom in protons:
                # Get atom from current model using full_id
                # Get hierarchy components from the original atom
                chain_id = atom.parent.parent.id
                residue_id = atom.parent.id
                atom_name = atom.name

                chain = model[chain_id]
                residue = chain[residue_id]
                model_atom = residue[atom_name]
                coords.append(model_atom.coord)

            coords = np.array(coords)
            dist_matrix = cdist(coords, coords)
            print(dist_matrix)

            for i, j in combinations(range(len(protons)), 2):
                if dist_matrix[i, j] > dist_threshold:
                    continue

                a1, a2 = protons[i], protons[j]

                ctx_a1 = self._get_atom_context(a1)
                ctx_a2 = self._get_atom_context(a2)
                ctx = (ctx_a1, ctx_a2)

                if ctx in contexts:
                    noe = dist_matrix[i, j]**self.NOE_POWER
                    noe_dict[(a1, a2)] += noe

        # Average across ALL models
        num_models = len(list(self.ensembleModels))
        return {key: value / num_models for key, value in noe_dict.items()}

    def _valid_context(self, proton, neighbor, contexts):
        """
        Check if the bonding context of two protons matches any allowed context.

        Args:
            proton (Atom): First hydrogen atom in the pair.
            neighbor (Atom): Second hydrogen atom in the pair.
            contexts (tuple): Allowed bonding contexts (e.g., ('C','C'), ('N','C').

        Returns:
            bool: True if the pair matches any allowed context (order-agnostic).
        """
        # Get bonded heavy atoms for both protons
        ctx1 = self._get_atom_context(proton)
        ctx2 = self._get_atom_context(neighbor)

        if not ctx1 or not ctx2:
            return False  # Invalid if no bonded atom found

        # check if context in contexts
        return (ctx1, ctx2) in contexts

    def broken_get_proton_dists_conn(self, contexts, dist_threshold="4.8"):
        # dist_threshold: distance threshold, in the structure, below which resonance pairs are expected to result in a peak
        # contexts: atom contexts that are detectable in the experiment. Example = (('N', 'C'), ('N', 'C'))
        noe_dict = defaultdict(float)

        # Flatten and get unique elements
        unique_elements = list(set([element for context in contexts for element in context]))

        protons = self.get_bound_protons(search_element='H',
                                         target_elements=unique_elements,
                                         bond_cutoff=1.2)

        for model in self.ensembleModels:
            ns = NeighborSearch(list(model.get_atoms()))
            # Find all protons within dist_threshold
            for proton in protons:
                neighbors = ns.search(proton.coord, dist_threshold, level='A')
                for neighbor in neighbors:
                    if self._valid_context(proton, neighbor, contexts):
                        dist = np.linalg.norm(proton.coord - neighbor.coord)
                        noe_dict[(proton, neighbor)] += dist**self.NOE_POWER
        # Average across models
        return {k: v / len(self.ensembleModels) for k, v in noe_dict.items()}

    def handle_prochirals(self, shift_pairs):
        """Group prochiral protons and filter pairs"""
        prochiral_groups = defaultdict(list)

        # Group by parent heavy atom
        for pair in shift_pairs:
            a1_parent = (pair['atom1'][:3])  # (chain, resnum, resname)
            a2_parent = (pair['atom2'][:3])
            prochiral_groups[a1_parent].append(pair['shift1'])
            prochiral_groups[a2_parent].append(pair['shift2'])

        # Filter groups within tolerance
        filtered_pairs = []
        for pair in shift_pairs:
            a1_parent = (pair['atom1'][:3])
            a2_parent = (pair['atom2'][:3])

            # Check if both protons are in prochiral groups
            if (max(prochiral_groups[a1_parent]) - min(prochiral_groups[a1_parent]) < self.prochiral_tolerance and
                    max(prochiral_groups[a2_parent]) - min(prochiral_groups[a2_parent]) < self.prochiral_tolerance):
                continue  # Exclude prochiral pairs

            filtered_pairs.append(pair)

        return filtered_pairs

    def generate_matched_pkl(self, peaklist_pid):

        pkl = self.project.getByPid(peaklist_pid)
        matched_peaks = self.matched_peaks.get(peaklist_pid, [])

        spectrum = pkl.spectrum

        matched_list = spectrum.newPeakList(
                title=f"RPF_Matched-{pkl.pid}",
                isSynthetic=True
                )

        for i, pk in enumerate(matched_peaks):

            peak = {
                'position_F1': pk['peak'].position_F1,
                'position_F2': pk['peak'].position_F2,
                'position_F3': pk['peak'].position_F3,
                'volume'     : pk['peak'].volume,
                'comment'    : f'Matched: {pk["best_match"]["atom1"]}-{pk["best_match"]["atom2"]}-{pk["best_match"]["heteroKey"]}',
                }

            matched_list.newPeak(
                    ppmPositions=[peak[f'position_F{dim}'] for dim in range(1, spectrum.dimensionCount + 1)],
                    comment=peak['comment'], volume=pk['peak'].volume, height=pk['peak'].height
                    )


    def generate_synthetic_peaks(self, peakList_pid):
        """Create synthetic peaks for missing pairs with validation"""
        synth_peaks = []
        peakListObject = self.project.getByPid(peakList_pid)
        if not peakListObject:
            raise ValueError(f"Peak list with PID {peakList_pid} not found")

        spectrum = peakListObject.spectrum

        try:

            missing_pairs = self.missing_pairs.get(peakListObject.pid, [])

            print(f"🔄 Generating synthetic peaks for {peakList_pid} ({len(missing_pairs)} missing pairs)")

            for i, pair in enumerate(missing_pairs):
                if not isinstance(pair, dict):
                    print(f"⚠️ Invalid pair format at index {i}: {type(pair)}")
                    continue

                required_keys = {'shift1', 'shift2', 'heteroShift', 'atom1', 'atom2', 'heteroKey'}
                if not required_keys.issubset(pair.keys()):
                    print(f"⚠️ Missing keys in pair {i}: {pair.keys()}")
                    continue

                peak = {
                    'position_F1': pair['shift1'],
                    'position_F2': pair['shift2'],
                    'position_F3': pair['heteroShift'],
                    'volume'     : 0,
                    'comment'    : f'Synthetic: {pair["atom1"]}-{pair["atom2"]}-{pair["heteroKey"]}',
                    }

                synth_peaks.append(peak)

            print(f"✅ Generated {len(synth_peaks)} synthetic peaks")
            return pd.DataFrame(synth_peaks) if synth_peaks else pd.DataFrame()

        except Exception as e:
            print(f"❌ Error in generate_synthetic_peaks: {str(e)}")
            import traceback

            traceback.print_exc()
            return pd.DataFrame()

    def calculate_per_residue_scores(self, matched, shift_pairs):
        """Per-residue RPF analysis"""

        residue_scores = defaultdict(lambda: {'tp': 0, 'fp': 0, 'fn': 0})

        # Count true positives
        for match in matched:
            res1 = match['best_match']['atom1'][:3]
            res2 = match['best_match']['atom2'][:3]
            residue_scores[res1]['tp'] += 1
            residue_scores[res2]['tp'] += 1

        # Count false negatives
        for pair in shift_pairs:
            res1 = pair['atom1'][:3]
            res2 = pair['atom2'][:3]
            if not any(m['best_match'] == pair for m in matched):
                residue_scores[res1]['fn'] += 1
                residue_scores[res2]['fn'] += 1

        # Calculate scores per residue
        scores = {}
        for res, counts in residue_scores.items():
            tp = counts['tp']
            fn = counts['fn']
            fp = counts['fp']

            recall = tp / (tp + fn) if (tp + fn) else 0
            precision = tp / (tp + fp) if (tp + fp) else 0
            f_measure = 2 * (precision * recall) / (precision + recall) if (precision + recall) else 0

            scores[res] = {
                'recall'   : recall,
                'precision': precision,
                'f_measure': f_measure
                }

        return scores

    def calc_dp_score(self, f_measure, true_pos, false_neg, true_pos_noe, false_pos_noe):
        """Discrimination potential score calculation"""
        # Random coil expectations using free chain values (index 2)
        try:
            rc_recall = true_pos[2] / (true_pos[2] + false_neg[2]) if (true_pos[2] + false_neg[2]) > 0 else 0
            rc_precision = true_pos_noe[2] / (true_pos_noe[2] + false_pos_noe[2]) if (true_pos_noe[2] + false_pos_noe[
                2]) > 0 else 0
            rc_f = 2 * (rc_recall * rc_precision) / (rc_recall + rc_precision) if (rc_recall + rc_precision) > 0 else 0

            return (f_measure - rc_f) / (1 - rc_f) if rc_f < 1 else 0.0
        except ZeroDivisionError:
            return 0.0

    def random_chain_noe(self, bond_count):
        """Calculate random chain NOE expectation"""
        if bond_count == 0:
            return 0.0
        r2 = bond_count * (1.5**2) * ((1.333) - (0.666 * (1 - (0.333**bond_count))))
        return r2**-3 if r2 > 0 else 0.0

    def calculate_rpf_metrics(self, matched_peaks, unexplained_peaks, missing_pairs):
        """Calculate metrics with enhanced validation"""
        try:
            # Validate input types
            if not isinstance(matched_peaks, list):
                print(f"⚠️ matched_peaks is {type(matched_peaks)}, expected list")
                matched_peaks = []

            if not isinstance(unexplained_peaks, list):
                print(f"⚠️ unexplained_peaks is {type(unexplained_peaks)}, expected list")
                unexplained_peaks = []

            # Empty state handling
            if not matched_peaks and not unexplained_peaks:
                print("⚠️ No peaks found in calculate_rpf_metrics")
                return {
                    'dp_score'       : 0.0,
                    'recall'         : 0.0,
                    'precision'      : 0.0,
                    'f_measure'      : 0.0,
                    'true_positives' : 0,
                    'false_positives': 0,
                    'false_negatives': len(missing_pairs)
                    }
            # Initialize counters as lists [full, local, free]
            true_pos = [0, 0, 0]
            false_neg = [0, 0, 0]
            true_pos_noe = [0.0, 0.0, 0.0]
            false_pos_noe = [0.0, 0.0, 0.0]

            # Calculate metrics for each category
            for match in matched_peaks:
                # Full count
                true_pos[0] += 1
                try:
                    true_pos_noe[0] += match.get('distance', 0)**self.NOE_POWER
                except:
                    true_pos_noe[0]

                # Local count (1-3 bonds)
                if 1 < match.get('bond_count', 0) < 4:
                    true_pos[1] += 1
                    try:
                        true_pos_noe[1] += match.get('distance', 0)**self.NOE_POWER
                    except:true_pos_noe[1]

                # Free chain count
                true_pos[2] += 1
                true_pos_noe[2] += self.random_chain_noe(match.get('bond_count', 0))

            # Calculate false negatives (missing pairs)
            false_neg[0] = len(missing_pairs)
            # Add logic to calculate local/free chain false negatives if available

            # Calculate base metrics
            tp = true_pos[0]  # True positives (full count)
            fp = sum(len(m.get('matches', [])) - 1 for m in matched_peaks)  # False positives
            fn = false_neg[0]  # False negatives

            # Safeguard calculations
            precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
            f_measure = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0

            # Calculate DP score using original formula
            dp_score = self.calc_dp_score(
                    f_measure,
                    true_pos,  # Pass list [full, local, free]
                    false_neg,  # Pass list [full, local, free]
                    true_pos_noe,  # Pass list [full, local, free]
                    false_pos_noe  # Pass list [full, local, free]
                    )

            return {
                'dp_score'       : round(dp_score, 3),
                'recall'         : round(recall, 3),
                'precision'      : round(precision, 3),
                'f_measure'      : round(f_measure, 3),
                'true_positives' : tp,
                'false_positives': fp,
                'false_negatives': fn
                }

        except Exception as e:
            print(f"❌ Critical error in calculate_rpf_metrics: {str(e)}")
            import traceback

            traceback.print_exc()
            return {}

    def _pdbName_to_possible_NEF_names(self, name, residue_type):
        """Translate PDB atom names to NEF with residue-specific rules"""
        # Define residue type groups from your code


        original_name = name
        expanded = []

        if original_name == 'HN':
            expanded = ['H', 'HN']
        elif original_name == 'HA':
            expanded = ['HA']
        elif original_name == 'HA2' and residue_type == 'GLY':
            expanded = ['HA2', 'HAx', 'HA%']
        elif original_name == 'HA3' and residue_type == 'GLY':
            expanded = ['HA3', 'HAy', 'HA%']

        elif original_name == 'HB1' and residue_type in ['ALA']:
            expanded = ['HB1', 'HB%']
        elif original_name == 'HB2' and residue_type in ['ALA']:
            expanded = ['HB2', 'HB%']
        elif original_name == 'HB3' and residue_type in ['ALA']:
            expanded = ['HB3', 'HB%']

        elif original_name == 'HB2' and residue_type in ['SER', 'LEU', 'TRP']:
            expanded = ['HB2', 'HBx', 'HB%']
        elif original_name == 'HB3' and residue_type in ['SER', 'LEU', 'TRP']:
            expanded = ['HB3', 'HBy', 'HB%']

        elif original_name == 'HG12' and residue_type in ['ILE']:
            expanded = ['HG12', 'HG1x', 'HG1%']
        elif original_name == 'HG13' and residue_type in ['ILE']:
            expanded = ['HG13', 'HG1y', 'HG1%']

        elif original_name == 'HG21' and residue_type in ['THR']:
            expanded = ['HG21', 'HG2%']
        elif original_name == 'HG22' and residue_type in ['THR']:
            expanded = ['HG22', 'HG2%']
        elif original_name == 'HG23' and residue_type in ['THR']:
            expanded = ['HG23', 'HG2%']

        elif original_name == 'HG11' and residue_type in ['VAL']:
            expanded = ['HG11', 'HG1%', 'HGx%', 'HG%']
        elif original_name == 'HG12' and residue_type in ['VAL']:
            expanded = ['HG12', 'HG1%', 'HGx%', 'HG%']
        elif original_name == 'HG13' and residue_type in ['VAL']:
            expanded = ['HG13', 'HG1%', 'HGx%', 'HG%']
        elif original_name == 'HG21' and residue_type in ['VAL', 'ILE']:
            expanded = ['HG21', 'HG2%', 'HGy%', 'HG%']
        elif original_name == 'HG22' and residue_type in ['VAL', 'ILE']:
            expanded = ['HG22', 'HG2%', 'HGy%', 'HG%']
        elif original_name == 'HG23' and residue_type in ['VAL', 'ILE']:
            expanded = ['HG23', 'HG2%', 'HGy%', 'HG%']

        elif original_name == 'HD11' and residue_type in ['LEU', 'ILE']:
            expanded = ['HD11', 'HD1%', 'HDx%', 'HDG%']
        elif original_name == 'HD12' and residue_type in ['LEU', 'ILE']:
            expanded = ['HG12', 'HD1%', 'HDx%', 'HD%']
        elif original_name == 'HD13' and residue_type in ['LEU', 'ILE']:
            expanded = ['HD13', 'HD1%', 'HDx%', 'HD%']
        elif original_name == 'HD21' and residue_type in ['LEU']:
            expanded = ['HD21', 'HD2%', 'HDy%', 'HD%']
        elif original_name == 'HD22' and residue_type in ['LEU']:
            expanded = ['HD22', 'HD2%', 'HDy%', 'HD%']
        elif original_name == 'HD23' and residue_type in ['LEU']:
            expanded = ['HD23', 'HD2%', 'HDy%', 'HD%']

        elif original_name == 'HB2' and residue_type in ['MET', 'HIS', 'PHE', 'PRO', 'ASP', 'ASN', 'GLU', 'GLN','ARG','CYS']:
            expanded = ['HB2', 'HBx', 'HB%']
        elif original_name == 'HB3' and residue_type in ['MET', 'HIS', 'PHE', 'PRO','ASP', 'ASN', 'GLU', 'GLN','ARG','CYS']:
            expanded = ['HB3', 'HBy', 'HB%']

        elif original_name == 'HG2' and residue_type in ['MET','PRO', 'GLN', 'GLU','LYS','ARG']:
            expanded = ['HG2', 'HGx', 'HG%']
        elif original_name == 'HG3' and residue_type in ['MET','PRO', 'GLN', 'GLU','LYS','ARG']:
            expanded = ['HG3', 'HGy', 'HG%']
        elif original_name == 'HD2' and residue_type in ['PRO','LYS','ARG']:
            expanded = ['HD2', 'HDx', 'HD%']
        elif original_name == 'HD3' and residue_type in ['PRO','LYS','ARG']:
            expanded = ['HD3', 'HDy', 'HD%']
        elif original_name == 'HE2' and residue_type in ['LYS']:
            expanded = ['HE2', 'HEx', 'HE%']
        elif original_name == 'HE3' and residue_type in ['LYS']:
            expanded = ['HE3', 'HEy', 'HE%']

        elif original_name == 'HD21' and residue_type in ['PRO', 'ASP']:
            expanded = ['HD21', 'HD2x', 'HD%']
        elif original_name == 'HD22' and residue_type in ['PRO','ASP']:
            expanded = ['HD22', 'HD2y', 'HD%']

        elif original_name == 'HE21' and residue_type in ['GLN']:
            expanded = ['HE21', 'HE2x', 'HE%']
        elif original_name == 'HE22' and residue_type in ['GLN']:
            expanded = ['HE22', 'HE2y', 'HE%']

        elif original_name == 'HZ1' and residue_type in ['LYS']:
            expanded = ['HZ1', 'HZ%']
        elif original_name == 'HZ2' and residue_type in ['LYS']:
            expanded = ['HZ2', 'HZ%']
        elif original_name == 'HZ3' and residue_type in ['LYS']:
            expanded = ['HZ3', 'HZ%']

        elif original_name == 'HH11' and residue_type in ['ARG']:
            expanded = ['HH11', 'HH1x', 'HH1%']
        elif original_name == 'HH12' and residue_type in ['ARG']:
            expanded = ['HH12', 'HH2y', 'HH2%']

        elif original_name == 'HD1' and residue_type in ['PHE', 'TYR']:
            expanded = ['HD1', 'HDx', 'HD%']
        elif original_name == 'HD2' and residue_type in ['PHE', 'TYR']:
            expanded = ['HD2', 'HDy', 'HD%']
        elif original_name == 'HE1' and residue_type in ['PHE', 'TYR']:
            expanded = ['HE1', 'HEx', 'HE%']
        elif original_name == 'HE2' and residue_type in ['PHE', 'TYR']:
            expanded = ['HE2', 'HEy', 'HE%']

        else: expanded = [original_name]

        return expanded

    def _get_possible_shifts(self, atom_key, shift_map):
        """Expand possible shifts considering atom name ambiguity"""
        chain, resnum, resname, atom = atom_key
        possible_atoms = self._pdbName_to_possible_NEF_names(atom, residue_type=resname.upper())

        shifts = []

        for variant in possible_atoms:
            key = (chain, resnum, resname, variant)
            if key in shift_map:
                shifts.append([key, shift_map[key]])

        return shifts if shifts else [[None,None]]

    def get_unique_pairs(self, pairs_list):
        """Filter pairs to ensure uniqueness (A-B == B-A)"""
        seen = set()
        unique_pairs = []

        for pair in pairs_list:
            # Extract atoms directly from the tuple
            a1 = pair[0]  # First atom in the pair
            a2 = pair[1]  # Second atom in the pair

            # Get basic info from atom objects
            a1_info = (
                a1.parent.parent.id,  # Chain
                a1.parent.id[1],  # Residue number
                a1.name  # Atom name
                )

            a2_info = (
                a2.parent.parent.id,
                a2.parent.id[1],
                a2.name
                )

            key = frozenset({a1, a2})  # Order-agnostic

            if key not in seen:
                seen.add(key)
                unique_pairs.append(pair)

        return unique_pairs

    def get_valid_context_pairs(self, protons, contexts):
        """With pair debugging"""
        valid_pairs = []
        model = self.ensembleModels[0]

        ns = NeighborSearch(list(model.get_atoms()))
        print(f"  Atom count in model: {len(list(model.get_atoms()))}")

        for i, (h1, h2) in enumerate(combinations(protons, 2)):
            try:
                h1_neighbors = ns.search(h1.coord, 1.2, level='A')
                h1_neighbors = [a for a in h1_neighbors if a != h1]
                h2_neighbors = ns.search(h2.coord, 1.2, level='A')
                h2_neighbors = [a for a in h2_neighbors if a != h2]
                if not h1_neighbors or not h2_neighbors:
                    continue

                e1 = h1_neighbors[0].element
                e2 = h2_neighbors[0].element

                if (e1, e2) in contexts or (e2, e1) in contexts:
                    valid_pairs.append((h1, h2))
            except Exception as e:
                print(f"⚠️ Error with pair {i}: {str(e)}")
                continue

        return valid_pairs

    def map_to_nmr_atoms(self, project: Project, pairs_list):
        """Convert pairs to CCPN NmrAtom pairs with ambiguity handling"""
        # Precompute all NmrAtom PIDs
        pid_map = {atom.pid: atom for atom in project.nmrAtoms}

        result = []

        for pair in pairs_list:
            # Extract atoms directly from the tuple
            a1 = pair[0]  # First atom in the pair
            a2 = pair[1]  # Second atom in the pair

            # Get basic info from atom objects
            a1_info = (
                a1.parent.parent.id,  # Chain
                a1.parent.id[1],  # Residue number
                a1.name  # Atom name
                )

            a2_info = (
                a2.parent.parent.id,
                a2.parent.id[1],
                a2.name
                )

            # Generate possible PIDs for each atom
            def get_possible_pids(chain, resnum, atom_name):
                pids = []
                for variant in self._pdbName_to_possible_NEF_names(atom_name):
                    # CCPN PID format: NA:<chain>.<resnum>.<resname>.<atom>
                    resName = self.get_residueType_from_number(int(resnum))
                    pid_pattern = f"NA:{chain}.{resnum}.{resName}.{variant}"
                    compiled = re.compile(pid_pattern.replace('.', r'\.').replace('*', '.*'))
                    pids.extend([pid for pid in pid_map if compiled.match(pid)])
                return [pid_map[pid] for pid in pids if pid in pid_map]

            # Find possible matches for both atoms
            a1_matches = get_possible_pids(*a1_info)
            a2_matches = get_possible_pids(*a2_info)

            # Create all valid combinations
            for nmr_atom1, nmr_atom2 in product(a1_matches, a2_matches):
                # Avoid self-pairs and ensure chain consistency
                if nmr_atom1 != nmr_atom2:
                    result.append((nmr_atom1, nmr_atom2))

        return list(set(result))  # Remove duplicates

    def load_pdb_ensemble(self):
        """Load PDB file with potential multiple models"""
        parser = PDBParser(QUIET=True)
        return parser.get_structure("ensemble", self.pdb_path)

    def get_bound_protons(self, search_element = 'H', target_elements: list =[], bond_cutoff: float = 1.2):
        """Find protons bound to specified elements using spatial proximity"""
        model = self.ensembleModels[0]
        ns = NeighborSearch(list(model.get_atoms()))
        protons = []

        for atom in model.get_atoms():
            if atom.element == search_element:
                # Find nearest non-hydrogen atom within bond distance
                neighbors = ns.search(atom.coord, bond_cutoff, level='A')
                for neighbor in neighbors:
                    if neighbor.element in target_elements and neighbor.element != 'H':
                        protons.append(atom)
                        break
        return protons

    def calculate_min_distance(self, pair, models):
        """Calculate minimum distance with proper error handling"""
        h1, h2 = pair
        min_dist = float('inf')

        for model in models:
            try:
                # Get coordinates using full hierarchy
                m_h1 = model[h1.parent.parent.id][h1.parent.id[1]][h1.name]
                m_h2 = model[h2.parent.parent.id][h2.parent.id[1]][h2.name]  # Same residue assumption

                dist = np.linalg.norm(m_h1.coord - m_h2.coord)
                if dist < min_dist:
                    min_dist = dist

            except KeyError as e:
                print(
                        f"⚠️ Missing atom in model {model.id}: {e} {h1.parent.parent.id}{h1.parent.id[1]}{h1.name}{h2.parent.parent.id}{h2.parent.id[1]}{h2.name}")
                continue
            except Exception as e:
                print(f"❌ Unexpected error in model {model.id}: {str(e)}")
                continue

        return min_dist if min_dist != float('inf') else None

    def get_proton_pairs(self, distance_threshold: float, contexts: tuple, target_element: str):
        """Main function with debug prints"""
        print(f"⏳ Loading PDB ensemble from {self.pdb_path}")
        ensemble = self.ensemble
        models = list(ensemble.get_models())
        print(f"✅ Loaded {len(models)} models")

        if not models:
            print("❌ Error: No models found in PDB file")
            return []

        print(f"🔍 Finding protons bound to {[target_element]} in first model...")
        protons = self.get_bound_protons(target_elements=[target_element])
        print(f"🔵 Found {len(protons)} potential protons")

        if not protons:
            print("❌ No protons found - check your PDB's hydrogen atoms")
            return []

        print(f"🔄 Generating valid context pairs from {len(protons)} protons...")
        valid_pairs = self.get_valid_context_pairs(protons, contexts)
        print(f"🔷 Found {len(valid_pairs)} context-valid pairs")

        results = []
        print(f"📏 Calculating distances across {len(models)} models...")
        for i, pair in enumerate(valid_pairs):
            min_dist = self.calculate_min_distance(pair, models)
            if min_dist <= distance_threshold:
                results.append({
                    'atom1'   : pair[0],
                    'atom2'   : pair[1],
                    'distance': min_dist
                    })

        results = self.get_unique_pairs(results)

        print(f"🎉 Final results: {len(results)} pairs under {distance_threshold}Å")
        return results

    def create_shift_map(self):
        """Create a dictionary mapping atoms to their chemical shifts"""
        chemicalShift_list = self.project.getByPid(self.chemical_shift_pid)
        chemicalShift_list_df = chemicalShift_list.getAsDataFrame()

        shift_map = defaultdict(dict)

        for _, row in chemicalShift_list_df.iterrows():
            key = (
                row['chainCode'],
                row['sequenceCode'],
                row['residueType'],
                row['atomName'])
            shift_map[key] = row['value']

        self.shift_map = shift_map

    def get_unmatched_peaks_report(self, unexplained_peaks, missing_pairs):
        """Generate detailed report of unmatched peaks and missing pairs"""
        report = {
            'unexplained_peaks': [],
            'missing_pairs'    : []
            }

        # Unexplained peaks
        for peak in unexplained_peaks:
            report['unexplained_peaks'].append({
                'peak_id'    : peak['serial'],
                'position_F1': peak['position_F1'],
                'position_F2': peak['position_F2'],
                'volume'     : peak['volume'],
                'snr'        : peak['signalToNoiseRatio']
                })

        # Missing structural pairs
        for pair in missing_pairs:
            report['missing_pairs'].append({
                'atom1'          : f"{pair['atom1'][0]}{pair['atom1'][1]}{pair['atom1'][3]}",
                'atom2'          : f"{pair['atom2'][0]}{pair['atom2'][1]}{pair['atom2'][3]}",
                'distance'       : pair['distance'],
                'expected_shift1': pair['shift1'],
                'expected_shift2': pair['shift2']
                })

        return report

    def get_structural_shift_pairs(self, proton_dists_conn, spectrum, tolerance=0.04):
        """Convert structural pairs to shift-based pairs with ambiguity handling"""
        shift_pairs = []

        for (a1, a2), distance in proton_dists_conn.items():
            # for each atom get the atom it is attached to (heteroAtom1) its
            # chemical shift, a1_het, and its element (CNOS)
            a1_het, element1, heteroAtom1 = self.get_chemicalShift_of_attached_atom(a1)
            a2_het, element2, heteroAtom2 = self.get_chemicalShift_of_attached_atom(a2)

            # Get atom identifiers and build keys
            key1 = (
                a1.parent.parent.id,
                a1.parent.id[1],
                a1.parent.resname,
                a1.name
                )

            key2 = (
                a2.parent.parent.id,
                a2.parent.id[1],
                a2.parent.resname,
                a2.name
                )

            # possible_shifts must be filtered for N dimentyion
            results = []

            oneBondHeteroAtom  = 'N'

            #O go through dimensions - find the nonH: ('Hn', 'H', 'Nh')
            for e in spectrum.referenceExperimentDimensions:
                if len(e) != 1 and e[0] != 'H':
                    oneBondHeteroAtom = e[0]
                    break

            if element1 == oneBondHeteroAtom:
                # a1 is Hn/Hc
                possible_shifts1 = self._get_possible_shifts(key1, self.shift_map)
                possible_shifts2 = self._get_possible_shifts(key2, self.shift_map)

                heteroShift = a1_het

                heteroKey = (
                    heteroAtom1.parent.parent.id,
                    heteroAtom1.parent.id[1],
                    heteroAtom1.parent.resname,
                    heteroAtom1.name
                    )
                heteroKey2 = (
                    heteroAtom2.parent.parent.id,
                    heteroAtom2.parent.id[1],
                    heteroAtom2.parent.resname,
                    heteroAtom2.name
                    )

                for shift1 in possible_shifts1:
                    for shift2 in possible_shifts2:
                        results.append([shift1[1], shift2[1],
                                        heteroShift, key1, key2, heteroKey, a2_het, heteroKey2])

            elif element2 == oneBondHeteroAtom:
                # a2 is Hn
                possible_shifts1 = self._get_possible_shifts(key2, self.shift_map)
                possible_shifts2 = self._get_possible_shifts(key1, self.shift_map)
                heteroShift = a2_het

                heteroKey = (
                    heteroAtom2.parent.parent.id,
                    heteroAtom2.parent.id[1],
                    heteroAtom2.parent.resname,
                    heteroAtom2.name
                    )
                heteroKey2 = (
                    heteroAtom1.parent.parent.id,
                    heteroAtom1.parent.id[1],
                    heteroAtom1.parent.resname,
                    heteroAtom1.name
                    )

                for shift1 in possible_shifts1:
                    for shift2 in possible_shifts2:
                        results.append([shift1[1], shift2[1],
                                        heteroShift, key2, key1, heteroKey, a1_het, heteroKey2])

            for result in results:
                shift1 = result[0]
                shift2 = result[1]
                heteroShift = result[2]
                key1 = result[3]
                key2 = result[4]
                heteroKey = result[5]
                heteroShift2 = result[6]
                heteroKey2 = result[7]

                if shift1 is not None and shift2 is not None and abs(shift1 - shift2) > tolerance:  # Avoid diagonal
                    shift_pairs.append({
                        'shift1'  : shift1,
                        'shift2'  : shift2,
                        'heteroShift': heteroShift,
                        'distance': distance,
                        'atom1'   : key1,
                        'atom2'   : key2,
                        'heteroKey': heteroKey,
                        'heteroShift2': heteroShift2,#might need in future for 4d- not currently used
                        'heteroKey2': heteroKey2 #might need in future for 4d - not currently used
                        })

        return shift_pairs

    def match_peaks_to_structure(self, shift_pairs, peak_df):
        """Match peaks to structural pairs using chemical shifts
        shift_pairs: List of structural pairs with chemical shifts
        (expected format: [{'shift1': float, 'shift2': float, ...}, ...])

        peak_df: DataFrame of experimental peaks (requires columns
        position_F1 and position_F2 and position F3)

        """
        matched_peaks = []
        unexplained_peaks = []

        """
        Creates a spatial index (KDTree) of structural pair shifts for fast lookup

        For each experimental peak:
        Finds structural pairs within prochiral_tolerance distance
        Records matches and identifies the closest structural pair
        Tracks peaks with no matches
        """

        # Deduplicate shift_pairs - this needed because shift pairs are duplicated for
        # ambiguos atoms eg HG1, HG2, HG3 would give 3
        unique_shift_pairs = []
        seen_shifts = set()

        for sp in shift_pairs:
            try:
                # Create a unique key (e.g., tuple of shifts)
                shift_key = (float(sp['shift1']), float(sp['shift2']), float(sp['heteroShift']))
                if shift_key not in seen_shifts:
                    seen_shifts.add(shift_key)
                    unique_shift_pairs.append(sp)
            except (TypeError, ValueError, KeyError):
                continue  # Skip invalid entries

        shift_pairs = unique_shift_pairs  # Use deduplicated list

        filtered_pairs = []
        original_indices = []  # Maps expanded entries to original shift_pairs indices
        for idx, pair in enumerate(shift_pairs):
            # check if its possible for the shiftpair to be assignable (all dimensions have values)
            try:
                shift1 = float(pair['shift1'])
                shift2 = float(pair['shift2'])
                heteroShift = float(pair['heteroShift'])
                filtered_pairs.append((shift1, shift2, heteroShift))
                original_indices.append(idx)
            except (TypeError, ValueError):
                pass  # Skip if any value can't be converted to float

        if not filtered_pairs:
            return [], peak_df.to_dict('records'), []

        # Build KDTree from expanded pairs
        tree = KDTree(np.array(filtered_pairs))

        for _, peak in peak_df.iterrows():
            if any(pd.isna(peak.get(key)) for key in ['position_F1', 'position_F2', 'position_F3']):
                print(f"Unexplained peak (NaN detected): {peak}")
                unexplained_peaks.append(peak)
                continue

            # Query expanded pairs
            query_point = [peak['position_F1'], peak['position_F2'], peak['position_F3']]
            indices = tree.query_ball_point(query_point, self.prochiral_tolerance)

            if indices:
                # Map expanded indices to original indices
                matched_orig_idxs = set(original_indices[idx] for idx in indices)
                matches = [shift_pairs[i] for i in matched_orig_idxs]
                best_match = min(matches, key=lambda x: x['distance'])
                matched_peaks.append({
                    'peak'      : peak,
                    'matches'   : matches,
                    'best_match': best_match,
                    'indices'   : list(matched_orig_idxs)
                    })
            else:
                unexplained_peaks.append(peak)

        # Identify missing pairs (original indices not matched)
        all_matched = set()
        for mp in matched_peaks:
            all_matched.update(mp['indices'])

        missing_pairs = []
        print("All_matched:", len(all_matched))

        for idx, sp in enumerate(shift_pairs):
            if idx in all_matched:
                continue  # Skip matched entries
            try:
                shift1 = float(sp['shift1'])
                shift2 = float(sp['shift2'])
                heteroShift = float(sp['heteroShift'])
                missing_pairs.append(sp)  # Only append if all conversions succeed
            except (TypeError, ValueError, KeyError):
                continue  # Skip invalid entries

        return matched_peaks, unexplained_peaks, missing_pairs

    def get_spectrum_context(self, spectrum) -> tuple:
        """Determine heteroatom context with validation"""
        try:
            if not spectrum:
                print("⚠️ No spectrum provided")
                return ()

            carbon_dims = 0
            nitrogen_dims = 0

            print(f"🔍 Analyzing spectrum dimensions: {spectrum.dimensionTriples}")

            for dim in spectrum.dimensionTriples:
                if len(dim) < 2:
                    print(f"⚠️ Invalid dimension tuple: {dim}")
                    continue

                code = dim[1].upper()  # Normalize to uppercase

                # Improved element detection
                if code.startswith('C') and 'N' not in code:
                    carbon_dims += 1
                elif code.startswith('N') and 'C' not in code:
                    nitrogen_dims += 1

            print(f"🔬 Detected C-dims: {carbon_dims}, N-dims: {nitrogen_dims}")

            context = self.CN_DIM_CONTEXT_DICT.get((carbon_dims, nitrogen_dims), ())
            print(f"📦 Context: {context}")
            return context

        except Exception as e:
            print(f"❌ Error in get_spectrum_context: {str(e)}")
            return ()

    def merge_rpf_results(self, matched_dict, unexplained_dict, missing_dict):
        """Merge RPF results from multiple peak lists, removing duplicates

        Args:
            matched_dict (dict): {peaklist_pid: [matched_entries]}
            unexplained_dict (dict): {peaklist_pid: [unexplained_peaks]}
            missing_dict (dict): {peaklist_pid: [missing_pairs]}

        Returns:
            tuple: (merged_matched, merged_unexplained, merged_missing)
        """
        merged_matched = []
        merged_unexplained = []
        merged_missing = []

        # Track unique identifiers
        seen_peaks = set()  # (spectrum_id, peak_id)
        seen_pairs = set()  # (atom1_id, atom2_id)

        # Helper to get peak signature
        def _get_peak_id(peak):
            return

        # Helper to get pair signature
        def _get_pair_id(pair):
            return tuple(sorted([
                (pair['atom1']['chain'], pair['atom1']['resnum'], pair['atom1']['name']),
                (pair['atom2']['chain'], pair['atom2']['resnum'], pair['atom2']['name'])
                ]))

        # Process matched peaks
        for pl_pid in matched_dict:
            for m in matched_dict[pl_pid]:
                peak_id = m['peak']['pid']
                if peak_id not in seen_peaks:
                    merged_matched.append(m)
                    seen_peaks.add(peak_id)

        # Process unexplained peaks
        for pl_pid in unexplained_dict:
            for peak in unexplained_dict[pl_pid]:
                peak_id = peak['pid']
                if peak_id not in seen_peaks:
                    merged_unexplained.append(peak)
                    seen_peaks.add(peak_id)

        # Process missing pairs
        for pl_pid in missing_dict.keys():
            merged_missing = merged_missing + missing_dict[pl_pid]

        return merged_matched, merged_unexplained, merged_missing

    def get_chain(self):
        csl = self.project.getByPid(self.chemical_shift_pid)
        cs = csl.chemicalShifts[0]
        chain = cs.nmrAtom.nmrResidue.nmrChain.chain
        return chain

    def get_residueType_from_number(self, resiNum):
        resiName = None
        chain = self.get_chain()
        for resi in chain:
            num = int(resi.sequenceCode)
            if num == int(resiNum):
                resiName = resi.residueType
                break
        return resiName

    def create_synthetic_peaklist(self, peaklist_pid):
        """Create CCPN PeakList for missing peaks"""

        pkl = self.project.getByPid(peaklist_pid)

        spectrum = pkl.spectrum

        synthetic_list = spectrum.newPeakList(
                title=f"RPF_Missing-{pkl.pid}",
                isSynthetic=True
                )
        # Add peaks from syntheticPeaks DataFrame
        for _, row in self.syntheticPeaks[peaklist_pid].iterrows():
            synthetic_list.newPeak(
                    ppmPositions=[row[f'position_F{dim}'] for dim in range(1, spectrum.dimensionCount + 1)],
                    comment=row['comment']
                    )


    def run_rpf_analysis(self):
        # Load data
        self.peaklists = self.getPeakList_objects()
        self.ensemble = self.load_pdb_ensemble()
        self.ensembleModels = [model for model in self.ensemble]
        self.chain = self.get_chain()

        # Get chemical shifts as self.shift_map
        self.create_shift_map()

        for pkl in self.peaklists:
            peak_df = pkl.getAsDataFrame()
            # Get structural pairs with NOE weighting
            spectrum = pkl.spectrum
            contexts = self.get_spectrum_context(spectrum)
            proton_dists_conn = self.get_proton_dists_conn(contexts, self.distance_threshold)
            self.proton_dists_conn[pkl.pid] = proton_dists_conn

            # Convert to shift pairs with prochiral exclusion
            shift_pairs = self.get_structural_shift_pairs(proton_dists_conn,
                                                          spectrum= spectrum,
                                                          tolerance=0.04)

            #shift_pairs.append({
            #    'shift1'  : shift1,
            #    'shift2'  : shift2,
            #    'distance': distance,
            #    'atom1'   : key1,
            #    'atom2'   : key2
            #    })

            shift_pairs = self.handle_prochirals(shift_pairs)

            self.shift_pairs[pkl.pid] = shift_pairs

            # Get tolerance dictionary from spectrum
            tol_dict = {dim: (0.05, 0.05) for dim in spectrum.dimensions}  #
            # tol_dict = {1: (0.05, 0.05), 2: (0.05, 0.05), 3: (0.05, 0.05)}


            # Perform ambiguous matching
            # i'm here - i think that i need to do somehting with the peak_assingments - possibly combine with matched
            # check v2
            connections, peak_assignments = self.get_ambig_noe_conn(self.shift_pairs[pkl.pid],
                                                                    peak_df,
                                                                    tol_dict)

            """
            matched_peaks: List of dictionaries for matched peaks containing:

            {
              'peak':       Peak data (pd.Series),
              'matches':    All matching structural pairs (list),
              'best_match': Closest structural pair (dict),
              'indices':    Indices of matches in shift_pairs (list)
            }

            unexplained_peaks: List of peaks with no matching structural pairs

            The missing_pairs list contains structural pairs (from the 
            PDB-derived shift_pairs) that:
                Have valid chemical shifts (i.e., they are within the expected 
                spectral region) but which Are NOT matched to any experimental peaks 
                in the peak_df list 
                (even after considering the prochiral_tolerance).
            """
            matched_peaks, unexplained_peaks, missing_pairs = self.match_peaks_to_structure(self.shift_pairs[pkl.pid],
                                                                                            peak_df)
            self.matched_peaks[pkl.pid] = matched_peaks
            self.unexplained_peaks[pkl.pid] = unexplained_peaks
            self.missing_pairs[pkl.pid] = missing_pairs
            self.syntheticPeaks[pkl.pid] = self.generate_synthetic_peaks(pkl.pid)

            # Calculate per residue scores
            self.residue_scores[pkl.pid] = self.calculate_per_residue_scores(self.matched_peaks[pkl.pid],
                                                                             self.shift_pairs[pkl.pid])

            # Generate reports
            self.unmatched_peaks_report[pkl.pid] = self.get_unmatched_peaks_report(self.unexplained_peaks[pkl.pid],
                                                                                   self.missing_pairs[pkl.pid])
            self.metrics[pkl.pid] = self.calculate_rpf_metrics(self.matched_peaks[pkl.pid],
                                                               self.unexplained_peaks[pkl.pid],
                                                               self.missing_pairs[pkl.pid])

        # After processing all peak lists
        merged_matched, merged_unexplained, merged_missing = self.merge_rpf_results(
                self.matched_peaks,
                self.unexplained_peaks,
                self.missing_pairs
                )

        # Generate reports
        self.unmatched_peaks_report['overall'] = self.get_unmatched_peaks_report(merged_unexplained, merged_missing)
        self.metrics['overall'] = self.calculate_rpf_metrics(merged_matched,
                                                             merged_unexplained,
                                                             merged_missing)


from ccpn.ui.gui.widgets.PulldownListsForObjects import ChemicalShiftListPulldown
from ccpn.ui.gui.popups.Dialog import CcpnDialogMainWidget
from ccpn.ui.gui.widgets.ListWidget import ListWidgetPair
from ccpn.ui.gui.widgets.Button import Button
from ccpn.ui.gui.widgets.Label import Label
from ccpn.ui.gui.widgets import MessageDialog
from ccpn.ui.gui.widgets import Entry
from ccpn.ui.gui.widgets.FileDialog import OtherFileDialog
from ccpn.ui.gui.widgets.Frame import Frame
from ccpn.framework.Application import getApplication
from ccpn.util.Path import aPath

class RpfGui(CcpnDialogMainWidget):
    FIXEDWIDTH = True
    title = 'RPF Analysis'

    def __init__(self, parent=None, mainWindow=None, title=title, **kwds):
        super().__init__(parent, setLayout=True, windowTitle=title, **kwds)
        self.mainWindow = mainWindow
        self.application = getApplication()
        self.project = self.application.project if self.application else None

        self._createWidgets()

    def _createWidgets(self):
        row = 0
        # PDB File Input
        Label(self.mainWidget, text="PDB File:", grid=(row, 0))
        self.pdbLineEdit = Entry.Entry(self.mainWidget, grid=(row, 1))
        self.browseButton = Button(self.mainWidget, text="Browse", callback=self._browsePDB, grid=(row, 2))

        row += 1
        Label(self.mainWidget, text="Distance Tolerence:", grid=(row, 0))
        self.distTolerance = Entry.FloatEntry(self.mainWidget, grid=(row, 1), text=4.8)

        row += 1
        #self.cslLabel = Label(self.mainWidget, text='Chemical Shift List', grid=(row, 0), hAlign='right')
        self.chemicalShiftListWidget = ChemicalShiftListPulldown(self.mainWidget, grid=(row, 0), gridSpan=(1, 2))

        row += 1
        self.plsLabel = Label(self.mainWidget, text='PeakLists', grid=(row, 0), hAlign='right')
        self.peakListWidget = ListWidgetPair(self.mainWidget, grid=(row, 1), gridSpan=(1, 1))
        self.peakListWidget.setMaximumHeight(100)
        row += 1
        # Buttons
        self.ButtonFrame = Frame(parent=self.mainWidget, setLayout=True, grid=(row, 0), gridSpan=(1, 3))
        self.analyzeButton = Button(self.ButtonFrame, text="Analyze", callback=self._runAnalysis, grid=(0, 0))
        self.closeButton = Button(self.ButtonFrame, text="Close", callback=self.reject, grid=(0, 1))

        self._populateWidgets()

    def _populateWidgets(self):
        # default file browser to project path
        #self.pdbLineEdit.setText(self.project.path)

        # get noesy peakList's
        _peakLists = []
        for pkl in self.project.peakLists:
            transfers = [transfer.transferType for transfer in pkl.spectrum.magnetisationTransfers]
            if 'through-space' in transfers:
                _peakLists.append(pkl)

        self.peakListWidget.setListObjects(_peakLists)

    def _browsePDB(self):
        _currentPath = self.pdbLineEdit.getText()
        if _currentPath is not None:
            _directory = aPath(_currentPath).parent.asString()
        else:
            _directory = self.project.application.dataPath.asString()

        dialog = OtherFileDialog(parent=self.mainWindow, _useDirectoryOnly=False,
                                 directory=_directory, fileFilter="*.pdb")
        dialog._show()
        if (path := dialog.selectedFile()) is not None:
            self.pdbLineEdit.setText(str(path))


    def _runAnalysis(self):
        pdb_path = self.pdbLineEdit.text()
        if not pdb_path:
            MessageDialog.showWarning("No PDB File", "Please select a PDB file.")
            return

        chemical_shift_pid = self.chemicalShiftListWidget.getText()

        peak_list_pids = self.peakListWidget.getRightList()

        project = self.project

        distance_threshold = float(self.distTolerance.getText())

        # Instantiate the calculator
        self.rpf_calc = RpfCalculator(
                project=project,
                pdb_path=pdb_path,
                peak_list_pids=peak_list_pids,
                chemical_shift_pid=chemical_shift_pid,
                distance_threshold=distance_threshold,
                prochiral_tolerance=0.04
                )

        self.rpf_calc.run_rpf_analysis()

        df = pd.DataFrame.from_dict(rpf_calc.metrics['overall'], orient='index', columns=['value'])
        df.reset_index(inplace=True)
        df.columns = ['metric', 'value']

        self.project.newDataTable(name="RPF_Metrics_Overall", data=df)

        for peak_list_pid in peak_list_pids:
            self.rpf_calc.create_synthetic_peaklist(peaklist_pid=peak_list_pid)
            self.rpf_calc.generate_matched_pkl(peaklist_pid=peak_list_pid)

            df = pd.DataFrame.from_dict(self.rpf_calc.metrics['overall'], orient='index', columns=['value'])
            df.reset_index(inplace=True)
            df.columns = ['metric', 'value']
            self.project.newDataTable(name="RPF_Metrics_" + peak_list_pid, data=df)

            df = pd.DataFrame.from_dict(self.rpf_calc.residue_scores[peak_list_pid], orient='index')
            df.index = pd.MultiIndex.from_tuples(df.index, names=['chain', 'residue_number', 'residue_name'])
            df.reset_index(inplace=True)
            self.project.newDataTable(name="RPF_PerResidueScores_" + peak_list_pid, data=df)


if __name__ == '__main__':
    # from ccpn.ui.gui.widgets.Application import TestApplication
    #
    # app = TestApplication()
    popup = RpfGui()
    popup.exec_()