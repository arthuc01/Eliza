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
                 chemical_shift_pid: str, distance_threshold: float = 4.8, prochiral_tolerance: float =0.04,
                 h_tolerance: float =0.04, c_tolerance: float =0.3, n_tolerance: float =0.3):
        self.project = project
        self.pdb_path = pdb_path
        self.peak_list_pids = peak_list_pids
        self.chemical_shift_pid = chemical_shift_pid
        self.distance_threshold = distance_threshold
        self.prochiral_tolerance = prochiral_tolerance
        self.h_tolerance = h_tolerance
        self.c_tolerance = c_tolerance
        self.n_tolerance = n_tolerance
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

    # def get_ambig_noe_conn(self, shift_pairs, peak_df, tol_dict, diagonal_tol=0.3):
    #     """Enhanced ambiguous peak-resonance mapping with v2 features"""
    #     connections = defaultdict(list)
    #     peak_assignments = defaultdict(list)
    #
    #     # Validate and format points array
    #     if not shift_pairs:
    #         print("⚠️ No shift pairs available for KDTree construction")
    #         return connections, peak_assignments
    #
    #     try:
    #         # Create 2D array of shape (N, 2)
    #         points = np.array([[pair['shift1'], pair['shift2']] for pair in shift_pairs])
    #         print(f"✅ KDTree points shape: {points.shape}")  # Should be (N, 2)
    #
    #         # Build KDTree with periodic boundary conditions for aliasing
    #         tree = KDTree(points, boxsize=[360, 360])  # Assuming ±180ppm aliasing
    #
    #     except Exception as e:
    #         print(f"❌ KDTree creation failed: {str(e)}")
    #         print(f"Sample points: {points[:3] if len(points) > 3 else points}")
    #         return connections, peak_assignments
    #     for idx, peak in peak_df.iterrows():
    #         if pd.isna(peak['position_F1']) or pd.isna(peak['position_F2']):
    #             continue
    #
    #         # Get dynamic tolerance from peak dimensions
    #         # there is an error here..
    #         #tol1 = tol_dict[peak['dim_F1']][0]
    #         #tol2 = tol_dict[peak['dim_F2']][0]
    #         # ... easy temporary fix....
    #         tol1 = tol_dict[1][0]
    #         tol2 = tol_dict[2][0]
    #         query_tol = np.sqrt(tol1**2 + tol2**2)
    #
    #         # Find matches considering aliasing
    #         matches = tree.query_ball_point([peak['position_F1'], peak['position_F2']], query_tol)
    #
    #         # Exclude diagonal peaks
    #         if abs(peak['position_F1'] - peak['position_F2']) < diagonal_tol:
    #             continue
    #
    #         # Store ambiguous connections
    #         for i in matches:
    #             pair = shift_pairs[i]
    #             connections[(pair['atom1'], pair['atom2'])].append(peak)
    #             peak_assignments[peak['serial']].append(pair)
    #
    #     return connections, peak_assignments

    def get_ambig_noe_conn(self, shift_pairs, peak_df, tol_dict, diagonal_tol=0.3):
        """
        Enhanced ambiguous NOE peak–resonance mapping that works for 2D/3D/4D:
          • Always uses the first two dimensions (F1, F2) as the homonuclear NOE axes.
          • Builds a KDTree over (F1_shift, F2_shift) from shift_pairs.
          • For each experimental peak, uses (position_F1, position_F2) to find all
            structural pairs within sqrt(tol1^2 + tol2^2), where tol1/tol2 come from tol_dict.
          • Ignores peaks near the diagonal (|F1–F2| < diagonal_tol).
          • Returns:
              connections      : dict mapping (atom1, atom2) → list of pd.Series peaks
              peak_assignments : dict mapping peak.serial → list of matching shift_pair dicts
        """

        from scipy.spatial import KDTree
        import numpy as np
        import pandas as pd
        from collections import defaultdict

        connections = defaultdict(list)
        peak_assignments = defaultdict(list)

        # 1) Build a list of 2D points (F1, F2) from shift_pairs
        points = []
        valid_pairs = []
        for pair in shift_pairs:
            try:
                f1 = float(pair["F1_shift"])
                f2 = float(pair["F2_shift"])
            except (KeyError, TypeError, ValueError):
                # Skip if missing or non‐numeric F1/F2
                print(pair)
                continue
            points.append((f1, f2))
            valid_pairs.append(pair)

        if not points:
            print("⚠️ No valid (F1, F2) shifts available for KDTree construction")
            return connections, peak_assignments

        points_arr = np.array(points)
        try:
            # Use boxsize of [360,360] if you want ±180 ppm aliasing; omit boxsize if not needed
            tree = KDTree(points_arr, boxsize=[360, 360])
        except Exception as e:
            print(f"❌ KDTree creation failed: {e}")
            return connections, peak_assignments

        # 2) Iterate over experimental peaks
        for _, peak in peak_df.iterrows():
            # Skip if no F1 or F2 position
            if pd.isna(peak.get("position_F1")) or pd.isna(peak.get("position_F2")):
                continue

            pos1 = float(peak["position_F1"])
            pos2 = float(peak["position_F2"])

            # Exclude diagonal NOEs
            if abs(pos1 - pos2) < diagonal_tol:
                continue

            # Determine per-dimension tolerances from tol_dict
            # tol_dict keys should be dimension indices: e.g. {1: (tol1, …), 2: (tol2, …), …}
            tol1 = tol_dict.get(1, (0.0,))[0]
            tol2 = tol_dict.get(2, (0.0,))[0]
            query_tol = np.sqrt(tol1**2 + tol2**2)

            # Query KDTree
            try:
                idxs = tree.query_ball_point([pos1, pos2], query_tol)
            except Exception:
                continue

            if not idxs:
                continue

            # Map matched indices back to valid_pairs
            matched_pairs = [valid_pairs[i] for i in idxs]

            # Record connections and peak assignments
            for sp in matched_pairs:
                key = (sp["atom1"], sp["atom2"])
                connections[key].append(peak)
                peak_assignments[peak["serial"]].append(sp)

        return connections, peak_assignments


    def get_chemicalShift_of_atom(self, atom):
        model = self.ensembleModels[0]
        chainCode = atom.parent.parent.id  # Chain
        resiNum = atom.parent.id[1]  # Residue number
        atomName = atom.name  # Atom name
        chain = model[chainCode]  # Get chain A
        residue = chain[(" ", resiNum, " ")]  # (hetflag, seqnum, insertion_code)
        resiName = residue.resname.upper()

        NmrAtomPID = f"NA:{chainCode}.{resiNum}.{resiName}.{atomName}"

        try:
            nmrAtom = self.project.getByPid(NmrAtomPID)
            chemicalShift = nmrAtom.chemicalShifts[0]
            return chemicalShift.value, NmrAtomPID
        except:
            return None, None

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

        chemicalShift, NmrAtomPID = self.get_chemicalShift_of_atom( h1_neighbors[0] )

        return chemicalShift, e1, h1_neighbors[0], NmrAtomPID


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
            #print(dist_matrix)

            for i, j in combinations(range(len(protons)), 2):
                if dist_matrix[i, j] > dist_threshold:
                    continue

                a1, a2 = protons[i], protons[j]

                ctx_a1 = self._get_atom_context(a1)
                ctx_a2 = self._get_atom_context(a2)
                ctx = (ctx_a1, ctx_a2)

                if ctx in contexts:
                    noe = dist_matrix[i, j]**self.NOE_POWER

                    atom1_key = (
                        a1.parent.parent.id,
                        a1.parent.id[1],
                        a1.parent.resname,
                        a1.name
                        )
                    atom2_key = (
                        a2.parent.parent.id,
                        a2.parent.id[1],
                        a2.parent.resname,
                        a2.name
                        )

                    noe_dict[(a1, a2,atom1_key,atom2_key)] += noe



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

    def handle_prochirals(self, shift_pairs):
        """
        Remove ambiguous restraints involving prochiral protons—i.e., protons attached
        to the same heavy atom whose chemical shifts are very similar (within a specified tolerance).

        The function performs filtering in two passes:
          1. It removes ambiguous entries along atom1.
          2. Then repeats the same filtering along atom2.

        A shift pair is removed if:
          - It duplicates another entry (same atom1 and atom2), or
          - The proton (atom1 or atom2) belongs to a prochiral group—defined by:
            all protons attached to the same residue and heavy atom have similar shifts.

        Parameters:
            shift_pairs: list of dicts with keys:
                - 'atom1' and 'atom2': (chainID, resNum, resName, atomName)
                - 'F1_shift' and 'F2_shift': float shifts (proton dimensions)

        Returns:
            A filtered list of shift_pairs where redundant and prochiral duplicates are removed.
        """

        def filter_exact_duplicates(pairs):
            """
            Remove exact duplicates based on (atom1, atom2).
            Keeps the first occurrence and ignores repeats.
            """
            seen = set()
            filtered = []
            for entry in pairs:
                key = (entry['atom1'], entry['atom2'])
                if key not in seen:
                    seen.add(key)
                    filtered.append(entry)
            return filtered

        def filter_prochiral(pairs, atom_key, shift_key):
            """
            Remove prochiral duplicates for the given atom (atom1 or atom2).
            Groups entries by (chainID, resNum, resName) and filters entries whose
            proton shifts are within the prochiral_tolerance.
            Keeps only one representative from each prochiral group.

            Parameters:
                pairs: atom shift pair table
                atom_key: 'atom1' or 'atom2'
                shift_key: 'F1_shift' or 'F2_shift' (matching the atom)
            """
            grouped = defaultdict(list)

            # Group protons by residue (chain, resNum, resName)
            for entry in pairs:
                atom = entry.get(atom_key)
                shift = entry.get(shift_key)
                if atom and shift is not None:
                    key = atom[:3]  # Residue ID
                    grouped[key].append((entry, atom[3], shift))  # (entry, atomName, shift)

            filtered = []
            for key, entries in grouped.items():
                used = set()
                entries.sort(key=lambda x: x[2])  # sort by shift value
                for i, (entry_i, name_i, shift_i) in enumerate(entries):
                    if i in used:
                        continue
                    for j in range(i + 1, len(entries)):
                        _, name_j, shift_j = entries[j]
                        if abs(shift_i - shift_j) < self.prochiral_tolerance:
                            used.add(j)  # mark as redundant
                    filtered.append(entry_i)

            # Retain entries that weren’t part of the prochiral groups (e.g., missing shifts)
            for entry in pairs:
                atom = entry.get(atom_key)
                shift = entry.get(shift_key)
                if atom is None or shift is None or atom[:3] not in grouped:
                    filtered.append(entry)

            return filtered

        def filter_all(pairs):
            """
            Perform full prochiral ambiguity filtering on both atom1 and atom2.
            This includes:
              - Deduplicating entries based on atom1–atom2
              - Filtering prochiral ambiguities for atom1
              - Flipping atoms and repeating the above for atom2
            """
            # Step 1: Remove exact (atom1, atom2) duplicates
            step1 = filter_exact_duplicates(pairs)

            # Step 2: Filter out prochiral ambiguities on atom1
            step2 = filter_prochiral(step1, 'atom1', 'F1_shift')

            # Step 3: Flip atom1 and atom2 so we can treat atom2 like atom1
            step3 = filter_exact_duplicates([{
                **entry,
                'atom1'   : entry['atom2'],
                'atom2'   : entry['atom1'],
                'F1_shift': entry['F2_shift'],
                'F2_shift': entry['F1_shift']
                } for entry in step2])

            # Step 4: Filter prochiral ambiguities on what was originally atom2
            step4 = filter_prochiral(step3, 'atom1', 'F1_shift')

            # Step 5: Flip atom1 and atom2 back to original orientation
            final = [{
                **entry,
                'atom1'   : entry['atom2'],
                'atom2'   : entry['atom1'],
                'F1_shift': entry['F2_shift'],
                'F2_shift': entry['F1_shift']
                } for entry in step4]

            return final

        filtered_pairs = filter_all(shift_pairs)
        return filtered_pairs

    def generate_matched_pkl(self, peaklist_pid):
        """
        Create a new synthetic peak list in the same spectrum, containing only peaks
        that were successfully matched to structural pairs. The comment for each new
        peak will display the matched structural proton pair and hetero-atom key.
        """

        pkl = self.project.getByPid(peaklist_pid)
        matched_peaks = self.matched_peaks.get(peaklist_pid, [])

        spectrum = pkl.spectrum
        ndim = spectrum.dimensionCount  # number of F‐dimensions in this spectrum

        # Create a new synthetic peak list
        matched_list = spectrum.newPeakList(
                title=f"RPF_Matched-{pkl.pid}",
                isSynthetic=True,
                comment='Peaks matched based on assignment and structure'
                )

        for mp in matched_peaks:
            # 'mp' is a dict with keys: 'peak', 'matches', 'best_match', 'indices'
            exp_peak = mp['peak']  # this is a CIP‐peak object with .position_F1, .position_F2, …

            best = mp['best_match']  # this is one of our shift_pair dicts

            # Build a comment that shows atom1–atom2–heteroKey
            # In the new shift_pair format, these keys are named 'atom1', 'atom2', and 'heteroKey1'
            # …inside generate_matched_pkl, after `best = mp['best_match']`…
            prot1 = best.get('atom1')
            prot2 = best.get('atom2')
            hetero1 = best.get('heteroKey1')
            hetero2 = best.get('heteroKey2')

            # Build a list of “F#=value” strings for every dimension
            shift_strings = []
            for dim in range(1, ndim + 1):
                key = f"F{dim}_shift"
                val = best.get(key)
                if val is None:
                    shift_strings.append(f"F{dim}=NA")
                else:
                    shift_strings.append(f"F{dim}={val:.3f}")
            shifts_joined = ", ".join(shift_strings)

            # Build hetero part of comment: always include hetero1; include hetero2 only if not None
            if hetero2:
                hetero_part = f"hetero1={hetero1}; hetero2={hetero2}"
            else:
                hetero_part = f"hetero={hetero1}"

            comment_text = (
                f"Matched: {prot1}-{prot2}; "
                f"{shifts_joined}; "
                f"{hetero_part}"
            )

            # Collect the ppmPositions array in the correct order (F1, F2, …)
            ppm_positions = []
            for dim in range(1, ndim + 1):
                attr = f"position_F{dim}"
                # Some older peaks might not have all dims; default to None or 0
                val = getattr(exp_peak, attr, None)
                ppm_positions.append(val)

            # Create the new matched peak in the synthetic list
            matched_list.newPeak(
                    ppmPositions=ppm_positions,
                    comment=comment_text,
                    volume=exp_peak.volume,
                    height=exp_peak.height
                    )

    def generate_synthetic_peaks(self, peakList_pid):
        """Create synthetic peaks for missing pairs with validation across 2D, 3D, and 4D experiments."""
        synth_peaks = []
        peakListObject = self.project.getByPid(peakList_pid)
        if not peakListObject:
            raise ValueError(f"Peak list with PID {peakList_pid} not found")

        spectrum = peakListObject.spectrum
        ndim = spectrum.dimensionCount  # Number of F‐dimensions (2, 3, or 4)

        try:
            missing_pairs = self.missing_pairs.get(peakListObject.pid, [])
            print(f"🔄 Generating synthetic peaks for {peakList_pid} ({len(missing_pairs)} missing pairs)")

            # Build set of required keys: F1_shift, F2_shift, …, F{ndim}_shift, plus atom1/atom2
            required_keys = {f"F{d}_shift" for d in range(1, ndim + 1)}
            required_keys.update({'atom1', 'atom2'})

            # In 3D or 4D there may be heteroKey and heteroKey2; include those if ndim >= 3
            if ndim >= 3:
                required_keys.add('heteroKey1')
            if ndim == 4:
                required_keys.add('heteroKey2')

            for i, pair in enumerate(missing_pairs):
                if not isinstance(pair, dict):
                    print(f"⚠️ Invalid pair format at index {i}: {type(pair)}")
                    continue

                # Check presence of all required keys
                if not required_keys.issubset(pair.keys()):
                    missing = required_keys - pair.keys()
                    print(f"⚠️ Missing keys in pair {i}: {missing}")
                    continue

                # Build the synthetic peak dictionary
                peak = {}
                # Fill position_F1 … position_F{ndim} from F?_shift
                for d in range(1, ndim + 1):
                    peak[f"position_F{d}"] = pair.get(f"F{d}_shift", None)

                peak['volume'] = 0

                # Build comment text based on dimensionality
                prot1 = pair.get('atom1')
                prot2 = pair.get('atom2')
                hetero1 = pair.get('heteroKey1')
                hetero2 = pair.get('heteroKey2') if ndim == 4 else None

                if ndim == 2:
                    comment_text = f"Missing: {prot1}-{prot2}"
                elif ndim == 3:
                    comment_text = f"Missing: {prot1}-{prot2}; hetero={hetero1}"
                else:  # ndim == 4
                    # hetero2 should always be present in 4D missing pairs
                    comment_text = f"Missing: {prot1}-{prot2}; hetero1={hetero1}; hetero2={hetero2}"

                peak['comment'] = comment_text
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
        elif original_name == 'HH11' and residue_type in ['ARG']:
            expanded = ['HZ1%', 'HZ11', 'HZ%']
        elif original_name == 'HH12' and residue_type in ['ARG']:
            expanded = ['HZ1%', 'HZ12', 'HZ%']
        elif original_name == 'HH21' and residue_type in ['ARG']:
            expanded = ['HZ2%', 'HZ21', 'HZ%']
        elif original_name == 'HH22' and residue_type in ['ARG']:
            expanded = ['HZ2%', 'HZ22', 'HZ%']

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
        for col in chemicalShift_list_df.columns:
            chemicalShift_list_df[col] = pd.to_numeric(chemicalShift_list_df[col], errors='ignore')

        shift_map = defaultdict(dict)

        for _, row in chemicalShift_list_df.iterrows():
            key = (
                row['chainCode'],
                row['sequenceCode'],
                row['residueType'],
                row['atomName'])
            shift_map[key] = row['value']

        self.shift_map = shift_map

    import re

    def get_unmatched_peaks_report(self, unexplained_peaks, missing_pairs):
        """Generate detailed report of unmatched peaks and missing pairs, handling 2D/3D/4D."""
        report = {
            'unexplained_peaks': [],
            'missing_pairs'    : []
            }

        def _sorted_dim_keys(prefix, data_keys):
            """
            Find all keys in data_keys that begin with `prefix` + digits,
            e.g. prefix='F' matches 'F1_shift', extracts 1, etc.
            Returns a list of keys sorted by the extracted integer.
            """
            matches = []
            pattern = re.compile(rf'^{re.escape(prefix)}(\d+)')
            for k in data_keys:
                m = pattern.match(k)
                if m:
                    idx = int(m.group(1))
                    matches.append((idx, k))
            matches.sort(key=lambda x: x[0])
            return [k for (_, k) in matches]

        # 1) Unexplained peaks
        for peak in unexplained_peaks:
            # peak is a pandas Series with keys like 'position_F1', 'position_F2', ...
            dim_keys = _sorted_dim_keys('position_F', peak.keys())
            serial = peak.get('serial', None)
            entry = {
                'peak_id': int(serial) if serial is not None else None,
                'volume' : peak.get('volume', None),
                'snr'    : peak.get('signalToNoiseRatio', None)
                }
            # Add each position_F# dynamically
            for key in dim_keys:
                entry[key] = peak.get(key, None)
            report['unexplained_peaks'].append(entry)

        # 2) Missing structural pairs
        for pair in missing_pairs:
            atom1 = pair.get('atom1')
            atom2 = pair.get('atom2')
            atom1_str = (
                f"{atom1[0]}{atom1[1]}{atom1[3]}"
                if isinstance(atom1, (tuple, list)) else str(atom1)
            )
            atom2_str = (
                f"{atom2[0]}{atom2[1]}{atom2[3]}"
                if isinstance(atom2, (tuple, list)) else str(atom2)
            )

            # Collect all 'F#_shift' keys in this pair and sort them
            dim_keys = _sorted_dim_keys('F', pair.keys())
            # Keep only those ending with '_shift'
            dim_keys = [k for k in dim_keys if k.endswith('_shift')]

            entry = {
                'atom1'   : atom1_str,
                'atom2'   : atom2_str,
                'distance': pair.get('distance', None)
                }
            # Add expected shifts for each dimension, renamed to 'expected_F#'
            for k in dim_keys:
                # Extract the dimension number from 'F{n}_shift'
                try:
                    idx = int(re.match(r'^F(\d+)', k).group(1))
                except Exception:
                    continue
                entry[f'expected_F{idx}'] = pair.get(k, None)

            report['missing_pairs'].append(entry)

        return report

    def build_shift_pair_from_mapping(self, mapping, dims, distance):
        """
        Given:
          • mapping: a dict where each dimension label (from dims) maps to [atom_key, shift_value].
               Example: {
                 'Hn': [(chain, resnum, resname, 'HN'), 8.996],
                 'Nh': [(chain, resnum, resname, 'N'), 121.358],
                 'H' : [(chain2, resnum2, resname2, 'H'), 8.996]
               }
          • dims:    a tuple of dimension labels, e.g. ('Hn', 'H', 'Nh').
          • distance: a float for the inter‐proton distance.

        Returns:
          A shift_pair dict with keys:
            'distance',
            'F1_shift', 'F2_shift', …, 'F{ndim}_shift',
            'atom1', 'atom2',
            'heteroKey1', 'heteroKey2'
          where:
            • atom1 = the first proton (first 'H' in dims),
            • atom2 = the second proton (second 'H'), or None if absent,
            • heteroKey1 = hetero bound to atom1 (matched by residue), or None,
            • heteroKey2 = hetero bound to atom2, or None.
        """
        ndim = len(dims)

        # 1) Identify H dims in order
        h_dims = [d for d in dims if d.startswith('H')]
        atom1 = mapping[h_dims[0]][0] if len(h_dims) >= 1 and h_dims[0] in mapping else None
        atom2 = mapping[h_dims[1]][0] if len(h_dims) >= 2 and h_dims[1] in mapping else None

        # 2) Identify hetero dims (non‐H)
        non_h_dims = [d for d in dims if not d.startswith('H')]

        # Determine heteroKey1 by matching residue number to atom1
        heteroKey1 = mapping[dims[2]][0] if len(dims) >= 3 else None
        heteroKey2 = mapping[dims[3]][0] if len(dims) >= 4 else None

        # 3) Build the F‐shift dictionary
        shifts = {}
        for i, dim_label in enumerate(dims):
            shift_val = mapping.get(dim_label, [None, None])[1]
            shifts[f"F{i + 1}_shift"] = shift_val

        # 4) Assemble final dict
        shift_pair = {
            'distance'  : distance,
            'atom1'     : atom1,
            'atom2'     : atom2,
            'heteroKey1': heteroKey1,
            'heteroKey2': heteroKey2,
            **shifts
            }

        return shift_pair

    def get_structural_shift_pairs(self, proton_dists_conn, spectrum):
        """
        For 2D/3D/4D, return shift_pair dicts with exactly these keys:
          'atom1', 'atom2', 'distance',
          'F1_shift', 'F2_shift', …, 'F{ndim}_shift',
          'heteroKey1', 'heteroKey2',
          'hetResidueMismatch1', 'hetResidueMismatch2'.

        In 3D, if one proton has no assigned shift, we still produce one entry
        so that the other proton’s shift goes into its correct H‐axis. We only skip
        a combination if both s1 and s2 are None—or if they lie on the same H‐axis
        and |s1 – s2| ≤ self.h_tolerance (diagonal).
        """
        shift_pairs = []

        dims = spectrum.referenceExperimentDimensions
        ndim = len(dims)

        # 1) Build a list of all onebond channels: [(Hdim_name, HetDim_name), ...]
        onebond_pairs = []
        for mt in spectrum.magnetisationTransfers:
            if mt.transferType != 'onebond':
                continue
            d1 = dims[mt.dimension1 - 1]
            d2 = dims[mt.dimension2 - 1]
            if d1.startswith('H') and not d2.startswith('H'):
                onebond_pairs.append((d1, d2))
            elif d2.startswith('H') and not d1.startswith('H'):
                onebond_pairs.append((d2, d1))
            # else: ignore any (H,H) or (X,X) or through-space, etc.
        #print(onebond_pairs)

        for (a1, a2, key1, key2), distance in proton_dists_conn.items():

            atom1_key = [
                a1.parent.parent.id,
                a1.parent.id[1],
                a1.parent.resname,
                a1.name
                ]
            atom2_key = [
                a2.parent.parent.id,
                a2.parent.id[1],
                a2.parent.resname,
                a2.name
                ]

            possibleProton1_shifts = self._get_possible_shifts(atom1_key, self.shift_map)
            possibleProton2_shifts = self._get_possible_shifts(atom2_key, self.shift_map)

            # If no assignments exist, replace with a single (None,None) so we still generate one entry
            if not possibleProton1_shifts:
                possibleProton1_shifts = [(None, None)]
            if not possibleProton2_shifts:
                possibleProton2_shifts = [(None, None)]

            #possibleProton1_shifts
            #Out[126]: [[('A', 2, 'MET', 'H'), 8.996]]

            atom1_shift = None
            if len(possibleProton1_shifts) == 1:
                atom1_shift = possibleProton1_shifts[0][1]
                atom1_key[3] = possibleProton1_shifts[0][0][3] # set atom name to NMRatomName
            elif len(possibleProton1_shifts) > 1:
                print('likely assignment error', possibleProton1_shifts[0][1])
            else:
                atom1_shift = None

            #possibleProton2_shifts
            #Out[127]: [[('A', 4, 'GLN', 'H'), 8.384]]
            if len(possibleProton2_shifts) == 1:
                atom2_shift = possibleProton2_shifts[0][1]
                atom2_key[3] = possibleProton2_shifts[0][0][3]  # set atom name to NMRatomName
            elif len(possibleProton2_shifts) > 1:
                print('likely assingment error', possibleProton2_shifts[0][1])
            else:
                atom2_shift = None

            # 2.1) Fetch attached hetero info: (hetero_shift, element, heteroAtomObj)
            a1_het_shift, element1, heteroAtom1, _ = self.get_chemicalShift_of_attached_atom(a1)
            a2_het_shift, element2, heteroAtom2, _ = self.get_chemicalShift_of_attached_atom(a2)

            #print(a1_het_shift, element1, heteroAtom1, a2_het_shift, element2, heteroAtom2)
            #113.996 N <Atom N> 121.358 N <Atom N>

            if heteroAtom1:
                heteroKey1 = (heteroAtom1.parent.parent.id,
                              heteroAtom1.parent.id[1],
                              heteroAtom1.parent.resname,
                              heteroAtom1.name)
            else:
                heteroKey1 = None

            if heteroAtom2:
                heteroKey2 = (heteroAtom2.parent.parent.id,
                              heteroAtom2.parent.id[1],
                              heteroAtom2.parent.resname,
                              heteroAtom2.name)
            else:
                heteroKey2 = None

            mappedAtoms = []
            if len(onebond_pairs) == 0 and len(dims) == 2:  #2d proton to proton
                AtomShiftMapper_dict = {}
                AtomShiftMapper_dict[dims[0]] = [atom1_key, atom1_shift]
                AtomShiftMapper_dict[dims[1]] = [atom2_key, atom2_shift]
                mappedAtoms.append(AtomShiftMapper_dict)
                # do reciprocal cross peak
                AtomShiftMapper_dict = {}
                AtomShiftMapper_dict[dims[0]] = [atom2_key, atom2_shift]
                AtomShiftMapper_dict[dims[1]] = [atom1_key, atom1_shift]
                mappedAtoms.append(AtomShiftMapper_dict)

            # handle 3d/4d noesy containing onebond connected protons
            for attachedPair in onebond_pairs:
                AtomShiftMapper_dict = {}
                #print(attachedPair)
                missing_keys = [d for d in dims if d not in attachedPair]
                heteroAtom_of_pair = attachedPair[1][0]

                if heteroAtom_of_pair == element1:
                    AtomShiftMapper_dict[attachedPair[0]] = [atom1_key, atom1_shift]
                    AtomShiftMapper_dict[attachedPair[1]] = [heteroKey1, a1_het_shift]

                    for missing_key in missing_keys:
                        if missing_key[0] == 'H':
                            AtomShiftMapper_dict[missing_key] = [atom2_key, atom2_shift]
                        else:
                            AtomShiftMapper_dict[missing_key] = [heteroKey2, a2_het_shift]
                    mappedAtoms.append(AtomShiftMapper_dict)

                if heteroAtom_of_pair == element2:
                    AtomShiftMapper_dict[attachedPair[0]] = [atom2_key, atom2_shift]
                    AtomShiftMapper_dict[attachedPair[1]] = [heteroKey2, a2_het_shift]

                    for missing_key in missing_keys:
                        if missing_key[0] == 'H':
                            AtomShiftMapper_dict[missing_key] = [atom1_key, atom1_shift]
                        else:
                            AtomShiftMapper_dict[missing_key] = [heteroKey1, a1_het_shift]
                    mappedAtoms.append(AtomShiftMapper_dict)

            for mappedAtomSet in mappedAtoms:
                sp = self.build_shift_pair_from_mapping(mappedAtomSet, dims, distance)
                shift_pairs.append(sp)

        return shift_pairs

    def match_peaks_to_structure(self, peaklist, peak_df):
        """
        Match peaks to structural pairs using per-dimension tolerances:
          • For any F‐dimension whose label starts with 'H', use self.h_tolerance.
          • For any F‐dimension whose label starts with 'C', use self.c_tolerance.
          • For any F‐dimension whose label starts with 'N', use self.n_tolerance.

        shift_pairs: List of dicts (one or more per structural proton–proton pair), each must include:
          - 'atom1', 'atom2', 'distance'
          - 'F1_shift', 'F2_shift', ..., up to the max F‐dimension
          - 'heteroKey1', 'heteroKey2', etc.

        peak_df: DataFrame of experimental peaks with columns:
          - 'position_F1', 'position_F2', ..., up to same number of F‐dimensions.

        Returns:
          matched_peaks     : list of dicts with keys:
             'peak'         : original pd.Series for that peak
             'matches'      : list of shift_pairs dicts that match within per‐dim tolerances
             'best_match'   : the one shift_pair (among matches) with smallest 'distance'
             'indices'      : indices (into deduplicated shift_pairs) that matched
          unexplained_peaks : list of pd.Series peaks that had NaN or no matches
          missing_pairs     : list of shift_pairs that were never matched to any peak
        """

        shift_pairs = self.shift_pairs[peaklist.pid]

        spectrum = peaklist.spectrum

        # 1) Figure out how many dimensions based on peak_df columns
        pos_cols = [col for col in peak_df.columns if col.startswith("position_F")]
        pos_cols = sorted(pos_cols, key=lambda c: int(c.split("_")[1][1:]))
        ndim = len(pos_cols)
        if ndim == 0:
            raise ValueError("No position_F* columns found in peak_df.")

        # 2) Deduplicate shift_pairs by their F‐shift tuple
        unique_shift_pairs = []
        seen = set()
        for sp in shift_pairs:
            try:
                tup = tuple(float(sp[f"F{i + 1}_shift"]) for i in range(ndim))
            except (KeyError, TypeError, ValueError):
                continue
            if tup not in seen:
                seen.add(tup)
                unique_shift_pairs.append(sp)
        shift_pairs = unique_shift_pairs

        # 3) Pre‐compute each shift_pair’s F‐shift tuple and store index
        pair_shifts = []
        for idx, sp in enumerate(shift_pairs):
            try:
                vec = [float(sp[f"F{i + 1}_shift"]) for i in range(ndim)]
            except (KeyError, TypeError, ValueError):
                continue
            pair_shifts.append((idx, vec))

        if not pair_shifts:
            # No valid structural points → all peaks unexplained
            return [], peak_df.to_dict("records"), []

        # 4) Determine per‐dimension labels from any shift_pair’s keys (same dims apply to all)
        #    We need the experiment’s dimension names in order:
        dims = spectrum.referenceExperimentDimensions  # assumed stored somewhere accessible
        # e.g. dims = ('H', 'H1', 'C', 'N1', ...)
        # We truncate/pad to ndim, just in case:
        dims = list(dims[:ndim])

        matched_peaks = []
        unexplained_peaks = []

        # 5) For each experimental peak, check per‐dimension tolerance
        for _, peak in peak_df.iterrows():
            # 5.a) Skip if any NaN in required columns
            if any(pd.isna(peak[col]) for col in pos_cols):
                unexplained_peaks.append(peak)
                continue

            # 5.b) Build numeric list of peak positions
            try:
                pvec = [float(peak[col]) for col in pos_cols]
            except (KeyError, TypeError, ValueError):
                unexplained_peaks.append(peak)
                continue

            # 5.c) For each shift_pair, check all dims
            matches = []
            for idx, svec in pair_shifts:
                ok = True
                for i in range(ndim):
                    dim_label = dims[i]  # e.g. 'H', 'H1', 'C', 'N1'
                    peak_val = pvec[i]
                    shift_val = svec[i]
                    if dim_label.startswith('H'):
                        tol = self.h_tolerance
                    elif dim_label.startswith('C'):
                        tol = self.c_tolerance
                    elif dim_label.startswith('N'):
                        tol = self.n_tolerance
                    else:
                        tol = self.h_tolerance  # default if unknown
                    if abs(peak_val - shift_val) > tol:
                        ok = False
                        break
                if ok:
                    matches.append(idx)

            if not matches:
                unexplained_peaks.append(peak)
                continue

            # 5.d) Build list of matching shift_pair dicts
            matched_idxs = set(matches)
            match_dicts = [shift_pairs[i] for i in matched_idxs]
            best = min(match_dicts, key=lambda sp: sp.get("distance", float("inf")))

            matched_peaks.append({
                "peak"      : peak,
                "matches"   : match_dicts,
                "best_match": best,
                "indices"   : list(matched_idxs)
                })

        # 6) Determine which structural pairs were never matched
        all_matched = set()
        for mp in matched_peaks:
            all_matched.update(mp["indices"])

        missing_pairs = [
            sp for i, sp in enumerate(shift_pairs) if i not in all_matched
            ]

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
                isSynthetic=True,
                comment='Missing peaks based on assignment and structure'
                )
        # Add peaks from syntheticPeaks DataFrame
        for _, row in self.syntheticPeaks[peaklist_pid].iterrows():
            synthetic_list.newPeak(
                    ppmPositions=[row[f'position_F{dim}'] for dim in range(1, spectrum.dimensionCount + 1)],
                    comment=row['comment']
                    )

    def run_rpf_analysis(self):
        # 1) Load data
        self.peaklists = self.getPeakList_objects()
        print(f"Loaded {len(self.peaklists)} peak lists")

        self.ensemble = self.load_pdb_ensemble()
        self.ensembleModels = [model for model in self.ensemble]
        print(f"PDB ensemble models: {len(self.ensembleModels)}")

        self.chain = self.get_chain()
        print(f"Working on chain: {self.chain}")

        # 2) Build shift_map
        self.create_shift_map()
        print(f"Shift map contains {len(self.shift_map)} entries")

        # 3) Process each peak list
        for pkl in self.peaklists:
            print("\n--- Processing PeakList:", pkl.pid, "---")
            peak_df = pkl.getAsDataFrame()
            print(f"PeakList {pkl.pid}: {len(peak_df)} experimental peaks (first few):")
            print(peak_df.head())

            spectrum = pkl.spectrum
            print("Spectrum dimensions:", spectrum.referenceExperimentDimensions)

            # 3a) Build proton_dists_conn
            contexts = self.get_spectrum_context(spectrum)
            proton_dists_conn = self.get_proton_dists_conn(
                    contexts, self.distance_threshold
                    )
            print(f"Found {len(proton_dists_conn)} structural proton–proton pairs")
            # Show a couple of entries:
            for i, ((a1, a2, k1,k2), dist) in enumerate(proton_dists_conn.items()):
                if i >= 3:
                    break
                key1 = (a1.parent.parent.id, a1.parent.id[1], a1.parent.resname, a1.name)
                key2 = (a2.parent.parent.id, a2.parent.id[1], a2.parent.resname, a2.name)
                print(f"  Pair {i}: {key1} – {key2}, dist = {dist:.2f}")

            self.proton_dists_conn[pkl.pid] = proton_dists_conn

            # 3b) Convert to shift_pairs
            shift_pairs = self.get_structural_shift_pairs(
                    proton_dists_conn, spectrum=spectrum
                    )
            print(f"After get_structural_shift_pairs: {len(shift_pairs)} shift_pairs")
            if shift_pairs:
                print("  Example shift_pair[0]:", shift_pairs[0])

            # 3c) Handle prochirals
            # (drop the extra spectrum argument)
            shift_pairs_filtered = self.handle_prochirals(shift_pairs)
            print(f"After handle_prochirals: {len(shift_pairs_filtered)} shift_pairs remain")
            if shift_pairs_filtered:
                print("  Example filtered shift_pair[0]:", shift_pairs_filtered[0])

            self.shift_pairs[pkl.pid] = shift_pairs_filtered

            # 3d) Build a (dummy) tolerance dictionary
            # We'll still use get_ambig_noe_conn later, but for debugging, show it here
            tol_dict = {}
            for dimObj in spectrum.dimensions:
                # dimObj is a DataDim; assign a tuple for (tol, dummy)
                tol_dict[dimObj] = (0.05, 0.05)
            print("Tolerance dictionary keys:", list(tol_dict.keys()))

            # 3e) Perform ambiguous NOE mapping (just to see connections)
            connections, peak_assignments = self.get_ambig_noe_conn(
                    self.shift_pairs[pkl.pid], peak_df, tol_dict
                    )
            print(f"get_ambig_noe_conn → {len(connections)} distinct (atom1,atom2) keys")
            # Show a couple of keys:
            for i, key in enumerate(connections):
                if i >= 3:
                    break
                print(f"  Conn key[{i}]: {key}, #peaks={len(connections[key])}")

            # 3f) Match peaks to structure
            matched_peaks, unexplained_peaks, missing_pairs = self.match_peaks_to_structure(
                    pkl, peak_df
                    )
            print(f"match_peaks_to_structure → matched {len(matched_peaks)} peaks, "
                  f"unexplained {len(unexplained_peaks)}, missing {len(missing_pairs)}")

            if matched_peaks:
                print("  Example matched_peaks[0]['peak'].serial =", matched_peaks[0]["peak"]["serial"])
                print("  Example matched_peaks[0]['best_match']:", matched_peaks[0]["best_match"])

            if unexplained_peaks:
                print("  Example unexplained_peaks[0]['serial'] =", unexplained_peaks[0]["serial"])

            if missing_pairs:
                print("  Example missing_pairs[0]:", missing_pairs[0])

            self.matched_peaks[pkl.pid] = matched_peaks
            self.unexplained_peaks[pkl.pid] = unexplained_peaks
            self.missing_pairs[pkl.pid] = missing_pairs

            # 3g) Generate synthetic peaks
            synth_df = self.generate_synthetic_peaks(pkl.pid)
            print(f"Generated synthetic peaks DataFrame with shape {synth_df.shape}")
            self.syntheticPeaks[pkl.pid] = synth_df

            # 3h) Per-residue scores
            if matched_peaks:
                scores = self.calculate_per_residue_scores(matched_peaks, shift_pairs_filtered)
                print(f"Per-residue scores: {len(scores)} entries")
            else:
                scores = {}
                print("No matched peaks → skipping per-residue scores")
            self.residue_scores[pkl.pid] = scores

            # 3i) Unmatched‐peaks report
            report = self.get_unmatched_peaks_report(unexplained_peaks, missing_pairs)
            print("Unmatched‐peaks report keys:", report.keys())
            self.unmatched_peaks_report[pkl.pid] = report

            # 3j) RPF metrics
            metrics_dict = self.calculate_rpf_metrics(matched_peaks, unexplained_peaks, missing_pairs)
            print("RPF metrics:", metrics_dict)
            self.metrics[pkl.pid] = metrics_dict

        # 4) Merge results across all peak lists
        merged_matched, merged_unexplained, merged_missing = self.merge_rpf_results(
                self.matched_peaks, self.unexplained_peaks, self.missing_pairs
                )
        print(
            f"\nOverall: matched {len(merged_matched)}, unexplained {len(merged_unexplained)}, missing {len(merged_missing)}")

        # 5) Generate overall reports
        overall_report = self.get_unmatched_peaks_report(merged_unexplained, merged_missing)
        print("Overall unmatched‐peaks report:", overall_report.keys())
        self.unmatched_peaks_report["overall"] = overall_report

        overall_metrics = self.calculate_rpf_metrics(merged_matched, merged_unexplained, merged_missing)
        print("Overall RPF metrics:", overall_metrics)
        self.metrics["overall"] = overall_metrics


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

import pandas as pd


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
        self.browseButton = Button(
            self.mainWidget,
            text="Browse",
            callback=self._browsePDB,
            grid=(row, 2),
        )

        row += 1
        # Distance Tolerance
        Label(self.mainWidget, text="Distance Tolerance:", grid=(row, 0))
        self.distTolerance = Entry.FloatEntry(
            self.mainWidget, grid=(row, 1), text=4.8
        )

        row += 1
        # H / C / N tolerances grouped in a Frame
        tolFrame = Frame(parent=self.mainWidget, setLayout=True, grid=(row, 0), gridSpan=(1, 3))
        # Inside tolFrame, place the three labels and entries
        tolRow = 0
        Label(tolFrame, text="Tolerances H:", grid=(tolRow, 0))
        self.hTolerance = Entry.FloatEntry(tolFrame, grid=(tolRow, 1), text=0.04)

        #tolRow += 1
        Label(tolFrame, text="C:", grid=(tolRow, 2))
        self.cTolerance = Entry.FloatEntry(tolFrame, grid=(tolRow, 3), text=0.3)

        #tolRow += 1
        Label(tolFrame, text="N:", grid=(tolRow, 4))
        self.nTolerance = Entry.FloatEntry(tolFrame, grid=(tolRow, 5), text=0.3)

        row += 1
        # Chemical Shift List
        self.chemicalShiftListWidget = ChemicalShiftListPulldown(
            self.mainWidget, grid=(row, 0), gridSpan=(1, 2)
        )

        row += 1
        # PeakLists
        Label(self.mainWidget, text="PeakLists", grid=(row, 0), hAlign="right")
        self.peakListWidget = ListWidgetPair(self.mainWidget, grid=(row, 1), gridSpan=(1, 1))
        self.peakListWidget.setMaximumHeight(100)

        row += 1
        # Buttons
        self.ButtonFrame = Frame(
            parent=self.mainWidget, setLayout=True, grid=(row, 0), gridSpan=(1, 3)
        )
        self.analyzeButton = Button(
            self.ButtonFrame, text="Analyze", callback=self._runAnalysis, grid=(0, 0)
        )
        self.closeButton = Button(
            self.ButtonFrame, text="Close", callback=self.reject, grid=(0, 1)
        )

        self._populateWidgets()

    def _populateWidgets(self):
        # Get NOESY PeakLists
        _peakLists = []
        for pkl in self.project.peakLists:
            transfers = [transfer.transferType for transfer in pkl.spectrum.magnetisationTransfers]
            if "through-space" in transfers:
                _peakLists.append(pkl)

        self.peakListWidget.setListObjects(_peakLists)

    def _browsePDB(self):
        _currentPath = self.pdbLineEdit.getText()
        if _currentPath:
            _directory = aPath(_currentPath).parent.asString()
        else:
            _directory = self.project.application.dataPath.asString()

        dialog = OtherFileDialog(
            parent=self.mainWindow,
            _useDirectoryOnly=False,
            directory=_directory,
            fileFilter="*.pdb",
        )
        dialog._show()
        if (path := dialog.selectedFile()) is not None:
            self.pdbLineEdit.setText(str(path))

    def _runAnalysis(self):
        pdb_path = self.pdbLineEdit.getText()
        if not pdb_path:
            MessageDialog.showWarning("No PDB File", "Please select a PDB file.")
            return

        chemical_shift_pid = self.chemicalShiftListWidget.getText()
        peak_list_pids = self.peakListWidget.getRightList()
        project = self.project

        distance_threshold = float(self.distTolerance.getText())
        h_tol = float(self.hTolerance.getText())
        c_tol = float(self.cTolerance.getText())
        n_tol = float(self.nTolerance.getText())

        # Instantiate the calculator with per-element tolerances
        self.rpf_calc = RpfCalculator(
            project=project,
            pdb_path=pdb_path,
            peak_list_pids=peak_list_pids,
            chemical_shift_pid=chemical_shift_pid,
            distance_threshold=distance_threshold,
            h_tolerance=h_tol,
            c_tolerance=c_tol,
            n_tolerance=n_tol,
        )

        self.rpf_calc.run_rpf_analysis()

        # Create overall metrics table
        df_overall = pd.DataFrame.from_dict(
            self.rpf_calc.metrics["overall"], orient="index", columns=["value"]
        )
        df_overall.reset_index(inplace=True)
        df_overall.columns = ["metric", "value"]
        self.project.newDataTable(name="RPF_Metrics_Overall", data=df_overall)

        # Process each PeakList
        for peak_list_pid in peak_list_pids:
            self.rpf_calc.create_synthetic_peaklist(peaklist_pid=peak_list_pid)
            self.rpf_calc.generate_matched_pkl(peaklist_pid=peak_list_pid)

            df_metrics = pd.DataFrame.from_dict(
                self.rpf_calc.metrics["overall"], orient="index", columns=["value"]
            )
            df_metrics.reset_index(inplace=True)
            df_metrics.columns = ["metric", "value"]
            self.project.newDataTable(name="RPF_Metrics_" + peak_list_pid, data=df_metrics)

            df_res = pd.DataFrame.from_dict(
                self.rpf_calc.residue_scores[peak_list_pid], orient="index"
            )
            df_res.index = pd.MultiIndex.from_tuples(
                df_res.index, names=["chain", "residue_number", "residue_name"]
            )
            df_res.reset_index(inplace=True)
            self.project.newDataTable(name="RPF_PerResidueScores_" + peak_list_pid, data=df_res)

            columns = ['atom1', 'atom2', 'heteroKey1', 'heteroKey2', 'F1_shift', 'F2_shift', 'F3_shift','distance']
            df = pd.DataFrame(self.rpf_calc.shift_pairs[peak_list_pid],
                                                        columns=columns)
            df.index = range(len(df))
            self.project.newDataTable(name="shift_pairs_" + peak_list_pid,
                                      data=pd.DataFrame(df,
                                                        columns=columns))


if __name__ == '__main__':
    # from ccpn.ui.gui.widgets.Application import TestApplication
    #
    # app = TestApplication()
    popup = RpfGui()
    popup.exec_()