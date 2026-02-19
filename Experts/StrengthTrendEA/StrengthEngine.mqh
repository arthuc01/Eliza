#ifndef __STRENGTH_ENGINE_MQH__
#define __STRENGTH_ENGINE_MQH__

#define CURR_COUNT 8

struct StrengthCandidate
{
   string symbol;
   bool   is_long;
   string strong_ccy;
   string weak_ccy;
   double divergence;
   double divergence_slope;
   double score;
};

class CStrengthEngine
{
private:
   string m_currencies[CURR_COUNT];
   string m_pairs[];
   double m_strength[CURR_COUNT];
   double m_dispersion;
   int    m_top_idx[2];
   int    m_bot_idx[2];

   int CurrencyIndex(const string c) const
   {
      for(int i=0;i<CURR_COUNT;i++) if(m_currencies[i]==c) return i;
      return -1;
   }

   bool PairCurrencies(const string symbol, string &base, string &quote) const
   {
      if(StringLen(symbol) < 6) return false;
      base  = StringSubstr(symbol,0,3);
      quote = StringSubstr(symbol,3,3);
      if(CurrencyIndex(base)<0 || CurrencyIndex(quote)<0) return false;
      return true;
   }

   bool FindTradablePair(const string base_ccy, const string quote_ccy, string &chosen, bool &is_long) const
   {
      chosen="";
      is_long=true;
      for(int i=0;i<ArraySize(m_pairs);i++)
      {
         string s=m_pairs[i];
         if(!SymbolInfoInteger(s,SYMBOL_EXIST)) continue;
         string b,q;
         if(!PairCurrencies(s,b,q)) continue;
         if(b==base_ccy && q==quote_ccy) { chosen=s; is_long=true; return true; }
         if(b==quote_ccy && q==base_ccy) { chosen=s; is_long=false; return true; }
      }
      return false;
   }

   double CurrencyStrengthTFShift(const string currency,
                                  const ENUM_TIMEFRAMES tf,
                                  const int lookback,
                                  const int atr_period,
                                  const int shift)
   {
      double sum = 0.0;
      int count = 0;
      for(int i=0;i<ArraySize(m_pairs);i++)
      {
         string s = m_pairs[i];
         if(!SymbolInfoInteger(s,SYMBOL_EXIST)) continue;
         if(!SymbolSelect(s,true)) continue;

         MqlRates rates[];
         int need = shift + lookback + atr_period + 5;
         if(CopyRates(s, tf, 0, need, rates) < need) continue;
         ArraySetAsSeries(rates,true);

         int hATR = iATR(s, tf, atr_period);
         if(hATR==INVALID_HANDLE) continue;
         double atrv[];
         ArraySetAsSeries(atrv,true);
         if(CopyBuffer(hATR,0,0,shift+3,atrv)<shift+3)
         {
            IndicatorRelease(hATR);
            continue;
         }
         IndicatorRelease(hATR);

         if(atrv[shift] <= 0.0) continue;

         double mom = (rates[shift].close - rates[shift+lookback].close) / atrv[shift];
         string b,q;
         if(!PairCurrencies(s,b,q)) continue;

         if(b==currency) { sum += mom; count++; }
         else if(q==currency) { sum -= mom; count++; }
      }
      if(count==0) return 0.0;
      return sum / count;
   }

public:
   void Init(const string symbol_list)
   {
      m_currencies[0]="USD"; m_currencies[1]="EUR"; m_currencies[2]="GBP"; m_currencies[3]="JPY";
      m_currencies[4]="CHF"; m_currencies[5]="CAD"; m_currencies[6]="AUD"; m_currencies[7]="NZD";
      ArrayResize(m_pairs,0);

      string raw = symbol_list;
      StringReplace(raw," ","");
      string tokens[];
      int n=StringSplit(raw,',',tokens);
      for(int i=0;i<n;i++)
      {
         if(StringLen(tokens[i])<6) continue;
         int new_size=ArraySize(m_pairs)+1;
         ArrayResize(m_pairs,new_size);
         m_pairs[new_size-1]=tokens[i];
      }
   }

   bool Compute(const int lookback, const int atr_period,
                const ENUM_TIMEFRAMES tf_h4, const ENUM_TIMEFRAMES tf_h1, const ENUM_TIMEFRAMES tf_m15,
                const double w_h4, const double w_h1, const double w_m15)
   {
      for(int i=0;i<CURR_COUNT;i++)
      {
         double s4 = CurrencyStrengthTFShift(m_currencies[i], tf_h4, lookback, atr_period, 0);
         double s1 = CurrencyStrengthTFShift(m_currencies[i], tf_h1, lookback, atr_period, 0);
         double sM = CurrencyStrengthTFShift(m_currencies[i], tf_m15, lookback, atr_period, 0);
         m_strength[i] = w_h4*s4 + w_h1*s1 + w_m15*sM;
      }

      int idx[];
      ArrayResize(idx,CURR_COUNT);
      for(int i=0;i<CURR_COUNT;i++) idx[i]=i;
      for(int i=0;i<CURR_COUNT-1;i++)
         for(int j=i+1;j<CURR_COUNT;j++)
            if(m_strength[idx[j]] > m_strength[idx[i]])
            {
               int t=idx[i]; idx[i]=idx[j]; idx[j]=t;
            }

      m_top_idx[0]=idx[0]; m_top_idx[1]=idx[1];
      m_bot_idx[0]=idx[CURR_COUNT-1]; m_bot_idx[1]=idx[CURR_COUNT-2];
      m_dispersion = m_strength[idx[0]] - m_strength[idx[CURR_COUNT-1]];
      return true;
   }

   double Dispersion() const { return m_dispersion; }
   double StrengthByCurrency(const string c) const { int i=CurrencyIndex(c); if(i<0) return 0.0; return m_strength[i]; }
   string TopCurrency(const int rank) const { return (rank>=0 && rank<2)?m_currencies[m_top_idx[rank]]:""; }
   string BottomCurrency(const int rank) const { return (rank>=0 && rank<2)?m_currencies[m_bot_idx[rank]]:""; }

   string RankingString() const
   {
      int idx[];
      ArrayResize(idx,CURR_COUNT);
      for(int i=0;i<CURR_COUNT;i++) idx[i]=i;
      for(int i=0;i<CURR_COUNT-1;i++)
         for(int j=i+1;j<CURR_COUNT;j++)
            if(m_strength[idx[j]] > m_strength[idx[i]])
            {
               int t=idx[i]; idx[i]=idx[j]; idx[j]=t;
            }
      string out="";
      for(int k=0;k<CURR_COUNT;k++)
      {
         if(k>0) out += " | ";
         out += m_currencies[idx[k]] + ":" + DoubleToString(m_strength[idx[k]],2);
      }
      return out;
   }

   int BuildCandidates(StrengthCandidate &out_candidates[],
                       const int max_candidates,
                       const int max_spread_points,
                       const double atr_ratio_min,
                       const ENUM_TIMEFRAMES confirm_tf,
                       const ENUM_TIMEFRAMES slope_tf,
                       const int slope_lookback_bars,
                       const int strength_lookback,
                       const int strength_atr_period,
                       const double div_min,
                       const double div_slope_min,
                       const double div_slope_weight)
   {
      ArrayResize(out_candidates,0);
      string strongs[2] = {m_currencies[m_top_idx[0]], m_currencies[m_top_idx[1]]};
      string weaks[2] = {m_currencies[m_bot_idx[0]], m_currencies[m_bot_idx[1]]};

      for(int i=0;i<2;i++)
      {
         for(int j=0;j<2;j++)
         {
            string A = strongs[i], B = weaks[j];
            string chosen = "";
            bool is_long = true;
            if(!FindTradablePair(A,B,chosen,is_long)) continue;
            if(!SymbolSelect(chosen,true)) continue;

            long spread = SymbolInfoInteger(chosen, SYMBOL_SPREAD);
            if(spread > max_spread_points) continue;

            int hATR = iATR(chosen, confirm_tf, 14);
            if(hATR==INVALID_HANDLE) continue;
            double atrbuf[];
            ArraySetAsSeries(atrbuf,true);
            if(CopyBuffer(hATR,0,0,25,atrbuf)<25)
            {
               IndicatorRelease(hATR);
               continue;
            }
            IndicatorRelease(hATR);
            double atr_sma=0.0;
            for(int k=1;k<=20;k++) atr_sma += atrbuf[k];
            atr_sma /= 20.0;
            if(atr_sma<=0.0) continue;
            double atr_ratio = atrbuf[0] / atr_sma;
            if(atr_ratio < atr_ratio_min) continue;

            double div_now = StrengthByCurrency(A)-StrengthByCurrency(B);
            if(div_now < div_min) continue;

            double a_prev = CurrencyStrengthTFShift(A,slope_tf,strength_lookback,strength_atr_period,slope_lookback_bars);
            double b_prev = CurrencyStrengthTFShift(B,slope_tf,strength_lookback,strength_atr_period,slope_lookback_bars);
            double div_prev = a_prev - b_prev;
            double div_slope = (div_now - div_prev) / MathMax(1,slope_lookback_bars);
            if(div_slope < div_slope_min) continue;

            int n=ArraySize(out_candidates);
            ArrayResize(out_candidates,n+1);
            out_candidates[n].symbol=chosen;
            out_candidates[n].is_long=is_long;
            out_candidates[n].strong_ccy=A;
            out_candidates[n].weak_ccy=B;
            out_candidates[n].divergence = div_now;
            out_candidates[n].divergence_slope = div_slope;
            out_candidates[n].score = div_now + div_slope_weight*div_slope - 0.02*spread + 0.2*atr_ratio;
         }
      }

      for(int a=0;a<ArraySize(out_candidates)-1;a++)
         for(int b=a+1;b<ArraySize(out_candidates);b++)
            if(out_candidates[b].score > out_candidates[a].score)
            {
               StrengthCandidate t=out_candidates[a];
               out_candidates[a]=out_candidates[b];
               out_candidates[b]=t;
            }

      if(ArraySize(out_candidates) > max_candidates)
         ArrayResize(out_candidates,max_candidates);
      return ArraySize(out_candidates);
   }
};

#endif
