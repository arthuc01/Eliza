#ifndef __REGIME_DETECTOR_MQH__
#define __REGIME_DETECTOR_MQH__

enum RegimeState
{
   REGIME_RANGE = 0,
   REGIME_TREND = 1
};

class CRegimeDetector
{
private:
   RegimeState m_state;
   int m_trend_count;
   int m_range_count;

public:
   void Init()
   {
      m_state = REGIME_RANGE;
      m_trend_count = 0;
      m_range_count = 0;
   }

   RegimeState State() const { return m_state; }
   string StateString() const { return m_state==REGIME_TREND ? "TREND" : "RANGE"; }

   bool Update(const string symbol, const ENUM_TIMEFRAMES tf,
               const int ema_period, const int slope_bars, const int atr_period,
               const double slope_thresh, const double atr_ratio_thresh, const double dispersion, const double dispersion_thresh)
   {
      int hMA = iMA(symbol, tf, ema_period, 0, MODE_EMA, PRICE_CLOSE);
      int hATR = iATR(symbol, tf, atr_period);
      if(hMA==INVALID_HANDLE || hATR==INVALID_HANDLE) return false;

      double mabuf[]; ArraySetAsSeries(mabuf,true);
      double atrbuf[]; ArraySetAsSeries(atrbuf,true);

      if(CopyBuffer(hMA,0,0,slope_bars+5,mabuf)<slope_bars+5 || CopyBuffer(hATR,0,0,25,atrbuf)<25)
      {
         IndicatorRelease(hMA);
         IndicatorRelease(hATR);
         return false;
      }
      IndicatorRelease(hMA);
      IndicatorRelease(hATR);

      double atr_sma=0.0;
      for(int i=1;i<=20;i++) atr_sma += atrbuf[i];
      atr_sma /= 20.0;
      if(atrbuf[0]<=0.0 || atr_sma<=0.0) return false;

      double slope_norm = (mabuf[0] - mabuf[slope_bars]) / atrbuf[0];
      double atr_ratio = atrbuf[0] / atr_sma;

      bool trend_now = (MathAbs(slope_norm)>slope_thresh && atr_ratio>atr_ratio_thresh && dispersion>dispersion_thresh);
      if(trend_now)
      {
         m_trend_count++;
         m_range_count=0;
      }
      else
      {
         m_range_count++;
         m_trend_count=0;
      }

      if(m_trend_count>=2) m_state=REGIME_TREND;
      if(m_range_count>=2) m_state=REGIME_RANGE;
      return true;
   }

   bool TrendQuality(const string symbol, const ENUM_TIMEFRAMES tf,
                     const bool is_long,
                     const int ema_period,
                     const double adx_min,
                     const bool use_adx_filter) const
   {
      int hMA = iMA(symbol, tf, ema_period, 0, MODE_EMA, PRICE_CLOSE);
      if(hMA==INVALID_HANDLE) return false;
      double ma[]; ArraySetAsSeries(ma,true);
      if(CopyBuffer(hMA,0,0,5,ma)<5)
      {
         IndicatorRelease(hMA);
         return false;
      }
      IndicatorRelease(hMA);

      MqlRates rates[];
      if(CopyRates(symbol, tf, 0, 5, rates)<5) return false;
      ArraySetAsSeries(rates,true);
      bool ema_ok = is_long ? (rates[0].close > ma[0] && ma[0]>=ma[3]) : (rates[0].close < ma[0] && ma[0]<=ma[3]);
      if(!ema_ok) return false;

      if(!use_adx_filter) return true;
      int hADX = iADX(symbol, tf, 14);
      if(hADX==INVALID_HANDLE) return false;
      double adx[]; ArraySetAsSeries(adx,true);
      bool ok = (CopyBuffer(hADX,0,0,2,adx)>=2 && adx[0]>=adx_min);
      IndicatorRelease(hADX);
      return ok;
   }

   double AtrRatio(const string symbol, const ENUM_TIMEFRAMES tf, const int atr_period) const
   {
      int hATR=iATR(symbol,tf,atr_period);
      if(hATR==INVALID_HANDLE) return 0.0;
      double b[]; ArraySetAsSeries(b,true);
      if(CopyBuffer(hATR,0,0,25,b)<25) { IndicatorRelease(hATR); return 0.0; }
      IndicatorRelease(hATR);
      double sma=0; for(int i=1;i<=20;i++) sma+=b[i]; sma/=20.0;
      if(sma<=0) return 0.0;
      return b[0]/sma;
   }
};

#endif
