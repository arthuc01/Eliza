#ifndef __TRADE_MANAGER_MQH__
#define __TRADE_MANAGER_MQH__

#include <Trade/Trade.mqh>

enum TrailingMode
{
   TRAIL_ATR = 0,
   TRAIL_SWING = 1
};

struct ManagedTrade
{
   ulong ticket;
   datetime entry_time;
   double entry;
   double initial_sl;
   double initial_risk_price;
   double divergence_entry;
   double dist_entry;
   double atr_ratio_entry;
   bool be_done;
   bool partial_done;
   double mfe_r;
   double mae_r;
};

class CTradeManager
{
private:
   CTrade m_trade;
   long m_magic;
   ManagedTrade m_book[];

   int FindIdx(const ulong ticket)
   {
      for(int i=0;i<ArraySize(m_book);i++) if(m_book[i].ticket==ticket) return i;
      return -1;
   }

   void EnsureEntry(const ulong ticket, const double entry, const double initial_sl,
                    const double div_e, const double dist_e, const double atr_ratio_e)
   {
      if(FindIdx(ticket)>=0) return;
      int n=ArraySize(m_book);
      ArrayResize(m_book,n+1);
      m_book[n].ticket=ticket;
      m_book[n].entry_time=TimeCurrent();
      m_book[n].entry=entry;
      m_book[n].initial_sl=initial_sl;
      m_book[n].initial_risk_price=MathAbs(entry-initial_sl);
      m_book[n].divergence_entry=div_e;
      m_book[n].dist_entry=dist_e;
      m_book[n].atr_ratio_entry=atr_ratio_e;
      m_book[n].be_done=false;
      m_book[n].partial_done=false;
      m_book[n].mfe_r=0.0;
      m_book[n].mae_r=0.0;
   }

public:
   void Init(const long magic)
   {
      m_magic=magic;
      m_trade.SetExpertMagicNumber(magic);
      ArrayResize(m_book,0);
   }

   CTrade &Trader() { return m_trade; }

   bool HasPosition(const string symbol)
   {
      if(!PositionSelect(symbol)) return false;
      return ((long)PositionGetInteger(POSITION_MAGIC)==m_magic);
   }

   bool OpenPosition(const string symbol, const bool is_long, const double lots, const double sl, const string comment)
   {
      if(lots<=0) return false;
      m_trade.SetDeviationInPoints(20);
      bool ok = is_long ? m_trade.Buy(lots,symbol,0.0,sl,0.0,comment)
                        : m_trade.Sell(lots,symbol,0.0,sl,0.0,comment);
      return ok;
   }

   void Manage(const ENUM_TIMEFRAMES m15_tf,
               const int atr_period,
               const double be_trigger_r,
               const double be_buffer_price,
               const double tp1_r,
               const double tp1_close_pct,
               const TrailingMode trail_mode,
               const double trail_k,
               const int stale_bars,
               const double stale_min_r,
               const bool exit_on_regime_flip,
               const bool regime_is_trend,
               const double min_div_for_hold,
               const double current_divergence)
   {
      for(int i=PositionsTotal()-1;i>=0;i--)
      {
         ulong ticket=PositionGetTicket(i);
         if(!PositionSelectByTicket(ticket)) continue;
         if((long)PositionGetInteger(POSITION_MAGIC)!=m_magic) continue;

         string symbol=PositionGetString(POSITION_SYMBOL);
         long type=PositionGetInteger(POSITION_TYPE);
         bool is_long=(type==POSITION_TYPE_BUY);
         double volume=PositionGetDouble(POSITION_VOLUME);
         double sl=PositionGetDouble(POSITION_SL);
         double price_open=PositionGetDouble(POSITION_PRICE_OPEN);

         int idx=FindIdx(ticket);
         if(idx<0)
         {
            EnsureEntry(ticket, price_open, sl, 0.0, 0.0, 0.0);
            idx=FindIdx(ticket);
         }

         double bid=SymbolInfoDouble(symbol,SYMBOL_BID);
         double ask=SymbolInfoDouble(symbol,SYMBOL_ASK);
         double px = is_long ? bid : ask;

         double r_price=m_book[idx].initial_risk_price;
         if(r_price<=0) continue;
         double profit_r = is_long ? (px-m_book[idx].entry)/r_price : (m_book[idx].entry-px)/r_price;

         m_book[idx].mfe_r = MathMax(m_book[idx].mfe_r, profit_r);
         m_book[idx].mae_r = MathMin(m_book[idx].mae_r, profit_r);

         if(!m_book[idx].be_done && profit_r>=be_trigger_r)
         {
            double be = is_long ? (m_book[idx].entry + be_buffer_price) : (m_book[idx].entry - be_buffer_price);
            if((is_long && (sl<be || sl==0.0)) || (!is_long && (sl>be || sl==0.0)))
            {
               m_trade.PositionModify(symbol, be, 0.0);
               m_book[idx].be_done=true;
            }
         }

         if(!m_book[idx].partial_done && profit_r>=tp1_r)
         {
            double close_vol = volume * (tp1_close_pct/100.0);
            double step=SymbolInfoDouble(symbol,SYMBOL_VOLUME_STEP);
            if(step>0) close_vol = MathFloor(close_vol/step)*step;
            double vmin=SymbolInfoDouble(symbol,SYMBOL_VOLUME_MIN);
            if(close_vol>=vmin) m_trade.PositionClosePartial(symbol, close_vol);
            m_book[idx].partial_done=true;
         }

         if(trail_mode==TRAIL_ATR)
         {
            int hATR=iATR(symbol,m15_tf,atr_period);
            if(hATR!=INVALID_HANDLE)
            {
               double a[]; ArraySetAsSeries(a,true);
               if(CopyBuffer(hATR,0,0,2,a)>=2)
               {
                  double trail = trail_k*a[0];
                  double new_sl = is_long ? (px-trail) : (px+trail);
                  if((is_long && new_sl>sl) || (!is_long && (sl==0 || new_sl<sl)))
                     m_trade.PositionModify(symbol,new_sl,0.0);
               }
               IndicatorRelease(hATR);
            }
         }

         MqlRates bars[];
         if(CopyRates(symbol,m15_tf,0,stale_bars+2,bars)>=stale_bars+2)
         {
            ArraySetAsSeries(bars,true);
            int mins=(int)((TimeCurrent()-m_book[idx].entry_time)/60);
            if(mins >= stale_bars*15 && profit_r < stale_min_r)
               m_trade.PositionClose(symbol);
         }

         if(current_divergence < min_div_for_hold)
            m_trade.PositionClose(symbol);
         if(exit_on_regime_flip && !regime_is_trend)
            m_trade.PositionClose(symbol);
      }
   }
};

#endif
