#ifndef __RISK_MANAGER_MQH__
#define __RISK_MANAGER_MQH__

class CRiskManager
{
private:
   long m_magic;

public:
   void Init(const long magic) { m_magic = magic; }

   double LotSizeByRisk(const string symbol, const double risk_pct, const double sl_distance_price)
   {
      if(sl_distance_price <= 0.0) return 0.0;

      double balance = AccountInfoDouble(ACCOUNT_BALANCE);
      double risk_money = balance * (risk_pct/100.0);

      double tick_value = SymbolInfoDouble(symbol, SYMBOL_TRADE_TICK_VALUE);
      double tick_size  = SymbolInfoDouble(symbol, SYMBOL_TRADE_TICK_SIZE);
      if(tick_value <= 0.0 || tick_size <= 0.0) return 0.0;

      double value_per_price_unit = tick_value / tick_size;
      double lots = risk_money / (sl_distance_price * value_per_price_unit);

      double vmin = SymbolInfoDouble(symbol, SYMBOL_VOLUME_MIN);
      double vmax = SymbolInfoDouble(symbol, SYMBOL_VOLUME_MAX);
      double vstep= SymbolInfoDouble(symbol, SYMBOL_VOLUME_STEP);
      if(vstep<=0) vstep=0.01;

      lots = MathFloor(lots/vstep)*vstep;
      lots = MathMax(vmin, MathMin(vmax, lots));
      return NormalizeDouble(lots, 2);
   }

   int OpenTradesTotal() const
   {
      int c=0;
      for(int i=PositionsTotal()-1;i>=0;i--)
      {
         ulong t=PositionGetTicket(i);
         if(!PositionSelectByTicket(t)) continue;
         if((long)PositionGetInteger(POSITION_MAGIC)!=m_magic) continue;
         c++;
      }
      return c;
   }

   int OpenTradesForCurrency(const string ccy) const
   {
      int c=0;
      for(int i=PositionsTotal()-1;i>=0;i--)
      {
         ulong t=PositionGetTicket(i);
         if(!PositionSelectByTicket(t)) continue;
         if((long)PositionGetInteger(POSITION_MAGIC)!=m_magic) continue;
         string s=PositionGetString(POSITION_SYMBOL);
         if(StringLen(s)<6) continue;
         string base=StringSubstr(s,0,3), quote=StringSubstr(s,3,3);
         if(base==ccy || quote==ccy) c++;
      }
      return c;
   }

   double TotalOpenRiskPct() const
   {
      double bal=AccountInfoDouble(ACCOUNT_BALANCE);
      if(bal<=0) return 0.0;
      double total_risk=0.0;
      for(int i=PositionsTotal()-1;i>=0;i--)
      {
         ulong t=PositionGetTicket(i);
         if(!PositionSelectByTicket(t)) continue;
         if((long)PositionGetInteger(POSITION_MAGIC)!=m_magic) continue;

         string s=PositionGetString(POSITION_SYMBOL);
         double vol=PositionGetDouble(POSITION_VOLUME);
         double open=PositionGetDouble(POSITION_PRICE_OPEN);
         double sl=PositionGetDouble(POSITION_SL);
         if(sl<=0 || vol<=0) continue;

         double tick_value = SymbolInfoDouble(s,SYMBOL_TRADE_TICK_VALUE);
         double tick_size = SymbolInfoDouble(s,SYMBOL_TRADE_TICK_SIZE);
         if(tick_size<=0) continue;
         double value_per_price = tick_value/tick_size;
         total_risk += MathAbs(open-sl)*value_per_price*vol;
      }
      return (total_risk/bal)*100.0;
   }
};

#endif
