#property strict
#property version   "1.001"
#property description "Strong-vs-weak currency strength trend EA (M15 pullback)"

#include <Trade/Trade.mqh>
#include "StrengthEngine.mqh"
#include "RegimeDetector.mqh"
#include "RiskManager.mqh"
#include "TradeManager.mqh"
#include "Logger.mqh"

input long   InpMagicNumber = 551177;
input string InpSymbolsForStrength = "EURUSD,GBPUSD,AUDUSD,NZDUSD,USDJPY,USDCHF,USDCAD,EURGBP,EURJPY,EURCHF,EURAUD,EURNZD,GBPJPY,GBPCHF,GBPAUD,GBPCAD,AUDJPY,AUDNZD,AUDCAD,NZDJPY,NZDCAD,CADJPY,CHFJPY";
input ENUM_TIMEFRAMES InpStrengthTF_H4 = PERIOD_H4;
input ENUM_TIMEFRAMES InpStrengthTF_H1 = PERIOD_H1;
input ENUM_TIMEFRAMES InpStrengthTF_M15 = PERIOD_M15;
input double InpWeightH4 = 0.5;
input double InpWeightH1 = 0.3;
input double InpWeightM15 = 0.2;
input int    InpStrengthLookback = 14;
input int    InpStrengthATRPeriod = 14;
input int    InpMaxCandidates = 3;
input double InpDivMin = 0.9;
input int    InpDivSlopeBars = 3;
input double InpDivSlopeMin = 0.03;
input double InpDivSlopeWeight = 4.0;

input int    InpRegimeEMA = 50;
input int    InpRegimeSlopeBars = 10;
input double InpSlopeThresh = 0.25;
input double InpAtrRatioThresh = 1.00;
input double InpDispersionThresh = 1.20;
input bool   InpUseADXFilter = true;
input double InpADXMin = 18.0;

input int    InpMaxSpreadPoints = 35;
input double InpDistMin = 0.2;
input double InpDistMax = 1.2;
input double InpADRPctMax = 0.85;
input bool   InpEntryUseSwingBreak = true;
input bool   InpEntryUseEngulfing = true;

input bool   InpUseStructureSL = true;
input double InpATRSLMult = 1.5;
input int    InpSwingLookback = 6;
input double InpSLBufferATR = 0.2;

input double InpRiskPerTradePct = 1.0;
input int    InpMaxOpenTradesTotal = 5;
input int    InpMaxOpenTradesPerCurrency = 2;
input double InpMaxTotalRiskOpenPct = 3.0;
input double InpDailyLossLimitPct = 3.0;
input int    InpCooldownAfterLosses = 2;
input int    InpCooldownMinutes = 120;

input double InpBETriggerR = 1.0;
input double InpTP1R = 2.0;
input double InpTP1ClosePct = 50.0;
input TrailingMode InpTrailingMode = TRAIL_ATR;
input double InpTrailATRk = 1.2;
input int    InpTimeStopBars = 12;
input double InpTimeStopMinR = 0.5;
input double InpMinDivForHold = 0.8;
input bool   InpExitOnRegimeFlip = true;

input bool   InpSelfCheckMode = true;
input bool   InpEnableCSVLog = true;

enum SymbolFSM { ST_IDLE=0, ST_SETUP=1, ST_IN_TRADE=2, ST_COOLDOWN=3 };
struct SymbolState { string symbol; SymbolFSM state; datetime cooldown_until; };

CStrengthEngine g_strength;
CRegimeDetector g_regime;
CRiskManager    g_risk;
CTradeManager   g_trade_mng;
CLogger         g_log;

SymbolState g_states[];
datetime g_last_h1=0, g_last_m15=0;
datetime g_last_entry_m15_bar=0;
StrengthCandidate g_candidates[];
int g_consecutive_losses=0;
datetime g_global_cooldown_until=0;
double g_daily_realized=0.0;
int g_day_of_year=-1;

int StateIndex(const string symbol)
{
   for(int i=0;i<ArraySize(g_states);i++) if(g_states[i].symbol==symbol) return i;
   int n=ArraySize(g_states);
   ArrayResize(g_states,n+1);
   g_states[n].symbol=symbol;
   g_states[n].state=ST_IDLE;
   g_states[n].cooldown_until=0;
   return n;
}

bool IsNewBar(const string symbol, const ENUM_TIMEFRAMES tf, datetime &last_time)
{
   datetime t[];
   if(CopyTime(symbol,tf,0,1,t)<1) return false;
   if(t[0]==last_time) return false;
   last_time=t[0];
   return true;
}

double GetATR(const string symbol, ENUM_TIMEFRAMES tf, const int p)
{
   int h=iATR(symbol,tf,p);
   if(h==INVALID_HANDLE) return 0;
   double b[]; ArraySetAsSeries(b,true);
   if(CopyBuffer(h,0,0,2,b)<2) { IndicatorRelease(h); return 0; }
   IndicatorRelease(h);
   return b[0];
}

bool PullbackAndTrigger(const string symbol, const bool is_long)
{
   MqlRates r[];
   if(CopyRates(symbol,PERIOD_M15,0,20,r)<20) return false;
   ArraySetAsSeries(r,true);

   int hMA=iMA(symbol,PERIOD_M15,20,0,MODE_EMA,PRICE_CLOSE);
   if(hMA==INVALID_HANDLE) return false;
   double ma[]; ArraySetAsSeries(ma,true);
   if(CopyBuffer(hMA,0,0,20,ma)<20) { IndicatorRelease(hMA); return false; }
   IndicatorRelease(hMA);

   double atr=GetATR(symbol,PERIOD_M15,14);
   if(atr<=0) return false;

   double dist=(r[0].close-ma[0])/atr;
   if(is_long && !(dist>=InpDistMin && dist<=InpDistMax)) return false;
   if(!is_long && !(-dist>=InpDistMin && -dist<=InpDistMax)) return false;

   double day_high=iHigh(symbol,PERIOD_D1,0), day_low=iLow(symbol,PERIOD_D1,0);
   double adr=0.0;
   for(int i=1;i<=20;i++) adr += (iHigh(symbol,PERIOD_D1,i)-iLow(symbol,PERIOD_D1,i));
   adr/=20.0;
   if(adr>0 && ((day_high-day_low)/adr)>InpADRPctMax) return false;

   long spread=SymbolInfoInteger(symbol,SYMBOL_SPREAD);
   if(spread>InpMaxSpreadPoints) return false;

   bool pullback=false;
   int counter=0;
   for(int i=3;i>=1;i--)
   {
      if(is_long && r[i].close<r[i].open) counter++;
      if(!is_long && r[i].close>r[i].open) counter++;
   }
   if(counter>=3) pullback=true;
   if(is_long && r[1].low<=ma[1]) pullback=true;
   if(!is_long && r[1].high>=ma[1]) pullback=true;
   if(!pullback) return false;

   bool trigger=false;
   if(InpEntryUseSwingBreak)
   {
      double sh=r[2].high, sl=r[2].low;
      for(int i=3;i<=8;i++) { sh=MathMax(sh,r[i].high); sl=MathMin(sl,r[i].low); }
      if(is_long && r[0].close>sh) trigger=true;
      if(!is_long && r[0].close<sl) trigger=true;
   }
   if(!trigger && InpEntryUseEngulfing)
   {
      bool engulf = is_long ? (r[0].close>r[1].high && r[0].open<r[1].low) : (r[0].close<r[1].low && r[0].open>r[1].high);
      if(engulf) trigger=true;
   }
   return trigger;
}

bool ComputeSL(const string symbol, const bool is_long, double &sl)
{
   double atr=GetATR(symbol,PERIOD_M15,14);
   if(atr<=0) return false;
   double bid=SymbolInfoDouble(symbol,SYMBOL_BID), ask=SymbolInfoDouble(symbol,SYMBOL_ASK);

   if(!InpUseStructureSL)
   {
      sl = is_long ? ask-InpATRSLMult*atr : bid+InpATRSLMult*atr;
      return true;
   }

   MqlRates r[];
   if(CopyRates(symbol,PERIOD_M15,0,InpSwingLookback+2,r)<InpSwingLookback+2) return false;
   ArraySetAsSeries(r,true);
   double swing = is_long ? r[1].low : r[1].high;
   for(int i=2;i<=InpSwingLookback;i++)
      swing = is_long ? MathMin(swing,r[i].low) : MathMax(swing,r[i].high);
   sl = is_long ? (swing - InpSLBufferATR*atr) : (swing + InpSLBufferATR*atr);
   return true;
}

bool DirectionAllowed(const StrengthCandidate &c)
{
   return g_regime.TrendQuality(c.symbol, PERIOD_H1, c.is_long, 200, InpADXMin, InpUseADXFilter);
}

void EvaluateEntries()
{
   if(TimeCurrent() < g_global_cooldown_until) return;
   if(g_daily_realized <= -AccountInfoDouble(ACCOUNT_BALANCE)*(InpDailyLossLimitPct/100.0)) return;
   if(g_last_entry_m15_bar==g_last_m15) return; // single entry globally per M15 bar

   for(int i=0;i<ArraySize(g_candidates);i++)
   {
      StrengthCandidate c=g_candidates[i];
      int sidx=StateIndex(c.symbol);
      if(g_states[sidx].state==ST_COOLDOWN && TimeCurrent()<g_states[sidx].cooldown_until) continue;
      if(g_trade_mng.HasPosition(c.symbol)) { g_states[sidx].state=ST_IN_TRADE; continue; }
      g_states[sidx].state=ST_SETUP;

      if(g_risk.OpenTradesTotal()>=InpMaxOpenTradesTotal) { g_log.LogEvent("ENTRY_VETO",c.symbol,"max_open_total"); continue; }
      if(g_risk.OpenTradesForCurrency(c.strong_ccy)>=InpMaxOpenTradesPerCurrency) { g_log.LogEvent("ENTRY_VETO",c.symbol,"strong_ccy_exposure"); continue; }
      if(g_risk.OpenTradesForCurrency(c.weak_ccy)>=InpMaxOpenTradesPerCurrency) { g_log.LogEvent("ENTRY_VETO",c.symbol,"weak_ccy_exposure"); continue; }
      if(g_risk.TotalOpenRiskPct()>=InpMaxTotalRiskOpenPct) { g_log.LogEvent("ENTRY_VETO",c.symbol,"max_total_risk"); continue; }
      if(!DirectionAllowed(c)) { g_log.LogEvent("ENTRY_VETO",c.symbol,"trend_filter_fail"); continue; }
      if(!PullbackAndTrigger(c.symbol,c.is_long)) { g_log.LogEvent("ENTRY_VETO",c.symbol,"pullback_or_trigger_fail"); continue; }

      double sl=0;
      if(!ComputeSL(c.symbol,c.is_long,sl)) { g_log.LogEvent("ENTRY_VETO",c.symbol,"sl_fail"); continue; }
      double entry = c.is_long ? SymbolInfoDouble(c.symbol,SYMBOL_ASK) : SymbolInfoDouble(c.symbol,SYMBOL_BID);
      double sl_dist = MathAbs(entry-sl);
      double lots = g_risk.LotSizeByRisk(c.symbol,InpRiskPerTradePct,sl_dist);
      if(lots<=0) { g_log.LogEvent("ENTRY_VETO",c.symbol,"lot_zero"); continue; }

      string cm = StringFormat("STE v0.1 div=%.2f slope=%.3f",c.divergence,c.divergence_slope);
      if(g_trade_mng.OpenPosition(c.symbol,c.is_long,lots,sl,cm))
      {
         g_states[sidx].state=ST_IN_TRADE;
         g_last_entry_m15_bar=g_last_m15;
         g_log.LogEvent("ENTRY",c.symbol,StringFormat("dir=%s lots=%.2f sl=%.5f",c.is_long?"LONG":"SHORT",lots,sl));
         break;
      }
   }
}

int OnInit()
{
   g_strength.Init(InpSymbolsForStrength);
   g_regime.Init();
   g_risk.Init(InpMagicNumber);
   g_trade_mng.Init(InpMagicNumber);
   g_log.Init("StrengthTrendEA", InpEnableCSVLog);
   EventSetTimer(10);

   if(InpSelfCheckMode)
   {
      Print("[SelfCheck] Symbols list: ", InpSymbolsForStrength);
      PrintFormat("[SelfCheck] %s point=%.5f tick_value=%.5f", _Symbol, SymbolInfoDouble(_Symbol,SYMBOL_POINT), SymbolInfoDouble(_Symbol,SYMBOL_TRADE_TICK_VALUE));
      double nominal_lot=g_risk.LotSizeByRisk(_Symbol,InpRiskPerTradePct,100*SymbolInfoDouble(_Symbol,SYMBOL_POINT));
      PrintFormat("[SelfCheck] Nominal lot for 100 points SL=%.2f", nominal_lot);
   }
   return(INIT_SUCCEEDED);
}

void OnDeinit(const int reason)
{
   EventKillTimer();
}

void OnTimer()
{
   MqlDateTime now_struct;
   TimeToStruct(TimeCurrent(),now_struct);
   int doy=now_struct.day_of_year;
   if(doy!=g_day_of_year) { g_day_of_year=doy; g_daily_realized=0.0; }

   if(IsNewBar(_Symbol,PERIOD_H1,g_last_h1))
   {
      g_strength.Compute(InpStrengthLookback,InpStrengthATRPeriod,InpStrengthTF_H4,InpStrengthTF_H1,InpStrengthTF_M15,InpWeightH4,InpWeightH1,InpWeightM15);
      g_regime.Update(_Symbol,PERIOD_H1,InpRegimeEMA,InpRegimeSlopeBars,14,InpSlopeThresh,InpAtrRatioThresh,g_strength.Dispersion(),InpDispersionThresh);
      g_strength.BuildCandidates(g_candidates,InpMaxCandidates,InpMaxSpreadPoints,0.8,PERIOD_H1,InpStrengthTF_H1,InpDivSlopeBars,InpStrengthLookback,InpStrengthATRPeriod,InpDivMin,InpDivSlopeMin,InpDivSlopeWeight);
      g_log.LogEvent("REGIME",_Symbol,g_regime.StateString()+" disp="+DoubleToString(g_strength.Dispersion(),2));
      if(InpSelfCheckMode) Print("[SelfCheck] Strength ranking: ",g_strength.RankingString()," regime=",g_regime.StateString());
   }

   if(IsNewBar(_Symbol,PERIOD_M15,g_last_m15) && g_regime.State()==REGIME_TREND)
      EvaluateEntries();
}

void OnTick()
{
   double div = (ArraySize(g_candidates)>0)?g_candidates[0].divergence:0.0;
   g_trade_mng.Manage(PERIOD_M15,14,InpBETriggerR,0.0,InpTP1R,InpTP1ClosePct,InpTrailingMode,InpTrailATRk,
                      InpTimeStopBars,InpTimeStopMinR,InpExitOnRegimeFlip,(g_regime.State()==REGIME_TREND),InpMinDivForHold,div);
}

void OnTradeTransaction(const MqlTradeTransaction& trans,
                        const MqlTradeRequest& request,
                        const MqlTradeResult& result)
{
   if(trans.type!=TRADE_TRANSACTION_DEAL_ADD) return;
   if(!HistoryDealSelect(trans.deal)) return;
   if((long)HistoryDealGetInteger(trans.deal,DEAL_MAGIC)!=InpMagicNumber) return;

   long entry_type = HistoryDealGetInteger(trans.deal, DEAL_ENTRY);
   string symbol = HistoryDealGetString(trans.deal, DEAL_SYMBOL);
   double profit = HistoryDealGetDouble(trans.deal, DEAL_PROFIT) + HistoryDealGetDouble(trans.deal, DEAL_SWAP) + HistoryDealGetDouble(trans.deal, DEAL_COMMISSION);

   if(entry_type==DEAL_ENTRY_OUT)
   {
      g_daily_realized += profit;
      if(profit<0) g_consecutive_losses++; else g_consecutive_losses=0;
      if(g_consecutive_losses>=InpCooldownAfterLosses)
         g_global_cooldown_until = TimeCurrent() + InpCooldownMinutes*60;

      int sidx=StateIndex(symbol);
      g_states[sidx].state=ST_COOLDOWN;
      g_states[sidx].cooldown_until=TimeCurrent()+InpCooldownMinutes*60;

      string reason = profit<0 ? "SL" : "Trail";
      g_log.LogTradeClose(symbol, profit>=0?"WIN":"LOSS", TimeCurrent()-60, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                          g_regime.StateString(), (ArraySize(g_candidates)>0?g_candidates[0].divergence:0), 0, 0, reason);
   }
}
