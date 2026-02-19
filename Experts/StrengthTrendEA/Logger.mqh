#ifndef __STRENGTH_TREND_LOGGER_MQH__
#define __STRENGTH_TREND_LOGGER_MQH__

class CLogger
{
private:
   string m_trade_file;
   string m_event_file;
   bool   m_enabled;

   void WriteLine(const string file_name, const string line)
   {
      int h = FileOpen(file_name, FILE_READ|FILE_WRITE|FILE_CSV|FILE_ANSI|FILE_SHARE_READ|FILE_SHARE_WRITE);
      if(h == INVALID_HANDLE)
      {
         PrintFormat("[Logger] Failed opening %s, err=%d", file_name, GetLastError());
         return;
      }
      FileSeek(h, 0, SEEK_END);
      FileWriteString(h, line + "\n");
      FileClose(h);
   }

public:
   bool Init(const string prefix, const bool enabled=true)
   {
      m_enabled = enabled;
      m_trade_file = prefix + "_trades.csv";
      m_event_file = prefix + "_events.csv";

      int h1 = FileOpen(m_trade_file, FILE_READ|FILE_WRITE|FILE_CSV|FILE_ANSI|FILE_SHARE_READ|FILE_SHARE_WRITE);
      if(h1 == INVALID_HANDLE) return false;
      if(FileSize(h1) == 0)
         FileWriteString(h1, "timestamp_close,symbol,direction,entry_time,entry_price,exit_price,lots,initial_SL_price,initial_risk_money,R_result,MAE_R,MFE_R,duration_minutes,regime_at_entry,divergence_at_entry,dist_at_entry,atr_ratio_at_entry,reason_exit\n");
      FileClose(h1);

      int h2 = FileOpen(m_event_file, FILE_READ|FILE_WRITE|FILE_CSV|FILE_ANSI|FILE_SHARE_READ|FILE_SHARE_WRITE);
      if(h2 == INVALID_HANDLE) return false;
      if(FileSize(h2) == 0)
         FileWriteString(h2, "timestamp,event_type,symbol,details\n");
      FileClose(h2);

      return true;
   }

   void LogEvent(const string event_type, const string symbol, const string details)
   {
      if(!m_enabled) return;
      string ts = TimeToString(TimeCurrent(), TIME_DATE|TIME_MINUTES|TIME_SECONDS);
      string line = ts + "," + event_type + "," + symbol + "," + details;
      WriteLine(m_event_file, line);
   }

   void LogTradeClose(const string symbol,
                      const string direction,
                      const datetime entry_time,
                      const double entry_price,
                      const double exit_price,
                      const double lots,
                      const double initial_sl,
                      const double initial_risk_money,
                      const double r_result,
                      const double mae_r,
                      const double mfe_r,
                      const int duration_minutes,
                      const string regime_at_entry,
                      const double divergence_at_entry,
                      const double dist_at_entry,
                      const double atr_ratio_at_entry,
                      const string reason_exit)
   {
      if(!m_enabled) return;
      string ts = TimeToString(TimeCurrent(), TIME_DATE|TIME_MINUTES|TIME_SECONDS);
      string line = StringFormat("%s,%s,%s,%s,%.5f,%.5f,%.2f,%.5f,%.2f,%.3f,%.3f,%.3f,%d,%s,%.3f,%.3f,%.3f,%s",
                                 ts,
                                 symbol,
                                 direction,
                                 TimeToString(entry_time, TIME_DATE|TIME_MINUTES),
                                 entry_price,
                                 exit_price,
                                 lots,
                                 initial_sl,
                                 initial_risk_money,
                                 r_result,
                                 mae_r,
                                 mfe_r,
                                 duration_minutes,
                                 regime_at_entry,
                                 divergence_at_entry,
                                 dist_at_entry,
                                 atr_ratio_at_entry,
                                 reason_exit);
      WriteLine(m_trade_file, line);
   }
};

#endif
