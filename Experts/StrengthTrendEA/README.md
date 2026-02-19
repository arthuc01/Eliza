# StrengthTrendEA (Version 0.1)

A production-oriented MT5 Expert Advisor implementing a **strong-vs-weak currency strength trend-following strategy** with M15 pullback entries.

## Strategy overview
- **Strength engine** (H4/H1/M15 weighted) computes currency scores for USD, EUR, GBP, JPY, CHF, CAD, AUD, NZD.
- Builds top-2 vs bottom-2 combinations and maps to tradable symbols.
- Candidate ranking uses divergence level plus divergence slope (strengthening over recent H1 bars) to avoid late/stale entries.
- **Regime detector** allows entries only when market is in TREND state, based on:
  - EMA(50) slope normalized by ATR
  - ATR ratio (ATR / ATR-SMA)
  - Strength dispersion
  - 2-bar persistence to enable/disable
- **Entry timeframe M15**:
  - Pullback toward EMA20 or 3+ countertrend candles
  - Trigger by swing break and/or engulfing rejection
  - Anti-chase, spread, and ADR exhaustion filters
  - Single new entry maximum per M15 bar (global throttle)
- **Risk and portfolio controls**:
  - Fixed-fractional lot sizing from SL distance
  - Max open trades total/per currency
  - Max aggregate open risk
  - Daily loss lockout and cooldown after loss streak
- **Trade management**:
  - Initial SL at order placement (structure or ATR)
  - BE move at +1R
  - TP1 partial at +2R
  - Runner with ATR trail (swing trail hook reserved)
  - Time-stop for stale trades
  - Optional close on regime flip or divergence collapse
- **CSV logs** in terminal `MQL5/Files`:
  - `StrengthTrendEA_trades.csv`
  - `StrengthTrendEA_events.csv`

## Files
- `StrengthTrendEA.mq5` — main orchestration and state machine.
- `StrengthEngine.mqh` — multi-timeframe currency strength + candidate ranking.
- `RegimeDetector.mqh` — trend/range classification + trend quality filter.
- `RiskManager.mqh` — lot sizing + exposure calculations.
- `TradeManager.mqh` — in-trade management logic.
- `Logger.mqh` — structured CSV logging.
- `params_default.set` — starter parameter profile.

## Key inputs
- `InpSymbolsForStrength`: comma-separated symbols used to compute currency strength (edit here to match broker symbols).
- `InpWeightH4/H1/M15`: weighted blend for strength score.
- `InpDivMin`, `InpDivSlopeBars`, `InpDivSlopeMin`, `InpDivSlopeWeight`: minimum divergence + strengthening filter and ranking boost for rising divergence.
- `InpSlopeThresh`, `InpAtrRatioThresh`, `InpDispersionThresh`: regime thresholds.
- `InpDistMin/InpDistMax`, `InpADRPctMax`, `InpMaxSpreadPoints`: anti-chase and execution filters.
- `InpUseStructureSL`, `InpATRSLMult`: SL mode.
- `InpRiskPerTradePct`, `InpMaxOpenTradesTotal` (default 5), `InpMaxOpenTradesPerCurrency`, `InpMaxTotalRiskOpenPct`.
- `InpDailyLossLimitPct`, `InpCooldownAfterLosses`, `InpCooldownMinutes`.

## How to run
1. Copy folder into `MQL5/Experts/StrengthTrendEA`.
2. Compile `StrengthTrendEA.mq5` in MetaEditor.
3. Attach EA to an M15 chart of a liquid major (e.g., EURUSD).
4. Ensure required H1/H4 history is downloaded.
5. Keep `InpSelfCheckMode=true` for first launch verification.

## Backtest recommendations
- Tester mode: **Every tick based on real ticks**.
- Chart timeframe: **M15**.
- Ensure H1 and H4 data availability.
- Test across multiple years and symbols, then walk-forward.

## Limitations and next steps
- Broker suffix/prefix handling may require adapting symbol parser.
- Current close-log detail uses conservative placeholders for some fields; can be upgraded with per-ticket persistence for full MAE/MFE and exact exit reasons.
- Swing-based trailing mode can be expanded further.
- Add optimization profiles with robust out-of-sample validation.
