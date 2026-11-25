# Data used in the econ-forecasting page

- sp500_shiller_monthly.csv: Robert Shiller's monthly S&P composite series (price, dividends,
  earnings, CPI, 10-year Treasury yield, CAPE), 1871 onward, as packaged by the Frictionless
  "datasets/s-and-p-500" repository on GitHub (raw.githubusercontent.com). The rows used end
  January 2023 (later rows have no CPI or dividends).
- DJ_const.rda, SP500_const.rda, SP500.rda: daily adjusted closing prices of the Dow Jones and
  S&P 500 constituents (as of the package's 2016 snapshot), 1962 to 2015, from the R package
  qrmdata (McNeil, Frey, Embrechts, Quantitative Risk Management), GitHub mirror of CRAN.
  The constituent lists are survivors, so early years lean toward firms that later did well.

FRED, the Ken French data library and Shiller's own site were refused by this session's
network, so these mirrors were used instead.
