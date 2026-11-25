"""Shared loaders for the econ-forecasting checks. Data sources: data/SOURCES.md."""
import os, warnings
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

def shiller_monthly():
    df = pd.read_csv(os.path.join(DATA, "sp500_shiller_monthly.csv"), parse_dates=["Date"])
    df = df.rename(columns={"SP500": "P", "Dividend": "D", "Earnings": "E",
                            "Consumer Price Index": "CPI", "Long Interest Rate": "GS10"})
    df = df[(df.CPI > 0) & (df.GS10 > 0) & (df.D > 0)].set_index("Date")
    return df

def shiller_annual():
    """January-to-January annual series.
    infl: CPI inflation over the year; r_nom: 10-year yield in January (decimal);
    R_real: real total return of the S&P over the year (price change + dividends, deflated)."""
    m = shiller_monthly()
    a = m[m.index.month == 1].copy()
    a.index = a.index.year
    out = pd.DataFrame(index=a.index)
    out["infl"] = a.CPI.shift(-1) / a.CPI - 1
    out["r_nom"] = a.GS10 / 100
    # dividends: annual rate D accrues through the year; use the average of monthly D/12
    dsum = m.D.groupby(m.index.year).mean()
    out["R_nom"] = (a.P.shift(-1) + dsum.reindex(a.index)) / a.P - 1
    out["R_real"] = (1 + out.R_nom) / (1 + out.infl) - 1
    out["pe10"] = m.PE10.groupby(m.index.year).first().reindex(a.index)
    out["real_price"] = a.P / a.CPI
    return out.dropna(subset=["infl"])

def xts(name):
    """Read an xts object from a qrmdata .rda file: returns a DataFrame with dates."""
    import rdata
    warnings.filterwarnings("ignore")
    path = os.path.join(DATA, name + ".rda")
    conv = rdata.read_rda(path)
    arr = conv[name]
    p = rdata.parser.parse_file(path)
    a = p.object.value[0].attributes
    idx = None
    while a is not None and a.value:
        try: tag = a.tag.value.value
        except Exception: tag = None
        if tag == b"index": idx = np.asarray(a.value[0].value); break
        a = a.value[1] if len(a.value) > 1 else None
    dates = pd.to_datetime(idx, unit="s").normalize()
    return pd.DataFrame(np.asarray(arr.values), index=dates, columns=[str(c) for c in arr.coords["dim_1"].values])
