# Data and generated files left out

Kept out of git to keep the repo small: Python caches, `.npy` and `.npz` arrays, run logs, Lean build output (`lean-library/build/`, `.olean`), and any file over 1 MB. The scripts in each topic regenerate them, and `lean-library/build.sh` rebuilds the Lean output. Downloaded data can be fetched again from the sources in each topic's `data/SOURCES.md`.

Largest files left out:
- `topics/dark-sector/data/pantheon/PantheonPlusSH0ES_STAT+SYS.cov` (33 MB, Pantheon+ covariance; see dark-sector/data/SOURCES.md)
- `topics/econ-forecasting/data/SP500_const.rda` (4 MB)
- `topics/ns-rungs/*.npz` (shell-model runs)
- `topics/navier-stokes/*.npy`, `*.npz`, `*.log`
- `topics/primes/zeta_zeros_1000.npy`, `wheel_predict.npy`
