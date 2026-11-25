# Data used in the actuarial page

All three files are R data files copied from the GitHub mirror of CRAN (github.com/cran).

- danish.rda, danishuni.rda: Danish fire insurance losses, 1980 to 1990. 2,167 claims of
  at least 1 million Danish kroner, adjusted for inflation to 1985 values (Copenhagen Re;
  studied by McNeil 1997, "Estimating the tails of loss severity distributions using extreme
  value theory"). From the R packages evir and fitdistrplus. danishuni has the claim dates
  (days since 1 January 1970).
- demoUsa.rda: US Social Security Area period life tables (l_x, survivors out of 100,000)
  for 1990, 2000 and 2007, male and female. From the R package lifecontingencies.
- soa08Act.rda: the Society of Actuaries Illustrative Life Table (l_x, interest 6%), the
  teaching table from the SOA exam syllabus. From the R package lifecontingencies.

Read them with the Python package rdata (python3 -m pip install rdata).
