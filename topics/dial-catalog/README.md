# Dial catalog: the law groups every dial survey sorts (DRAFT)

**DRAFT, Nov 26 2026.** Shared by the dial surveys asked for by Brennan ("I'd like similar surveys for the other dials as
well", after the shape-dial survey).

- [`catalog.py`](catalog.py): 65 law groups. Groups 1–49 are the shape-dial survey's list, copied unchanged; groups 50–65
  cover the Lean module families that list left out. `python3 catalog.py` checks that every module family in
  `lean-library/theorems.tsv` is named by some group (it prints `[]` when none is missing).
- [`page_build.py`](page_build.py): builds a survey's self-contained page from its README.

The surveys that use it: [`homothety-dial-survey/`](../homothety-dial-survey/) (k),
[`curvature-dial-survey/`](../curvature-dial-survey/) (d), [`order-dial-survey/`](../order-dial-survey/) (q),
[`odds-dial-survey/`](../odds-dial-survey/) (t) and
[`ratio-dial-survey/`](../ratio-dial-survey/) (r). The roundness dial p has its own,
[`shape-dial-survey/`](../shape-dial-survey/). The cubic form C has none: Brennan settled that it is derived from k and d,
not a free dial (see [`dial-map/`](../dial-map/)).

## Assumptions

- A group's Lean modules are named by family; a group can hold laws from several modules, and the sort is a judgement.

## Known vs new

Nothing here is a result: it is a list.
