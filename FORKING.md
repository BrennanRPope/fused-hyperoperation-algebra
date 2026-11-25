# Forking this project

You are welcome to fork Fused Hyperoperation Algebra, reuse parts of it, or build on it. This page explains what that
involves.

**This is a plain-language summary, not legal advice.** The license texts themselves, [LICENSE](LICENSE) (Apache 2.0)
and [LICENSE-docs](LICENSE-docs) (CC BY 4.0), are what apply. Where this page and a license differ, the license wins.

## Two licenses, depending on the file

- **Code** (every `.lean` file, the Python checks and scripts, and the page-building code): Apache License 2.0.
- **Everything else** (write-ups, READMEs, notes, pages, figures, generated data): Creative Commons Attribution 4.0
  International (CC BY 4.0).

## If you reuse or fork the code (Apache 2.0)

- **Keep the license and notices.** Give recipients a copy of the license, and keep the copyright, patent, trademark
  and attribution notices in the source (section 4(a) and 4(c)). Keep the [NOTICE](NOTICE) file's attribution in what
  you distribute (section 4(d)).
- **Say what you changed.** Modified files must carry prominent notices stating that you changed them (section 4(b)).
- **Patent grant, and how it ends.** Contributors grant you a patent license for their contributions (section 3). It
  ends for you if you start patent litigation claiming the work infringes a patent (also section 3).
- **No trademark rights.** The license gives no permission to use the project's or author's names, except to describe
  where the work came from and to reproduce the NOTICE file (section 6).
- **No warranty.** The code comes as is (sections 7 and 8).

## If you reuse or fork the write-ups, figures or data (CC BY 4.0)

- **Give credit.** Name the creator (Brennan R. Pope), keep the copyright and license notices and the disclaimer, and
  link to the source where you reasonably can (section 3(a)(1)).
- **Say if you changed it,** and keep any earlier notes of changes (section 3(a)(1)(B)).
- **No endorsement implied.** Don't suggest that the author endorses or sponsors your version (section 2(a)(6)).
- **No added restrictions.** Don't add terms or technical locks that stop others using the material under CC BY
  (section 2(a)(5)(B)).
- If the author asks you to remove the credit from your version, you must do so where reasonably practicable
  (section 3(a)(3)).

## How to cite

Use [CITATION.cff](CITATION.cff); GitHub's "Cite this repository" button reads it. The concept DOI
[10.5281/zenodo.24539305](https://doi.org/10.5281/zenodo.24539305) cites all versions, and Zenodo gives each release
(such as v0.1.0) its own version DOI. Cite the specific version you used when you can, because names and notation
change at 1.0 (see [NAMING-COLLISIONS.md](NAMING-COLLISIONS.md)).

## Name and attribution

- The project name "Fused Hyperoperation Algebra" and the authorship credit "Brennan R. Pope" stay with the original
  project.
- Don't present a fork as the original project, or as endorsed by its author. Say near the top of your README that it
  is a fork, and link back.
- If your fork changes definitions, notation or conclusions, give it a different name, so readers can tell the two
  apart. Small forks made to send changes back don't need a new name.

## Staying compatible with the original

If you want your work to merge back cleanly, or to stay readable alongside this project, keep its conventions:

- **The design principle:** build on published work first, cite it, and label anything new
  ([CONTRIBUTING.md](CONTRIBUTING.md)).
- **The naming-collision log:** use the names in [NAMING-COLLISIONS.md](NAMING-COLLISIONS.md), and add a row for any
  new clash.
- **Assumptions:** every write-up has an Assumptions section and a Known vs new section next to it.
- **Lean rules:** no `sorry`, standard axioms only (propext, Classical.choice, Quot.sound). `lean-library/build.sh` checks the core library, the CI workflow runs it and every topic's checks on each push, and
  the Mathlib build workflow checks the mathlib edition. [ASSUMPTIONS.md](ASSUMPTIONS.md) lists which axioms each theorem
  uses.
- **Versioning:** after 1.0, changing a definition, notation or theorem name means a new major version. A fork that
  changes these is no longer compatible with that major version.

## Sending work back

Improvements are welcome upstream. Open a pull request against the original repository and fill in the template:
what existing work it builds on, what is new, assumptions, naming collisions checked, and pictures. See
[CONTRIBUTING.md](CONTRIBUTING.md). Code you contribute comes under Apache 2.0 unless you say otherwise (section 5), and the project asks that
documentation you contribute comes under CC BY 4.0.
