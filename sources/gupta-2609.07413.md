# DRAFT — Gupta 2026: the 2PN two-body Lagrangian for Einstein–Maxwell–dilaton binaries (transcribed reference)

**Status: DRAFT. Transcribed reference, plus one clearly marked section that is not in the paper.**
Nothing in the paper's part of this file is our own result.

## Citation

- Pawan Kumar Gupta, "Binary Black Holes in Einstein-Maxwell-Dilaton Theory: Second
  Post-Newtonian Dynamics from Effective Field Theory", arXiv:2609.07413 [gr-qc], version 1,
  submitted 7 Sep 2026 (25 pages, 8 figures). No journal reference is listed on arXiv as of
  Nov 27, 2026.
- Read from the arXiv TeX source (`main.tex`), downloaded from
  `https://arxiv.org/e-print/2609.07413` on Nov 27, 2026. SHA-256 of the downloaded archive:
  `2129264cd72f5c9723270f673a063e924fda83c578c184b3ab0b43b4b75303ec`.
- The paper and its TeX are **not** in this repository; copyright stays with the author. This
  file records the formulas (coefficients, signs, conventions, equation numbers) in a layout of
  its own, not the paper's text.

## How this file was made, and what was checked

- The three Lagrangian equations, (4.2), (5.9) and (7.1), and the test-body energy (7.8) were
  read from the TeX by a script into a computer-algebra system (sympy). Every table below was
  printed from that parsed form. No coefficient was typed by hand.
- Because the tables are regenerated, the terms are **expanded and regrouped**: one row per
  kinematic structure, with its whole coefficient. The paper groups some terms differently
  (for example it writes a common factor $\Lambda_1[v_1^2 - 2(v_1 v_2)]$). The values are the
  same.
- Equation numbers are those printed in the arXiv PDF (version 1).
- Checks run on the parsed form:
  1. With all charges and scalar couplings set to zero, (5.9) plus its $1\leftrightarrow2$ image
     equals the Einstein–Infeld–Hoffmann 1PN Lagrangian exactly.
  2. The static test-body limit of (4.2)+(5.9)+(7.1) equals minus the paper's Eq. (7.8), term
     by term, at orders $1/r$, $1/r^2$ and $1/r^3$. This is the paper's own statement (7.10).
     It confirms that the static terms were read correctly and that the $1\leftrightarrow2$
     rule below is applied the way the paper applies it.
- **Not checked:** the velocity- and acceleration-dependent terms of (7.1) have no independent
  check here beyond the reading itself. Nothing was compared with a journal version (none
  exists yet), and the paper's diagram-by-diagram results (its Sec. 6) were not read into the
  script or summed.

## Conventions (the paper's)

- Units: $c = 1$, except where powers of $1/c^2$ are shown to mark the PN order. Gaussian-type
  charge units: the Coulomb term is $-q_1 q_2/r$ with no $4\pi\epsilon_0$.
- Metric signature $(-,+,+,+)$. Einstein frame. Harmonic gauge for gravity, Feynman gauge for
  the vector field. Riemann convention, Eq. (1.1):
  $R^\mu{}_{\nu\alpha\beta} = \Gamma^\mu{}_{\nu\beta,\alpha} - \Gamma^\mu{}_{\nu\alpha,\beta} + \Gamma^\rho{}_{\nu\beta}\Gamma^\mu{}_{\rho\alpha} - \Gamma^\rho{}_{\nu\alpha}\Gamma^\mu{}_{\rho\beta}$.
- Bodies $A = 1, 2$ with positions $\mathbf{x}_A$, velocities $\mathbf{v}_A$, accelerations
  $\mathbf{a}_A$, masses $m_A$, charges $q_A$.
- $\mathbf{r} \equiv \mathbf{x}_1 - \mathbf{x}_2$, $r = |\mathbf{r}|$,
  $\mathbf{n} \equiv \mathbf{r}/r$. So $\mathbf{n}$ points **from body 2 to body 1**.
- "$+(1\leftrightarrow2)$" means: add the same expression with the body labels exchanged
  **and** $\mathbf{n} \to -\mathbf{n}$.
- The Lagrangian is the EFT (Kaluza–Klein potential-mode) one. At 2PN it **depends on the
  accelerations** $\mathbf{a}_1, \mathbf{a}_2$; they have not been removed.
- Sign of $L$: $L = T - V$, so the Newtonian term is $+G_{12} m_1 m_2/r$.

In the tables, $(v_1 v_2) = \mathbf{v}_1\cdot\mathbf{v}_2$, $(n v_1) = \mathbf{n}\cdot\mathbf{v}_1$,
$(v_1 a_2) = \mathbf{v}_1\cdot\mathbf{a}_2$, $v_1^2 = \mathbf{v}_1\cdot\mathbf{v}_1$, and so on.

## The theory and the body parameters

Action, Eqs. (2.1)–(2.5), with $F_{\mu\nu} = \nabla_\mu A_\nu - \nabla_\nu A_\mu$ and
$\Gamma^\mu \equiv \Gamma^\mu_{\rho\sigma} g^{\rho\sigma}$:

$$S = S_g + S_{em} + S_\varphi + S_m \qquad (2.1)$$

$$S_g = \frac{1}{16\pi G}\int d^4x \sqrt{-g}\, R - \frac{1}{32\pi G}\int d^4x \sqrt{-g}\, g_{\mu\nu}\Gamma^\mu\Gamma^\nu \qquad (2.2)$$

$$S_{em} = -\frac{1}{16\pi}\int d^4x \sqrt{-g}\; e^{-2a\varphi}\left(F^2 + 2(\nabla_\mu A^\mu)^2\right) \qquad (2.3)$$

$$S_\varphi = -\frac{1}{16\pi G}\int d^4x \sqrt{-g}\,\left(2 g^{\mu\nu}\partial_\mu\varphi\,\partial_\nu\varphi\right) \qquad (2.4)$$

$$S_m = -\sum_A \int dt \left[\mathfrak{m}_A(\varphi)\sqrt{-g_{\mu\nu} v_A^\mu v_A^\nu} - q_A A_\mu v_A^\mu\right] \qquad (2.5)$$

$a$ is the dilaton coupling constant. Each body has a field-dependent mass $\mathfrak{m}_A(\varphi)$,
with, Eq. (2.6),

$$\alpha(\varphi) \equiv \frac{d\ln\mathfrak{m}}{d\varphi}, \qquad \beta(\varphi) \equiv \frac{d\alpha}{d\varphi}, \qquad \beta'(\varphi) \equiv \frac{d\beta}{d\varphi},$$

and, Eq. (2.7), with $m \equiv \mathfrak{m}(\varphi_0)$, $\delta\varphi = \varphi - \varphi_0$ and the
background value set to $\varphi_0 = 0$:

$$\mathfrak{m}(\varphi) = m\left[1 + \alpha\,\delta\varphi + \tfrac12(\alpha^2+\beta)\,\delta\varphi^2 + \tfrac{1}{3!}\left[\alpha(\alpha^2+3\beta)+\beta'\right]\delta\varphi^3 + O(\delta\varphi^4)\right].$$

$\alpha_A, \beta_A, \beta'_A$ below are the values at $\varphi_0$. For black holes with $a = 1$,
Eq. (2.8) gives

$$\alpha = \frac{q^2}{2Gm^2}, \qquad \beta = \frac{q^2}{Gm^2} - \frac{q^4}{2G^2m^4}, \qquad \beta' = \frac{2q^2}{Gm^2} - \frac{3q^4}{G^2m^4} + \frac{q^6}{G^3m^6}.$$

### Combinations used in the Lagrangian

| Symbol | Definition | Eq. |
|---|---|---|
| $G_{12}$ | $G(1 + \alpha_1\alpha_2)$ | after (4.2) |
| $\bar\gamma_{12}$ | $-\dfrac{2\alpha_1\alpha_2}{1+\alpha_1\alpha_2}$ | (5.10) |
| $\bar\beta_1$ | $\dfrac12\,\dfrac{\beta_1\alpha_2^2}{(1+\alpha_1\alpha_2)^2}$ | (5.11) |
| $\hat m_A$ | $m_A(1 + a\alpha_A)$ | (5.12) |
| $\Lambda_1$ | $7 + 10\bar\gamma_{12} + \tfrac72\bar\gamma_{12}^2 - 2\delta_2 + 2\bar\beta_1$ | (7.2) |
| $\Lambda_2$ | $8 + 10\bar\gamma_{12} + \tfrac72\bar\gamma_{12}^2 - 2\delta_2 + 4\bar\beta_1$ | (7.3) |
| $\delta_1$ | $\dfrac{\alpha_1^2}{(1+\alpha_1\alpha_2)^2}$ | (7.4) |
| $\epsilon_1$ | $\dfrac{\beta'_1\alpha_2^3}{(1+\alpha_1\alpha_2)^3}$ | (7.4) |
| $\zeta_{12}$ | $\dfrac{\beta_1\beta_2\alpha_1\alpha_2}{(1+\alpha_1\alpha_2)^3}$ | (7.5) |

Points to read carefully:

- The paper says $\bar\beta_2$ is $\bar\beta_1$ with $1\leftrightarrow2$. It defines only
  $\delta_1$ and $\epsilon_1$, yet (7.1) uses $\delta_2$ and (7.8) uses $\epsilon_2$. This file
  takes $\delta_2 = \alpha_2^2/(1+\alpha_1\alpha_2)^2$ and
  $\epsilon_2 = \beta'_2\alpha_1^3/(1+\alpha_1\alpha_2)^3$, by the same exchange. Check 2 above
  passes with this reading but would pass with any reading that is consistent between (7.1)
  and (7.8), so it does not test it.
- Note the index pattern: $\bar\beta_1$ and $\epsilon_1$ carry body 1's $\beta_1$, $\beta'_1$
  but body **2**'s $\alpha_2$; $\delta_1$ carries body 1's $\alpha_1$.
- $\Lambda_2$ is **not** $\Lambda_1$ with the labels exchanged. Both are built from $\delta_2$
  and $\bar\beta_1$. When forming the $1\leftrightarrow2$ image, this file first writes
  $\Lambda_1$, $\Lambda_2$ out with (7.2), (7.3) and then exchanges labels, so their images
  contain $\delta_1$ and $\bar\beta_2$. The paper does not spell this step out.
- $G_{12}$, $\bar\gamma_{12}$ and $\zeta_{12}$ are symmetric under $1\leftrightarrow2$.
- Some terms of (5.9) and (7.1) carry the bare $G$ and others $G_{12}$. That is how they are
  printed; the tables keep the distinction.

## The Lagrangian

$$L = L_{0\rm PN} + L_{1\rm PN} + L_{2\rm PN}.$$

Each of (5.9) and (7.1) is printed as a "half" $\ell$ followed by $+(1\leftrightarrow2)$:
$L_{n\rm PN} = \ell_n + \ell_n\big|_{1\leftrightarrow2,\ \mathbf{n}\to-\mathbf{n}}$. The tables give
the **half** $\ell_n$, exactly the part the paper prints. A term that is symmetric under the
exchange therefore appears in $L$ with twice the coefficient shown.

Under the exchange: $m_1\leftrightarrow m_2$, $q_1\leftrightarrow q_2$,
$\alpha_1\leftrightarrow\alpha_2$, $\beta_1\leftrightarrow\beta_2$,
$\bar\beta_1\leftrightarrow\bar\beta_2$, $\delta_1\leftrightarrow\delta_2$,
$\epsilon_1\leftrightarrow\epsilon_2$, $v_1^2\leftrightarrow v_2^2$,
$(n v_1)\leftrightarrow -(n v_2)$, $(n a_1)\leftrightarrow -(n a_2)$,
$(v_1 a_2)\leftrightarrow(a_1 v_2)$, $(v_1 a_1)\leftrightarrow(v_2 a_2)$; $(v_1 v_2)$ and
$(a_1 a_2)$ are unchanged.

### 0PN, Eq. (4.2)

$$L_{0\rm PN} = \frac{G_{12} m_1 m_2}{r} - \frac{q_1 q_2}{r}.$$

This is written in full; there is no $+(1\leftrightarrow2)$. The paper does not print the
Newtonian kinetic term $\tfrac12 m_1 v_1^2 + \tfrac12 m_2 v_2^2$ in (4.2). The $v^4$ and $v^6$
kinetic terms in (5.9) and (7.1) are those of $-m\sqrt{1-v^2}$, whose $v^2$ term is
$\tfrac12 m v^2$.

### 1PN, Eq. (5.9): the half $\ell_1$

| Structure | Power of $r$ | Coefficient |
|---|---|---|
| $v_1^4$ | $r^0$ | $\frac{m_1}{8}$ |
| $v_1^2$ | $r^{-1}$ | $\frac{G_{12} m_1 m_2 \left(2 \bar\gamma_{12} + 3\right)}{2}$ |
| $(v_1 v_2)$ | $r^{-1}$ | $\frac{- 4 G_{12} \bar\gamma_{12} m_1 m_2 - 7 G_{12} m_1 m_2 + q_1 q_2}{4}$ |
| $(n v_1) (n v_2)$ | $r^{-1}$ | $\frac{- G_{12} m_1 m_2 + q_1 q_2}{4}$ |
| $1$ | $r^{-2}$ | $\frac{m_1 \left(2 G a \alpha_1 q_1 q_2 - G a \alpha_1 q_2^{2} + 2 G q_1 q_2 - G q_2^{2} - 2 G_{12}^{2} \bar\beta_2 m_1 m_2 - G_{12}^{2} m_1 m_2\right)}{2}$ |

5 rows. Here $\hat m_1 = m_1(1+a\alpha_1)$ has been written out.

### 2PN, Eq. (7.1): the half $\ell_2$

$\Lambda_1$ and $\Lambda_2$ are left as symbols, as in the paper.

**Velocity terms (no acceleration).**

| Structure | Power of $r$ | Coefficient |
|---|---|---|
| $v_1^6$ | $r^0$ | $\frac{m_1}{16}$ |
| $v_1^4$ | $r^{-1}$ | $\frac{G_{12} m_1 m_2 \left(4 \bar\gamma_{12} + 7\right)}{8}$ |
| $v_1^2 v_2^2$ | $r^{-1}$ | $- \frac{- 3 G_{12} m_1 m_2 + q_1 q_2}{16}$ |
| $v_1^2 (v_1 v_2)$ | $r^{-1}$ | $- \frac{G_{12} m_1 m_2 \left(2 \bar\gamma_{12} + 5\right)}{4}$ |
| $v_1^2 (n v_1) (n v_2)$ | $r^{-1}$ | $- \frac{G_{12} m_1 m_2 \left(2 \bar\gamma_{12} + 3\right)}{4}$ |
| $v_2^2 (n v_1)^2$ | $r^{-1}$ | $\frac{- G_{12} m_1 m_2 + q_1 q_2}{8}$ |
| $(v_1 v_2)^2$ | $r^{-1}$ | $\frac{G_{12} m_1 m_2 + q_1 q_2}{8}$ |
| $(v_1 v_2) (n v_1) (n v_2)$ | $r^{-1}$ | $\frac{G_{12} m_1 m_2 \left(2 \bar\gamma_{12} + 3\right)}{4}$ |
| $(n v_1)^2 (n v_2)^2$ | $r^{-1}$ | $- \frac{3 \left(- G_{12} m_1 m_2 + q_1 q_2\right)}{16}$ |
| $v_1^2$ | $r^{-2}$ | $- \frac{m_1 \left(2 G a \alpha_1 q_1 q_2 - G a \alpha_1 q_2^{2} + 2 G q_1 q_2 + G q_2^{2} - G_{12}^{2} \Lambda_1 m_2^{2}\right)}{4}$ |
| $v_2^2$ | $r^{-2}$ | $\frac{m_1 \left(2 G a \alpha_1 q_2^{2} + G_{12}^{2} \Lambda_2 m_2^{2}\right)}{4}$ |
| $(v_1 v_2)$ | $r^{-2}$ | $- \frac{m_1 \left(G a \alpha_1 q_2^{2} - G q_2^{2} + G_{12}^{2} \Lambda_1 m_2^{2}\right)}{2}$ |
| $(n v_1)^2$ | $r^{-2}$ | $\frac{m_1 \left(8 G a \alpha_1 q_1 q_2 + 8 G q_1 q_2 - 4 G q_2^{2} + 4 G_{12}^{2} \delta_2 m_2^{2} + G_{12}^{2} \bar\gamma_{12}^{2} m_2^{2} + 4 G_{12}^{2} \bar\gamma_{12} m_2^{2} + 4 G_{12}^{2} m_2^{2}\right)}{8}$ |
| $(n v_1) (n v_2)$ | $r^{-2}$ | $- \frac{m_1 \left(8 G a \alpha_1 q_1 q_2 - 4 G a \alpha_1 q_2^{2} + 8 G q_1 q_2 - 8 G q_2^{2} - 8 G_{12}^{2} \bar\beta_1 m_2^{2} + 4 G_{12}^{2} \delta_2 m_2^{2} + G_{12}^{2} \bar\gamma_{12}^{2} m_2^{2} + 4 G_{12}^{2} \bar\gamma_{12} m_2^{2}\right)}{4}$ |
| $(n v_2)^2$ | $r^{-2}$ | $- \frac{m_1 \left(4 G a \alpha_1 q_2^{2} + 8 G q_2^{2} + 8 G_{12}^{2} \bar\beta_1 m_2^{2} - 4 G_{12}^{2} \delta_2 m_2^{2} - G_{12}^{2} \bar\gamma_{12}^{2} m_2^{2} - 4 G_{12}^{2} \bar\gamma_{12} m_2^{2}\right)}{8}$ |

**Acceleration terms.**

| Structure | Power of $r$ | Coefficient |
|---|---|---|
| $(n a_1) (n a_2)$ | $r$ | $\frac{- G_{12} m_1 m_2 + q_1 q_2}{16}$ |
| $(a_1 a_2)$ | $r$ | $- \frac{- 8 G_{12} \bar\gamma_{12} m_1 m_2 - 15 G_{12} m_1 m_2 + 3 q_1 q_2}{16}$ |
| $v_1^2 (n a_2)$ | $r^0$ | $- \frac{G_{12} m_1 m_2}{16}$ |
| $v_2^2 (n a_1)$ | $r^0$ | $- \frac{- G_{12} m_1 m_2 + 2 q_1 q_2}{16}$ |
| $(n v_1)^2 (n a_2)$ | $r^0$ | $\frac{G_{12} m_1 m_2}{16}$ |
| $(n v_1) (v_1 a_2)$ | $r^0$ | $- \frac{- 4 G_{12} \bar\gamma_{12} m_1 m_2 - 7 G_{12} m_1 m_2 + 2 q_1 q_2}{8}$ |
| $(n v_1) (v_2 a_2)$ | $r^0$ | $- \frac{G_{12} m_1 m_2 \left(2 \bar\gamma_{12} + 3\right)}{4}$ |
| $(n v_2)^2 (n a_1)$ | $r^0$ | $\frac{- G_{12} m_1 m_2 + 2 q_1 q_2}{16}$ |
| $(n v_2) (a_1 v_2)$ | $r^0$ | $- \frac{G_{12} m_1 m_2 \left(4 \bar\gamma_{12} + 7\right)}{8}$ |
| $(n v_2) (v_1 a_1)$ | $r^0$ | $\frac{G_{12} m_1 m_2 \left(2 \bar\gamma_{12} + 3\right)}{4}$ |

**Static terms.**

| Structure | Power of $r$ | Coefficient |
|---|---|---|
| $1$ | $r^{-3}$ | $G^{2} m_1^{2} \left(- \frac{q_2 \left(2 a^{2} \alpha_1^{2} q_1 - 3 a^{2} \alpha_1^{2} q_2 + 4 a \alpha_1 q_1 - 6 a \alpha_1 q_2 + \alpha_1^{2} q_1 + 3 q_1 - 3 q_2\right)}{3}\right)$ |
| $1$ | $r^{-3}$ | $G^{2} m_1 m_2 \frac{2 a^{2} \alpha_1 \alpha_2 q_1^{2} + 30 a^{2} \alpha_1 \alpha_2 q_1 q_2 - 6 a \alpha_1^{2} \alpha_2 q_1 q_2 + 3 a \alpha_1^{2} \alpha_2 q_2^{2} + 2 a \alpha_1 q_1^{2} - 18 a \alpha_1 q_1 q_2 + 5 a \alpha_1 q_2^{2} - 6 a \alpha_2 \beta_1 q_1 q_2 + 3 a \alpha_2 \beta_1 q_2^{2} + \alpha_1 \alpha_2 q_1^{2} - 18 \alpha_1 \alpha_2 q_1 q_2 + 3 q_1^{2} - 60 q_1 q_2}{6}$ |
| $1$ | $r^{-3}$ | $G \left(- \frac{q_1^{2} q_2 \left(2 a^{2} q_1 - 3 a^{2} q_2 - 9 q_2\right)}{6}\right)$ |
| $1$ | $r^{-3}$ | $G_{12}^{3} m_1^{2} m_2^{2} \frac{4 \bar\beta_1 + \bar\gamma_{12} + \zeta_{12} + 3}{2}$ |
| $1$ | $r^{-3}$ | $G_{12}^{3} m_1 m_2^{3} \frac{12 \bar\beta_1 + 4 \delta_2 + 2 \epsilon_1 + \bar\gamma_{12}^{2} + 4 \bar\gamma_{12} + 6}{12}$ |

15 + 10 + 5 rows.

### The same halves in machine-readable form

Printed from the same parsed expressions (sympy `sstr`). Names: `G12` $=G_{12}$, `gb` $=\bar\gamma_{12}$, `bb1` $=\bar\beta_1$, `d2` $=\delta_2$, `e1` $=\epsilon_1$, `z12` $=\zeta_{12}$, `L1`, `L2` $=\Lambda_{1,2}$, `al1` $=\alpha_1$, `be1` $=\beta_1$ (the un-barred $\beta_1$), `v1v1` $=v_1^2$, `v1v2` $=(v_1 v_2)$, `nv1` $=(n v_1)$, `v1a2` $=(v_1 a_2)$, `a1v2` $=(a_1 v_2)$, `v1a1` $=(v_1 a_1)$, `v2a2` $=(v_2 a_2)$, `na1` $=(n a_1)$, `a1a2` $=(a_1 a_2)$.

```python
ell_1PN = G*m1*(a*al1 + 1)*(q1*q2 - q2**2/2)/r**2 - G12**2*m1**2*m2*(2*bb2 + 1)/(2*r**2) + G12*m1*m2*(gb*(v1v1 - v1v2) - nv1*nv2/4 + 3*v1v1/2 - 7*v1v2/4)/r + m1*v1v1**2/8 + q1*q2*(nv1*nv2 + v1v2)/(4*r)

ell_2PN = G**2*m1**2*(-2*a**2*al1**2*q1*q2/3 + a**2*al1**2*q2**2 - 4*a*al1*q1*q2/3 + 2*a*al1*q2**2 - al1**2*q1*q2/3 - q1*q2 + q2**2)/r**3 + G**2*m1*m2*(a**2*al1*al2*q1**2/3 + 5*a**2*al1*al2*q1*q2 - a*al1**2*al2*q1*q2 + a*al1**2*al2*q2**2/2 + a*al1*q1**2/3 - 3*a*al1*q1*q2 + 5*a*al1*q2**2/6 - a*al2*be1*q1*q2 + a*al2*be1*q2**2/2 + al1*al2*q1**2/6 - 3*al1*al2*q1*q2 + q1**2/2 - 10*q1*q2)/r**3 + G*m1*(a*al1*q2**2*v2v2/2 + nv1**2*(q1*q2*(a*al1 + 1) - q2**2/2) + nv1*nv2*(-2*q1*q2*(a*al1 + 1) + q2**2*(a*al1 + 2)) - nv2**2*q2**2*(a*al1/2 + 1) - q1*q2*v1v1*(a*al1 + 1)/2 - q2**2*v1v1*(-a*al1 + 1)/4 + q2**2*v1v2*(-a*al1 + 1)/2)/r**2 + G*(-a**2*q1**3*q2/3 + a**2*q1**2*q2**2/2 + 3*q1**2*q2**2/2)/r**3 + G12**3*m1**2*m2**2*(4*bb1 + gb + z12 + 3)/(2*r**3) + G12**3*m1*m2**3*(2*bb1 + 2*d2/3 + e1/3 + gb**2/6 + 2*gb/3 + 1)/(2*r**3) + G12**2*m1*m2**2*(L1*(v1v1 - 2*v1v2) + L2*v2v2 + nv1**2*(2*d2 + gb**2/2 + 2*gb + 2) - nv1*nv2*(-8*bb1 + 4*d2 + gb**2 + 4*gb) + nv2**2*(-4*bb1 + 2*d2 + gb**2/2 + 2*gb))/(4*r**2) + G12*m1*m2*(3*nv1**2*nv2**2 - 2*nv1**2*v2v2 - nv1*nv2*v1v1*(8*gb + 12) + nv1*nv2*v1v2*(8*gb + 12) + r**2*(a1a2*(8*gb + 15) - na1*na2) + r*(-na1*nv2**2 + na1*v2v2 + na2*nv1**2 - na2*v1v1 + (8*gb + 12)*(-nv1*v2a2 + nv2*v1a1) + (8*gb + 14)*(-a1v2*nv2 + nv1*v1a2)) + v1v1**2*(8*gb + 14) - v1v1*v1v2*(8*gb + 20) + 3*v1v1*v2v2 + 2*v1v2**2)/(16*r) + m1*v1v1**3/16 + q1*q2*(-3*nv1**2*nv2**2 + 2*nv1**2*v2v2 + r**2*(-3*a1a2 + na1*na2) + r*(2*na1*nv2**2 - 2*na1*v2v2 - 4*nv1*v1a2) - v1v1*v2v2 + 2*v1v2**2)/(16*r)
```

## The paper's static test-body check, Eqs. (7.6)–(7.10)

Body 1 is an exact EMd black hole, body 2 a point particle in its field, in harmonic
coordinates. The exact static energy, Eq. (7.6), is
$V_{\rm TB}(r) = c^2\,\mathfrak{m}_2[\varphi_1(r)]\sqrt{-g^{(1)}_{00}(r)} - q_2 A^{(1)}_0(r)$.
The source black hole obeys, Eq. (7.7),

$$\frac{q_1^2}{G m_1^2} = \frac{2\alpha_1}{a} - \frac{1-a^2}{a^2}\,\alpha_1^2.$$

Expanded to static 2PN order, Eq. (7.8) gives $V_{\rm TB} - m_2 c^2$ as the terms below
(written here with $c = 1$; in the paper the $1/r^2$ terms carry $1/c^2$ and the $1/r^3$ terms
$1/c^4$, with remainder $O(c^{-6})$):

| Structure | Power of $r$ | Coefficient |
|---|---|---|
| $1$ | $r^{-1}$ | $- G_{12} m_1 m_2 + q_1 q_2$ |
| $1$ | $r^{-2}$ | $\frac{- 2 G a \alpha_1 m_1 q_1 q_2 + G a \alpha_2 m_2 q_1^{2} - 2 G m_1 q_1 q_2 + G m_2 q_1^{2} + 2 G_{12}^{2} \bar\beta_2 m_1^{2} m_2 + G_{12}^{2} m_1^{2} m_2}{2}$ |
| $1$ | $r^{-3}$ | $G^{2} m_1^{2} \frac{q_1 q_2 \left(2 a^{2} \alpha_1^{2} + 4 a \alpha_1 + \alpha_1^{2} + 3\right)}{3}$ |
| $1$ | $r^{-3}$ | $- \frac{q_1^{2} \left(2 a^{2} \alpha_1 \alpha_2 + 3 a \alpha_1 \alpha_2^{2} + 3 a \alpha_1 \beta_2 + 2 a \alpha_1 + 5 a \alpha_2 + \alpha_1 \alpha_2 + 3\right)}{6} G^{2} m_1 m_2$ |
| $1$ | $r^{-3}$ | $G \frac{a^{2} q_1^{3} q_2}{3}$ |
| $1$ | $r^{-3}$ | $G_{12}^{3} m_1^{3} m_2 \left(- \frac{12 \bar\beta_2 + 4 \delta_1 + 2 \epsilon_2 + \bar\gamma_{12}^{2} + 4 \bar\gamma_{12} + 6}{12}\right)$ |

The static interaction Lagrangian is $L^{\rm stat}_{\rm TB,int} = -(V_{\rm TB} - m_2 c^2)$,
Eq. (7.9). The test-body limit of the two-body result is taken by $m_2 \to \varepsilon m_2$,
$q_2 \to \varepsilon q_2$, keeping the terms linear in $\varepsilon$ with $q_2/m_2$,
$\alpha_2$, $\beta_2$, $\beta'_2$ fixed, at $\mathbf{v}_A = \mathbf{a}_A = 0$, and including the
$1\leftrightarrow2$ image. Eq. (7.10) states that this equals $L^{\rm stat}_{\rm TB,int}$ up to
$O(c^{-6})$. The script confirms it on the parsed forms (check 2 above).

## Limits the paper states

- $\alpha_A = \beta_A = \beta'_A = 0$ and $a = 0$: the 2PN Einstein–Maxwell Lagrangian of the
  paper's Ref. [Gupta:2022spq].
- $q_A = 0$: the 2PN scalar-tensor Lagrangian of its Ref. [Almeida:2024uph].
- Both: the 2PN general-relativity Lagrangian of Gilmore and Ross, its Ref. [Gilmore:2008gq].
- $G \to 0$: the 2PN electromagnetic two-body dynamics of its Ref. [Barker:1980kef].
- Agreement with a calculation using the EFTofPNG package, its Ref. [Levi:2017kzq].

These are the paper's statements. Only the 1PN general-relativity limit and the static
test-body limit were re-checked here.

# NOT IN THE PAPER — equal-mass reductions (our own, DRAFT)

**Everything from here to the "Known vs new" section was derived for this file. The paper has
no equal-mass reduction, no centre-of-mass reduction and no symmetric-mass-ratio form.** The
reductions below are plain substitutions into the parsed Lagrangian, with the
$1\leftrightarrow2$ image added, done by the same script. They are **full** Lagrangians, not
halves. They inherit any misreading of the paper and the interpretation of $\Lambda_{1,2}$,
$\delta_2$, $\epsilon_2$ described above. They have not been checked against any other source.

Not done: removing the accelerations (order reduction), passing to the centre-of-mass frame,
or reducing to circular orbits. A test that needs those must do them itself.

## Reduction A: general relativity, equal masses

Set $q_A = 0$, $\alpha_A = \beta_A = \beta'_A = 0$, $a = 0$ (so $G_{12} = G$ and
$\bar\gamma_{12} = \bar\beta_A = \delta_A = \epsilon_A = \zeta_{12} = 0$, $\Lambda_1 = 7$,
$\Lambda_2 = 8$) and $m_1 = m_2 = m$.

$$L_{0\rm PN} = \frac{G m^2}{r} \quad (\text{plus } \tfrac12 m v_1^2 + \tfrac12 m v_2^2, \text{ not printed in the paper}).$$

**Full $L_{1\rm PN}$:**

| Structure | Power of $r$ | Coefficient |
|---|---|---|
| $v_1^4$ | $r^0$ | $\frac{m}{8}$ |
| $v_2^4$ | $r^0$ | $\frac{m}{8}$ |
| $v_1^2$ | $r^{-1}$ | $\frac{3 G m^{2}}{2}$ |
| $v_2^2$ | $r^{-1}$ | $\frac{3 G m^{2}}{2}$ |
| $(v_1 v_2)$ | $r^{-1}$ | $- \frac{7 G m^{2}}{2}$ |
| $(n v_1) (n v_2)$ | $r^{-1}$ | $- \frac{G m^{2}}{2}$ |
| $1$ | $r^{-2}$ | $- G^{2} m^{3}$ |

**Full $L_{2\rm PN}$:**

| Structure | Power of $r$ | Coefficient |
|---|---|---|
| $(n a_1) (n a_2)$ | $r$ | $- \frac{G m^{2}}{8}$ |
| $(a_1 a_2)$ | $r$ | $\frac{15 G m^{2}}{8}$ |
| $v_1^6$ | $r^0$ | $\frac{m}{16}$ |
| $v_1^2 (n a_2)$ | $r^0$ | $- \frac{G m^{2}}{8}$ |
| $v_2^6$ | $r^0$ | $\frac{m}{16}$ |
| $v_2^2 (n a_1)$ | $r^0$ | $\frac{G m^{2}}{8}$ |
| $(n v_1)^2 (n a_2)$ | $r^0$ | $\frac{G m^{2}}{8}$ |
| $(n v_1) (v_1 a_2)$ | $r^0$ | $\frac{7 G m^{2}}{4}$ |
| $(n v_1) (v_2 a_2)$ | $r^0$ | $- \frac{3 G m^{2}}{2}$ |
| $(n v_2)^2 (n a_1)$ | $r^0$ | $- \frac{G m^{2}}{8}$ |
| $(n v_2) (a_1 v_2)$ | $r^0$ | $- \frac{7 G m^{2}}{4}$ |
| $(n v_2) (v_1 a_1)$ | $r^0$ | $\frac{3 G m^{2}}{2}$ |
| $v_1^4$ | $r^{-1}$ | $\frac{7 G m^{2}}{8}$ |
| $v_1^2 v_2^2$ | $r^{-1}$ | $\frac{3 G m^{2}}{8}$ |
| $v_1^2 (v_1 v_2)$ | $r^{-1}$ | $- \frac{5 G m^{2}}{4}$ |
| $v_1^2 (n v_1) (n v_2)$ | $r^{-1}$ | $- \frac{3 G m^{2}}{4}$ |
| $v_1^2 (n v_2)^2$ | $r^{-1}$ | $- \frac{G m^{2}}{8}$ |
| $v_2^4$ | $r^{-1}$ | $\frac{7 G m^{2}}{8}$ |
| $v_2^2 (v_1 v_2)$ | $r^{-1}$ | $- \frac{5 G m^{2}}{4}$ |
| $v_2^2 (n v_1)^2$ | $r^{-1}$ | $- \frac{G m^{2}}{8}$ |
| $v_2^2 (n v_1) (n v_2)$ | $r^{-1}$ | $- \frac{3 G m^{2}}{4}$ |
| $(v_1 v_2)^2$ | $r^{-1}$ | $\frac{G m^{2}}{4}$ |
| $(v_1 v_2) (n v_1) (n v_2)$ | $r^{-1}$ | $\frac{3 G m^{2}}{2}$ |
| $(n v_1)^2 (n v_2)^2$ | $r^{-1}$ | $\frac{3 G m^{2}}{8}$ |
| $v_1^2$ | $r^{-2}$ | $\frac{15 G^{2} m^{3}}{4}$ |
| $v_2^2$ | $r^{-2}$ | $\frac{15 G^{2} m^{3}}{4}$ |
| $(v_1 v_2)$ | $r^{-2}$ | $- 7 G^{2} m^{3}$ |
| $(n v_1)^2$ | $r^{-2}$ | $\frac{G^{2} m^{3}}{2}$ |
| $(n v_2)^2$ | $r^{-2}$ | $\frac{G^{2} m^{3}}{2}$ |
| $1$ | $r^{-3}$ | $4 G^{3} m^{4}$ |

For unequal masses in general relativity, the static part of the full $L_{2\rm PN}$ is
$\frac{G^{3} m_1 m_2 \left(m_1^{2} + 6 m_1 m_2 + m_2^{2}\right)}{2 r^{3}}$.

## Reduction B: two identical Einstein–Maxwell–dilaton bodies

Set $m_1 = m_2 = m$, $q_1 = q_2 = q$, $\alpha_1 = \alpha_2 = \alpha$,
$\beta_1 = \beta_2 = \beta$, $\beta'_1 = \beta'_2 = \beta'$. The combinations become

$$G_{12} = G(1+\alpha^2), \quad \bar\gamma_{12} = -\frac{2\alpha^2}{1+\alpha^2}, \quad \bar\beta \equiv \bar\beta_1 = \bar\beta_2 = \frac{\beta\alpha^2}{2(1+\alpha^2)^2},$$

$$\delta \equiv \delta_1 = \delta_2 = \frac{\alpha^2}{(1+\alpha^2)^2}, \quad \epsilon \equiv \epsilon_1 = \epsilon_2 = \frac{\beta'\alpha^3}{(1+\alpha^2)^3}, \quad \zeta_{12} = \frac{\beta^2\alpha^2}{(1+\alpha^2)^3}.$$

$\Lambda_1$ and $\Lambda_2$ are written out with (7.2), (7.3). The coefficients below keep
$G_{12}$, $\bar\gamma_{12}$, $\bar\beta$, $\delta$, $\epsilon$, $\zeta_{12}$ as symbols.

$$L_{0\rm PN} = \frac{G_{12} m^2 - q^2}{r}.$$

**Full $L_{1\rm PN}$:**

| Structure | Power of $r$ | Coefficient |
|---|---|---|
| $v_1^4$ | $r^0$ | $\frac{m}{8}$ |
| $v_2^4$ | $r^0$ | $\frac{m}{8}$ |
| $v_1^2$ | $r^{-1}$ | $\frac{G_{12} m^{2} \left(2 \bar\gamma_{12} + 3\right)}{2}$ |
| $v_2^2$ | $r^{-1}$ | $\frac{G_{12} m^{2} \left(2 \bar\gamma_{12} + 3\right)}{2}$ |
| $(v_1 v_2)$ | $r^{-1}$ | $\frac{- 4 G_{12} \bar\gamma_{12} m^{2} - 7 G_{12} m^{2} + q^{2}}{2}$ |
| $(n v_1) (n v_2)$ | $r^{-1}$ | $\frac{- G_{12} m^{2} + q^{2}}{2}$ |
| $1$ | $r^{-2}$ | $m \left(G a \alpha q^{2} + G q^{2} - 2 G_{12}^{2} \bar\beta m^{2} - G_{12}^{2} m^{2}\right)$ |

**Full $L_{2\rm PN}$, velocity terms:**

| Structure | Power of $r$ | Coefficient |
|---|---|---|
| $v_1^6$ | $r^0$ | $\frac{m}{16}$ |
| $v_2^6$ | $r^0$ | $\frac{m}{16}$ |
| $v_1^4$ | $r^{-1}$ | $\frac{G_{12} m^{2} \left(4 \bar\gamma_{12} + 7\right)}{8}$ |
| $v_1^2 v_2^2$ | $r^{-1}$ | $- \frac{- 3 G_{12} m^{2} + q^{2}}{8}$ |
| $v_1^2 (v_1 v_2)$ | $r^{-1}$ | $- \frac{G_{12} m^{2} \left(2 \bar\gamma_{12} + 5\right)}{4}$ |
| $v_1^2 (n v_1) (n v_2)$ | $r^{-1}$ | $- \frac{G_{12} m^{2} \left(2 \bar\gamma_{12} + 3\right)}{4}$ |
| $v_1^2 (n v_2)^2$ | $r^{-1}$ | $\frac{- G_{12} m^{2} + q^{2}}{8}$ |
| $v_2^4$ | $r^{-1}$ | $\frac{G_{12} m^{2} \left(4 \bar\gamma_{12} + 7\right)}{8}$ |
| $v_2^2 (v_1 v_2)$ | $r^{-1}$ | $- \frac{G_{12} m^{2} \left(2 \bar\gamma_{12} + 5\right)}{4}$ |
| $v_2^2 (n v_1)^2$ | $r^{-1}$ | $\frac{- G_{12} m^{2} + q^{2}}{8}$ |
| $v_2^2 (n v_1) (n v_2)$ | $r^{-1}$ | $- \frac{G_{12} m^{2} \left(2 \bar\gamma_{12} + 3\right)}{4}$ |
| $(v_1 v_2)^2$ | $r^{-1}$ | $\frac{G_{12} m^{2} + q^{2}}{4}$ |
| $(v_1 v_2) (n v_1) (n v_2)$ | $r^{-1}$ | $\frac{G_{12} m^{2} \left(2 \bar\gamma_{12} + 3\right)}{2}$ |
| $(n v_1)^2 (n v_2)^2$ | $r^{-1}$ | $- \frac{3 \left(- G_{12} m^{2} + q^{2}\right)}{8}$ |
| $v_1^2$ | $r^{-2}$ | $\frac{m \left(G a \alpha q^{2} - 3 G q^{2} + 6 G_{12}^{2} \bar\beta m^{2} - 4 G_{12}^{2} \delta m^{2} + 7 G_{12}^{2} \bar\gamma_{12}^{2} m^{2} + 20 G_{12}^{2} \bar\gamma_{12} m^{2} + 15 G_{12}^{2} m^{2}\right)}{4}$ |
| $v_2^2$ | $r^{-2}$ | $\frac{m \left(G a \alpha q^{2} - 3 G q^{2} + 6 G_{12}^{2} \bar\beta m^{2} - 4 G_{12}^{2} \delta m^{2} + 7 G_{12}^{2} \bar\gamma_{12}^{2} m^{2} + 20 G_{12}^{2} \bar\gamma_{12} m^{2} + 15 G_{12}^{2} m^{2}\right)}{4}$ |
| $(v_1 v_2)$ | $r^{-2}$ | $- \frac{m \left(2 G a \alpha q^{2} - 2 G q^{2} + 4 G_{12}^{2} \bar\beta m^{2} - 4 G_{12}^{2} \delta m^{2} + 7 G_{12}^{2} \bar\gamma_{12}^{2} m^{2} + 20 G_{12}^{2} \bar\gamma_{12} m^{2} + 14 G_{12}^{2} m^{2}\right)}{2}$ |
| $(n v_1)^2$ | $r^{-2}$ | $\frac{m \left(2 G a \alpha q^{2} - 2 G q^{2} - 4 G_{12}^{2} \bar\beta m^{2} + 4 G_{12}^{2} \delta m^{2} + G_{12}^{2} \bar\gamma_{12}^{2} m^{2} + 4 G_{12}^{2} \bar\gamma_{12} m^{2} + 2 G_{12}^{2} m^{2}\right)}{4}$ |
| $(n v_1) (n v_2)$ | $r^{-2}$ | $- \frac{m \left(4 G a \alpha q^{2} - 8 G_{12}^{2} \bar\beta m^{2} + 4 G_{12}^{2} \delta m^{2} + G_{12}^{2} \bar\gamma_{12}^{2} m^{2} + 4 G_{12}^{2} \bar\gamma_{12} m^{2}\right)}{2}$ |
| $(n v_2)^2$ | $r^{-2}$ | $\frac{m \left(2 G a \alpha q^{2} - 2 G q^{2} - 4 G_{12}^{2} \bar\beta m^{2} + 4 G_{12}^{2} \delta m^{2} + G_{12}^{2} \bar\gamma_{12}^{2} m^{2} + 4 G_{12}^{2} \bar\gamma_{12} m^{2} + 2 G_{12}^{2} m^{2}\right)}{4}$ |

**Full $L_{2\rm PN}$, acceleration terms:**

| Structure | Power of $r$ | Coefficient |
|---|---|---|
| $(n a_1) (n a_2)$ | $r$ | $\frac{- G_{12} m^{2} + q^{2}}{8}$ |
| $(a_1 a_2)$ | $r$ | $- \frac{- 8 G_{12} \bar\gamma_{12} m^{2} - 15 G_{12} m^{2} + 3 q^{2}}{8}$ |
| $v_1^2 (n a_2)$ | $r^0$ | $\frac{- G_{12} m^{2} + q^{2}}{8}$ |
| $v_2^2 (n a_1)$ | $r^0$ | $- \frac{- G_{12} m^{2} + q^{2}}{8}$ |
| $(n v_1)^2 (n a_2)$ | $r^0$ | $- \frac{- G_{12} m^{2} + q^{2}}{8}$ |
| $(n v_1) (v_1 a_2)$ | $r^0$ | $- \frac{- 4 G_{12} \bar\gamma_{12} m^{2} - 7 G_{12} m^{2} + q^{2}}{4}$ |
| $(n v_1) (v_2 a_2)$ | $r^0$ | $- \frac{G_{12} m^{2} \left(2 \bar\gamma_{12} + 3\right)}{2}$ |
| $(n v_2)^2 (n a_1)$ | $r^0$ | $\frac{- G_{12} m^{2} + q^{2}}{8}$ |
| $(n v_2) (a_1 v_2)$ | $r^0$ | $\frac{- 4 G_{12} \bar\gamma_{12} m^{2} - 7 G_{12} m^{2} + q^{2}}{4}$ |
| $(n v_2) (v_1 a_1)$ | $r^0$ | $\frac{G_{12} m^{2} \left(2 \bar\gamma_{12} + 3\right)}{2}$ |

**Full $L_{2\rm PN}$, static terms:**

| Structure | Power of $r$ | Coefficient |
|---|---|---|
| $1$ | $r^{-3}$ | $G^{2} m^{2} \frac{q^{2} \left(34 a^{2} \alpha^{2} - 3 a \alpha^{3} - 3 a \alpha \beta - 7 a \alpha - 19 \alpha^{2} - 57\right)}{3}$ |
| $1$ | $r^{-3}$ | $G \frac{q^{4} \left(a^{2} + 9\right)}{3}$ |
| $1$ | $r^{-3}$ | $G_{12}^{3} m^{4} \frac{36 \bar\beta + 4 \delta + 2 \epsilon + \bar\gamma_{12}^{2} + 10 \bar\gamma_{12} + 6 \zeta_{12} + 24}{6}$ |

## Reduction C: equal masses only

With $m_1 = m_2 = m$ but different charges or scalar couplings, nothing simplifies beyond
replacing $m_1, m_2$ by $m$ in the half tables and adding the $1\leftrightarrow2$ image. It is
not tabulated.

# Known vs new

- **Known (published, not ours):** everything above the line "NOT IN THE PAPER". It is the
  content of arXiv:2609.07413, an unrefereed preprint as of Nov 27, 2026.
- **Ours, DRAFT:** the equal-mass reductions A and B, and the reading of $\delta_2$,
  $\epsilon_2$ and of the $1\leftrightarrow2$ image of $\Lambda_{1,2}$.

# Assumptions

- Of the paper: Einstein–Maxwell–dilaton theory in the Einstein frame; non-spinning,
  monopolar point bodies with field-dependent mass (no finite-size terms); conservative
  near-zone (potential-mode) dynamics only, through 2PN; black-hole values of
  $\alpha, \beta, \beta'$ are known in closed form only for $a = 1$, Eq. (2.8).
- Of this file: the TeX source matches what the author intended; the script's reading of the
  equations is right (supported by the two checks listed, which do not cover the 2PN velocity
  and acceleration terms); the interpretation of the exchange rule stated above.
