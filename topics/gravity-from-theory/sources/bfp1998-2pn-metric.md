# Blanchet, Faye & Ponsot 1998: the 2PN (2.5PN) metric of two point masses in harmonic coordinates

**Transcribed reference. Nothing in this file is our own result.** It is a verbatim copy of
equations from the published literature, kept here so that sessions which cannot reach arXiv
can read them.

## Citation

- L. Blanchet, G. Faye and B. Ponsot, "Gravitational field and equations of motion of compact
  binaries to 5/2 post-Newtonian order", Phys. Rev. D **58**, 124002 (1998),
  doi:10.1103/PhysRevD.58.124002, arXiv:gr-qc/9804079.
- Transcribed from the arXiv LaTeX source of **version 2** (11 Aug 1998), file
  `change.grqc.tex`, downloaded from `https://arxiv.org/e-print/gr-qc/9804079` on Nov 27, 2026.
  SHA-256 of the downloaded archive: `e4555cfa49528c0b7dbd3ca72d6337c3a0d16568f48427781503001976249f83`.
- Backup source, quoted in the last section: L. Blanchet, "Gravitational Radiation from
  Post-Newtonian Sources and Inspiralling Compact Binaries", Living Rev. Relativity **17**, 2
  (2014), arXiv:1310.1528, **version 2** (17 Feb 2014).
  SHA-256 of the downloaded archive (`https://arxiv.org/e-print/1310.1528v2`): `9d023916dea7381c144053a20ad8cdc677fb7d4e0e5e91ae7efad8661bc3ed61`.

The papers themselves are not in this repository.

## How this was transcribed, and what was checked

- Every LaTeX block below was cut from the source file by line range with a script. Nothing was
  retyped, reordered or corrected. Only trailing spaces were removed.
- The LaTeX is the paper's own 1998 REVTeX (`eqnarray`, `mathletters`, `{a \over b}`), so the
  blocks are shown as code, not rendered.
- Equation numbers are the ones printed in the arXiv PDF. In the source they are the `\label`s:
  `\label{2.4b}` is Eq. (2.4b), `\label{B.1a}` is Eq. (B1a).
- **Not checked:** the transcription has not been compared with the journal version or with any
  erratum, and none of the formulas has been re-derived or tested numerically here.

## Where things are in the paper

The request that led to this file expected the two-body potentials in Section III. They are in
**Appendix B**; Section III only sets up the point-mass source and the regularization.

| What | Where in the paper |
|---|---|
| Source densities $\sigma, \sigma_i, \sigma_{ij}$ | Sec. II, Eq. (2.1) |
| Definitions of $V$, $V_i$, $\hat W_{ij}$, $\hat R_i$, $\hat X$ for any matter source | Sec. II, Eqs. (2.3), (2.4) |
| The metric in terms of the potentials | Sec. II, Eq. (2.5) |
| Harmonic-gauge identities; compact / quadratic / non-compact split | Sec. II, Eqs. (2.6)-(2.10) |
| Point-mass source, effective masses, Hadamard partie finie | Sec. III, Eqs. (3.1)-(3.9) |
| 1.5PN equations of motion used to reduce accelerations | Sec. III, Eq. (3.11) |
| **The potentials for two bodies, at any field point** | **Appendix B, Eq. (B1)** |
| **Their regularized values at body 1** | **Appendix B, Eq. (B2)** |
| Their regularized gradients at body 1 | Appendix B, Eq. (B3) |
| The explicit two-body metric at any field point | Sec. VII, Eq. (7.2) |
| The explicit metric regularized at body 1 | Sec. VII, Eq. (7.6) |

## Notation (collected from the paper)

- Body positions ${\bf y}_1(t), {\bf y}_2(t)$ and coordinate velocities ${\bf v}_1, {\bf v}_2$ in
  harmonic coordinates; $m_1, m_2$ are the constant Schwarzschild masses.
- $r_1 = |{\bf x}-{\bf y}_1|$, ${\bf n}_1 = ({\bf x}-{\bf y}_1)/r_1$, and the same with 2.
- $r_{12} = |{\bf y}_1-{\bf y}_2|$, ${\bf v}_{12} = {\bf v}_1-{\bf v}_2$, and ${\bf n}_{12}$ is
  the unit vector from body 2 to body 1. The source prints
  `${\bf n}_{12}=|{\bf y}_1-{\bf y}_2|/r_{12}$` after Eq. (3.11); the absolute-value bars there
  are a misprint for $({\bf y}_1-{\bf y}_2)/r_{12}$. It is left as printed in the block below.
- Scalar products are written with parentheses: $(n_{12}v_1) = {\bf n}_{12}\cdot{\bf v}_1$.
- $S = r_1 + r_2 + r_{12}$ (Eq. (5.6) of the paper, recalled after Eq. (7.2)).
- $O(n)$ means a term of order $\varepsilon^n$, with $\varepsilon \sim 1/c$.
- "$1 \leftrightarrow 2$" means the same terms with the labels 1 and 2 exchanged (under which
  ${\bf n}_{12}$ and ${\bf v}_{12}$ change sign).
- $(F)_1$ is the Hadamard partie finie of $F$ at body 1, Eq. (3.6).
- $n_1^{(i} n_2^{j)}$ is symmetrization over $ij$; $n_{12}^{ij}$ and $v_{12}^{ij}$ are the usual
  shorthand for $n_{12}^i n_{12}^j$ and $v_{12}^i v_{12}^j$. The paper does not spell these last
  two out in the passages quoted here.
- $\hat W_{ij} = W_{ij} - \delta_{ij} W_{kk}$, where $W_{ij}$ is the definition of the paper's
  Ref. [B95] (the paper's footnote N4).

**Orders.** The expansions below are 2.5PN, not 2PN: they include the odd, radiation-reaction
terms ($1/c^5$ and $1/c^7$ in $g_{00}$, $1/c^6$ in $g_{0i}$, $1/c^5$ in $g_{ij}$, and the
matching odd terms in the potentials). The 2PN metric is what remains when those odd terms are
dropped.

## Section II: definitions for a general matter source

Source densities, Eq. (2.1):

```latex
\begin{mathletters}
\label{2.1}
\begin{eqnarray}
 \sigma &\equiv& {T^{00}+T^{ii}\over c^2} \ , \\
 \sigma_i &\equiv& {T^{0i}\over c} \ , \\
 \sigma_{ij} &\equiv& T^{ij} \ .
\end{eqnarray}
\end{mathletters}
```

The retarded potential $V$, Eq. (2.3):

```latex
\begin{equation}
V({\bf x},t) = \Box^{-1}_R\left\{-4\pi G\sigma \right\} \equiv G
\int {d^3{\bf z}\over |{\bf x}-{\bf z}|}
\sigma( {\bf z}, t - |{\bf x} -{\bf z}|/c ) \ . \label{2.3}
\end{equation}
```

The other retarded potentials, Eqs. (2.4a)-(2.4e):

```latex
\begin{mathletters}
\label{2.4}
\begin{eqnarray}
V_i &=& \Box^{-1}_R\left\{-4\pi G \sigma_i\right\} \ ,  \label{2.4a}
\\
{\hat W}_{ij} &=& \Box^{-1}_R\left\{-4 \pi G (\sigma_{ij} - \delta_{ij}
\sigma_{kk}) - \partial_i V \partial_j V\right\} \ , \label{2.4b} \\
{\hat R}_i &=& \Box^{-1}_R\left\{ - 4\pi G (V\sigma_i - V_i \sigma) - 2
\partial_k V
\partial_i V_k - {3\over 2} \partial_t V \partial_i V \right\} \ ,
\label{2.4c}\\
{\hat X} &=& \Box^{-1}_R\biggl\{ -4\pi G V \sigma_{ii} + 2 V_i \partial_t
\partial_i V +V \partial_t^2 V  \nonumber \\
&&\quad\quad+{3\over 2} (\partial_t V)^2 - 2 \partial_i V_j \partial_j V_i +
\hat{W}_{ij} \partial^2_{ij} V \biggr\} \ . \label{2.4d}
\end{eqnarray}
In addition, we shall often consider the trace of the potential
${\hat W}_{ij}$, i.e.

\begin{equation}
{\hat W}_{ii} = \Box^{-1}_R\left\{8\pi G \sigma_{ii}-\partial_i V \partial_i
V \right\} \ . \label{2.4e}
\end{equation}
\end{mathletters}
```

The 2.5PN metric in harmonic coordinates, Eqs. (2.5a)-(2.5c):

```latex
\begin{mathletters}
\label{2.5}
\begin{eqnarray}
g_{00} &=& -1 + {2\over c^2}V - {2\over c^4}V^2 + {8\over c^6} \left[\hat{X}
+ V_iV_i + {V^3\over 6}\right] +O(8) \ , \label{2.5a}\\
g_{0i} &=& -{4\over c^3} V_i - {8\over c^5} \hat{R}_i + O(7) \ , \label{2.5b}\\
g_{ij} &=& \delta_{ij} \left( 1 + {2\over c^2} V + {2\over c^4} V^2
\right) + {4\over c^4} \hat{W}_{ij} + O(6) \ . \label{2.5c}
\end{eqnarray}
\end{mathletters}
```

The harmonic-coordinate conditions as identities between the potentials, Eq. (2.6):

```latex
\begin{mathletters}
\label{2.6}
\begin{eqnarray}
& & \partial_t\left\{V+{1\over c^2}\left[{1\over 2}\hat{W}_{ii}+
2V^2\right]\right\} +\partial_i\left\{V_i+{2\over c^2}
\left[\hat{R}_i+VV_i\right]\right\}= O(4) \ , \label{2.6a}\\
& & \partial_t V_i +\partial_j\left\{{\hat W}_{ij}-{1\over 2}\delta_{ij}
{\hat W}_{kk} \right\}=O(2) \ . \label{2.6b}
\end{eqnarray}
\end{mathletters}
```

The split into compact (C), quadratic ($\partial V\partial V$) and non-compact (NC) parts, Eqs. (2.7)-(2.10):

```latex
\begin{mathletters}
\label{2.7}
\begin{eqnarray}
\hat{W}_{ij} &=& \hat{W}_{ij}^{({\rm C})} + \hat{W}_{ij}^{(\partial V
\partial V)} \ , \label{2.7a}\\
\hat{R}_i &=& \hat{R}_i^{({\rm C})} + \hat{R}_i^{(\partial V\partial V)} \ ,
\label{2.7b}\\
\hat{X} &=& \hat{X}^{({\rm C})} + \hat{X}^{(\partial V\partial V)} +
\hat{X}^{({\rm NC})} \ . \label{2.7c}
\end{eqnarray}
\end{mathletters}
```

```latex
\begin{mathletters}
\label{2.8}
\begin{eqnarray}
\hat{W}_{ij}^{({\rm C})} &=& \Box^{-1}_R\left\{-4 \pi G (\sigma_{ij} -
\delta_{ij} \sigma_{kk})\right\} \ , \label{2.8a} \\
\hat{R}_i^{({\rm C})} &=& \Box^{-1}_R\left\{ - 4\pi G (V\sigma_i -
V_i \sigma)\right\} \ , \label{2.8b} \\
\hat{X}^{({\rm C})} &=& \Box^{-1}_R\left\{ -4\pi G V \sigma_{ii}\right\} \ .
\label{2.8c}
\end{eqnarray}
\end{mathletters}
```

```latex
\begin{mathletters}
\label{2.9}
\begin{eqnarray}
\hat{W}_{ij}^{(\partial V\partial V)} &=& \Box^{-1}_R\left\{- \partial_i V
\partial_j V\right\} \ , \label{2.9a} \\
\hat{R}_i^{(\partial V\partial V)} &=& \Box^{-1}_R\left\{ - 2
\partial_k V \partial_i V_k - {3\over 2} \partial_t V \partial_i V \right\}
\ , \label{2.9b} \\
\hat{X}^{(\partial V\partial V)} &=& \Box^{-1}_R\left\{ 2 V_i \partial_t
\partial_i V +V \partial_t^2 V  +{3\over 2}
(\partial_t V)^2 \right.\nonumber \\
&&\quad\quad -\left.2 \partial_i V_j \partial_j V_i + \hat{W}_{ij}^{({\rm C})}
\partial^2_{ij} V \right\} \ . \label{2.9c}
\end{eqnarray}
\end{mathletters}
```

```latex
\begin{equation}
\hat{X}^{({\rm NC})} = \Box^{-1}_R\left\{\hat{W}_{ij}^{(\partial V\partial V)}
\partial^2_{ij} V \right\} \ . \label{2.10}
\end{equation}
```

## Section III: two point masses and the regularization

Stress-energy tensor, Eq. (3.1):

```latex
\begin{equation}
T^{\mu\nu} ({\bf x},t) = \mu_1(t) v_1^\mu(t) v_1^\nu(t)
\delta({\bf x}-{\bf y}_1(t)) + 1\leftrightarrow 2 \ . \label{3.1}
\end{equation}
```

Effective masses, Eqs. (3.2), (3.3):

```latex
\begin{equation}
\mu_1(t) =\left( {m_1\over \sqrt{gg_{\rho\sigma} {\displaystyle
\frac{v_1^\rho v_1^\sigma}{c^2}}}} \right)_1 \ , \label{3.2}
\end{equation}
```

```latex
\begin{equation}
\tilde{\mu}_1(t) = \mu_1(t)\left[1+{v_1^2 \over c^2}\right] \ , \label{3.3}
\end{equation}
```

Source densities for two particles, Eq. (3.4):

```latex
\begin{mathletters}
\label{3.4}
\begin{eqnarray}
\sigma &=& \tilde{\mu}_1 \delta ({\bf x}-{\bf y}_1) + 1\leftrightarrow 2 \ ,
\label{3.4a} \\
\sigma_i &=& \mu_1 v_1^i \delta ({\bf x}-{\bf y}_1)+ 1\leftrightarrow 2 \ ,
\label{3.4b} \\
\sigma_{ij} &=& \mu_1 v_1^i v_1^j \delta ({\bf x}
-{\bf y}_1)+ 1\leftrightarrow 2 \ . \label{3.4c}
\end{eqnarray}
\end{mathletters}
```

Hadamard partie finie at body 1, Eqs. (3.5)-(3.7):

```latex
\begin{equation}
F({\bf x};{\bf y}_1,{\bf y}_2)=\sum_{\scriptstyle -k_0\leq k\leq 0}
r_1^k f_k({\bf n}_1;{\bf y}_1,{\bf y}_2)+O(r_1)  \label{3.5}
\end{equation}
```

```latex
\begin{equation}
(F)_1 \equiv F({\bf y}_1;{\bf y}_1,{\bf y}_2)\equiv\int {d\Omega({\bf n}_1)
\over 4\pi} f_0({\bf n}_1;{\bf y}_1,{\bf y}_2)\ . \label{3.6}
\end{equation}
```

```latex
\begin{equation}
\int d^3{\bf x} F({\bf x};{\bf y}_1,{\bf y}_2) \delta({\bf x}-{\bf y}_1)
\equiv (F)_1 \ , \label{3.7}
\end{equation}
```

Newtonian-order example, Eqs. (3.8), (3.9):

```latex
\begin{mathletters}
\label{3.8}
\begin{eqnarray}
V &=& {Gm_1\over r_1}+O(2)+1\leftrightarrow 2 \ , \label{3.8a}\\
V_i &=& {Gm_1\over r_1}v_1^i+O(2)+1\leftrightarrow 2 \ . \label{3.8b}
\end{eqnarray}
\end{mathletters}
```

```latex
\begin{mathletters}
\label{3.9}
\begin{eqnarray}
(V)_1 &=& {Gm_2\over r_{12}}+O(2) \ , \label{3.9a}\\
(V_i)_1 &=& {Gm_2\over r_{12}}v_2^i+O(2) \ ,  \label{3.9b}
\end{eqnarray}
\end{mathletters}
```

The 1.5PN equations of motion used to replace accelerations ("order reduction"), Eq. (3.11), with the sentence that follows it:

```latex
\begin{eqnarray}
{dv_1^i\over dt} =&-&{Gm_2\over r_{12}^2} n_{12}^i \left\{1+
{1\over c^2}\left[ - 5{Gm_1\over r_{12}} - 4{Gm_2\over r_{12}} +
v^2_1 + 2v^2_2 - 4(v_1v_2) -{3\over 2} (n_{12}v_2)^2 \right]\right\}\nonumber\\
 &+&{Gm_2\over c^2r_{12}^2} v^i_{12} \left[4(n_{12}v_1) -3(n_{12}v_2)\right]
+O(4) \ , \label{3.11}
\end{eqnarray}
with ${\bf n}_{12}=|{\bf y}_1-{\bf y}_2|/r_{12}$ and ${\bf v}_{12}={\bf v}_1-
{\bf v}_2$; scalar products are denoted with parenthesis, e.g. $(n_{12}v_1)=
{\bf n}_{12}.{\bf v}_1$. The acceleration of body 2
is obtained by exchanging the labels 1 and 2 (remembering that ${\bf n}_{12}$
and ${\bf v}_{12}$ change sign in this operation).
```

## Appendix B, Eq. (B1): the potentials of two bodies at any field point

In order: $V$ (B1a), $V_i$ (B1b), $\hat W_{ij}$ (B1c), $\hat R_i$ (B1d), $\hat X$ (B1e). The paper introduces them with: "We give first all the relevant potentials (valid all-over space-time) which are used in the obtention of the 2.5PN metric (7.2)".

```latex
\begin{mathletters}
\label{B.1}
\begin{eqnarray}
V &=& \frac{G m_1}{r_1} +\frac{G m_1}{c^2} \left(
-\frac{(n_1v_1)^2}{2r_1} +\frac{2v_1^2}{r_1} +
Gm_2 \left(-\frac{r_1}{4 r_{12}^3}-\frac{5}{4r_1r_{12}} +
\frac{r_2^2}{4r_1r_{12}^3}\right) \right) \nonumber \\
&+&\frac{2 G^2 m_1 m_2 (n_{12}v_{12})}{3c^3r_{12}^2}+
\frac{G m_1}{c^4r_1} \left( \frac{3 (n_1v_1)^4}{8}-
\frac{3(n_1v_1)^2 v_1^2}{2} + 2 v_1^4 \right) \nonumber \\
&+& \frac{G^2 m_1 m_2}{c^4} \left\{ v_1^2 \left( \frac{3 r_1^3}{16 r_{12}^5} -
\frac{37 r_1}{16 r_{12}^3} - \frac{1}{r_1 r_{12}} -
\frac{3 r_1 r_2^2}{16 r_{12}^5} + \frac{r_2^2}{r_1 r_{12}^3} \right)
\right. \nonumber \\ & & \qquad \qquad +
v_2^2 \left( \frac{3 r_1^3}{16 r_{12}^5} + \frac{3 r_1}{16 r_{12}^3} +
\frac{3}{2 r_1 r_{12}} -
  \frac{3 r_1 r_2^2}{16 r_{12}^5} + \frac{r_2^2}{2 r_1 r_{12}^3} \right)
\nonumber \\ & & \qquad \qquad +
(v_1v_2) \left(-\frac{3 r_1^3}{8 r_{12}^5} + \frac{13 r_1}{8 r_{12}^3} -
\frac{3}{r_1 r_{12}} +
  \frac{3 r_1 r_2^2}{8 r_{12}^5} - \frac{r_2^2}{r_1 r_{12}^3} \right)
\nonumber \\ & & \qquad \qquad
+(n_{12}v_1)^2 \left(-\frac{15 r_1^3}{16 r_{12}^5}+\frac{57 r_1}{16 r_{12}^3}+
  \frac{15 r_1 r_2^2}{16 r_{12}^5} \right)
\nonumber \\ & & \qquad \qquad
+(n_{12}v_2)^2 \left( -\frac{15 r_1^3}{16 r_{12}^5} -
\frac{33 r_1}{16 r_{12}^3} + \frac{7}{8 r_1 r_{12}} +
  \frac{15 r_1 r_2^2}{16 r_{12}^5} - \frac{3 r_2^2}{8 r_1 r_{12}^3} \right)
\nonumber \\ & & \qquad \qquad
+(n_{12}v_1)(n_{12}v_2) \left( \frac{15 r_1^3}{8 r_{12}^5} -
\frac{9 r_1}{8 r_{12}^3} - \frac{15 r_1 r_2^2}{8 r_{12}^5} \right)
\nonumber \\ & & \qquad \qquad+
(n_1v_1)(n_{12}v_1) \left( -\frac{3 r_1^2}{2 r_{12}^4} +
\frac{3}{4 r_{12}^2} + \frac{3 r_2^2}{4 r_{12}^4} \right)+
(n_{1}v_2)(n_{12}v_1) \left( \frac{3 r_1^2}{4 r_{12}^4}+
\frac{2}{r_{12}^2} \right) \nonumber \\ & & \qquad \qquad
+(n_1v_1)(n_{12}v_2) \left( \frac{3 r_1^2}{2 r_{12}^4} +
\frac{13}{4 r_{12}^2} - \frac{3 r_2^2}{4 r_{12}^4} \right)+
(n_{1}v_2)(n_{12}v_2) \left( -\frac{3 r_1^2}{4 r_{12}^4} -
\frac{3}{2 r_{12}^2} \right)
\nonumber \\ & & \qquad \qquad \left. +
(n_1v_1)^2 \left(-\frac{r_1}{8 r_{12}^3} +\frac{7}{8 r_1 r_{12}}
 - \frac{3 r_2^2}{8 r_1 r_{12}^3} \right)+
\frac{(n_1v_1) (n_1v_2) r_1}{2 r_{12}^3} \right\}
\nonumber \\ &+&
\frac{G^3 m_1^2 m_2}{c^4} \left( -\frac{r_1^3}{8 r_{12}^6} +
 \frac{5 r_1}{8 r_{12}^4} + \frac{3}{4 r_1 r_{12}^2} +
  \frac{r_1 r_2^2}{8 r_{12}^6} - \frac{5 r_2^2}{4 r_1 r_{12}^4} \right)
\nonumber \\ &+&
\frac{G^3 m_1 m_2^2}{c^4} \left(-\frac{r_1^3}{32 r_{12}^6} +
\frac{43 r_1}{16 r_{12}^4} + \frac{91}{32 r_1 r_{12}^2} -
  \frac{r_1 r_2^2}{16 r_{12}^6} - \frac{23 r_2^2}{16 r_1 r_{12}^4} +
  \frac{3 r_2^4}{32 r_1 r_{12}^6} \right)
\nonumber \\ &+&
\frac{G^2 m_1 m_2}{c^5} \left\{
(n_{12}v_{12})^3 \left( -\frac{7r_1^4}{4r_{12}^6} +
\frac{3 r_1^2}{2 r_{12}^4} - \frac{11}{8 r_{12}^2} +
 \frac{7 r_1^2 r_2^2}{4 r_{12}^6}\right) \right. \nonumber \\
& & \qquad \qquad
+(n_1v_{12}) (n_{12}v_{12})^2 \left( -\frac{2 r_1^3}{r_{12}^5}
-\frac{5 r_1}{r_{12}^3}
+ \frac{2 r_1 r_2^2}{r_{12}^5} \right)
-(n_1v_{12})^2 (n_{12}v_{12}) \frac{r_1^2}{r_{12}^4}
\nonumber \\ & & \qquad \qquad
+(n_{12}v_1)(n_{12}v_{12})^2 \left(\frac{5 r_1^2}{r_{12}^4}
+ \frac{5}{r_{12}^2} -\frac{5 r_2^2}{r_{12}^4} \right)
+(n_1v_{12})(n_{12}v_1)(n_{12}v_{12})\frac{6r_1}{r_{12}^3}
\nonumber \\ & & \qquad \qquad
+(n_1v_1)(n_{12}v_{12})^2 \frac{2r_1}{r_{12}^3}
-\frac{(n_{12}v_1)^2 (n_{12}v_{12})}{r_{12}^2}
-(n_1v_{12})v_1^2 \frac{8r_1}{3r_{12}^3}
+(n_{12}v_1) v_1^2 \frac{8}{3r_{12}^2}
\nonumber \\ & & \qquad \qquad
+ (n_{12}v_{12}) v_1^2 \left(-\frac{3 r_1^2}{r_{12}^4}
-\frac{8}{3 r_{12}^2} + \frac{3 r_2^2}{r_{12}^4} \right)
+(n_{12}v_{12})(v_1v_2) \frac{7}{3 r_{12}^2}
\nonumber \\ & & \qquad \qquad
+(n_1v_{12})(v_1v_2)\frac{8r_1}{3r_{12}^3}
-(n_{12}v_{1})(v_1v_2)\frac{8}{3r_{12}^2}
-(n_1v_1)v_{12}^2 \frac{2r_1}{3r_{12}^3}
\nonumber \\ & & \qquad \qquad
+(n_{12}v_{12})v_{12}^2 \left( \frac{3 r_1^4}{4 r_{12}^6} -
\frac{9 r_1^2}{10 r_{12}^4} - \frac{9}{8 r_{12}^2} -
  \frac{3 r_1^2 r_2^2}{4 r_{12}^6} \right)
\nonumber \\ & & \qquad \qquad \left.
+(n_1v_{12})v_{12}^2 \left(\frac{2 r_1^3}{5 r_{12}^5}+\frac{5 r_1}{3 r_{12}^3}
- \frac{2 r_1 r_2^2}{5 r_{12}^5} \right)
+(n_{12}v_1)v_{12}^2 \left( -\frac{r_1^2}{r_{12}^4} - \frac{5}{3 r_{12}^2}
+ \frac{r_2^2}{r_{12}^4} \right) \right\}
\nonumber \\ &+&\frac{G^3 m_1^2 m_2}{c^5} \left\{
(n_1v_{12}) \left( -\frac{4 r_1^3}{15 r_{12}^6}- \frac{8 r_1}{3r_{12}^4} +
  \frac{4 r_1 r_2^2}{15 r_{12}^6} \right)
+(n_{12}v_1) \left( -\frac{2 r_1^2}{3 r_{12}^5} + \frac{10}{3 r_{12}^3} +
\frac{2 r_2^2}{3 r_{12}^5} \right) \right.
\nonumber \\ & & \qquad \qquad \left.
+(n_1v_1) \frac{4 r_1}{3r_{12}^4}
+(n_{12}v_{12}) \left(-\frac{2r_1^2}{3 r_{12}^5} - \frac{20}{3 r_{12}^3}
+ \frac{8}{3 r_1 r_{12}^2} + \frac{34 r_2^2}{15 r_{12}^5} \right) \right\}
\nonumber \\ &+& O(6)+1\leftrightarrow 2  \ , \label{B.1a} \\
V_i &=&  \frac{G m_1 v_1^i}{r_1}+
       n^i_{12} \frac{G^2 m_1 m_2}{c^2 r_{12}^2} \left((n_1v_1)+
\frac{3 (n_{12}v_{12}) r_1}{2 r_{12}}\right)
\nonumber \\ &+& \frac{v_1^i}{c^2} \left\{
\frac{G m_1}{ r_1} \left(-\frac{(n_1v_1)^2}{2} + v_1^2 \right)
+ G^2 m_1 m_2 \left( -\frac{3 r_1}{4 r_{12}^3} +
\frac{r_2^2}{4 r_1 r_{12}^3} - \frac{5}{4 r_1 r_{12}} \right) \right\}
+v_2^i \frac{G^2 m_1 m_2 r_1}{2 c^2 r_{12}^3} \nonumber \\
&+& \frac{n_{12}^i}{c^3} G^2 m_1 m_2 \left[
(n_{12}v_{12})^2 \frac{5 r_1^2}{2 r_{12}^4}
+(n_{12}v_{12})(n_1v_{12}) \frac{3 r_1}{2 r_{12}^3}
-v_{12}^2 \frac{r_1^2}{2 r_{12}^4} +\frac{v_1^2}{2 r_{12}^2} \right.
\nonumber \\ & & \left. \qquad
+ G m_1 \left( -\frac{r_1^2}{6 r_{12}^5} + \frac{r_2^2}{6 r_{12}^5} +
\frac{1}{2 r_{12}^3} \right) \right]
-\frac{G^2 m_1 m_2 v_1^i}{3c^3 r_{12}^2} \nonumber \\
&+&\frac{G^2 m_1 m_2 v_{12}^i}{c^3 r_{12}^2} \left(
- \frac{3 r_1^2}{2 r_{12}^4}(n_{12}v_{12}) -
\frac{2 r_1}{3 r_{12}^3} (n_1v_{12}) +\frac{(n_{12}v_1)}{2 r_{12}^2} \right)
+ O(4)+1\leftrightarrow 2 \ , \label{B.1b} \\
{\hat W}_{ij} &=&
\delta^{ij} \left(  -\frac{G m_1 v_1^2}{r_1}-\frac{G^2 m_1^2}{4 r_1^2}+
\frac{G^2 m_1 m_2}{r_{12} S} \right)
+\frac{G m_1 v_1^i v_1^j}{r_1} +\frac{G^2 m_1^2 n_1^i n_1^j}{4 r_1^2}
\nonumber \\ &+&
G^2 m_1 m_2 \left\{ \frac{1}{S^2} \left(n_1^{(i} n_2^{j)}
+ 2n_1^{(i} n_{12}^{j)} \right)-
n_{12}^i n_{12}^j \left(\frac{1}{S^2}+\frac{1}{r_{12} S} \right) \right\}
\nonumber \\ &+&
\frac{G^2 m_1 m_2}{c r_{12}^2} \left(-\frac{(n_{12}v_{12})}{2} \delta^{ij}-
\frac{3 (n_{12}v_{12})}{2} n_{12}^i n_{12}^j +
2n_{12}^{(i} v_{12}^{j)} \right)
\nonumber \\ &+&O(2)+1\leftrightarrow 2 \ , \label{B.1c} \\
{\hat R}_i &=&
G^2 m_1 m_2 n_{12}^i \left\{ -\frac{(n_{12}v_1)}{2S} \left(\frac{1}{S}+
\frac{1}{r_{12}} \right)
-\frac{2 (n_2v_1)}{S^2}+\frac{3(n_2v_2)}{2S^2} \right\}
\nonumber \\ &+&
n_1^i \left\{ \frac{G^2 m_1^2 (n_1v_1)}{8 r_1^2}
+\frac{G^2 m_1 m_2}{S^2} \left( 2 (n_{12}v_1) - \frac{3(n_{12}v_2)}{2}
+2 (n_2v_1) - \frac{3(n_2v_2)}{2} \right) \right\}
\nonumber \\ &+&
v_1^i \left\{ -\frac{G^2 m_1^2}{8r_1^2}+
G^2 m_1 m_2 \left( \frac{1}{r_1 r_{12}} + \frac{1}{2 r_{12}S} \right) \right\}
-v_2^i\frac{G^2 m_1 m_2}{r_1 r_{12}}
+ v_1^i \frac{G^2 m_1 m_2 (n_{12}v_1)}{2c r_{12}^2}
\nonumber \\ &+& n_{12}^i \frac{G^2 m_1 m_2}{c r_{12}^2} \left(
-\frac{3(n_{12}v_1)^2}{4} + \frac{v_1^2}{4} \right)
+ O(2)+1\leftrightarrow 2 \ , \label{B.1d} \\
{\hat X} &=&
\frac{G^2 m_1^2}{8r_1^2} \left((n_1v_1)^2-v_1^2\right)
+G^2 m_1 m_2 v_1^2 \left(\frac{1}{r_1 r_{12}}
+ \frac{1}{r_1 S} + \frac{1}{r_{12} S} \right)
\nonumber \\ &+&
G^2 m_1 m_2 \left\{v_2^2 \left(-\frac{1}{r_1 r_{12}} + \frac{1}{r_1 S}
+ \frac{1}{r_{12} S} \right)
-\frac{(v_1v_2)}{S} \left(\frac{2}{r_1}+
\frac{3}{2 r_{12}} \right)
-\frac{(n_{12}v_1)^2}{S} \left(\frac{1}{S}+\frac{1}{r_{12}}\right) \right.
\nonumber \\ & & \qquad \qquad
-\frac{(n_{12}v_2)^2}{S} \left(\frac{1}{S} + \frac{1}{r_{12}} \right)
+\frac{3(n_{12}v_1)(n_{12}v_2)}{2S} \left(\frac{1}{S} +
\frac{1}{r_{12}}\right)+\frac{2(n_{12}v_1)(n_1v_1)}{S^2}
\nonumber \\ & & \qquad \qquad
-\frac{5(n_{12}v_2)(n_1v_1)}{S^2}
-\frac{(n_1v_1)^2}{S} \left(\frac{1}{S} + \frac{1}{r_1} \right)
+\frac{2(n_{12}v_2)(n_1v_2)}{S^2} \nonumber \\ & & \qquad \qquad
+\frac{2(n_1v_1)(n_1v_2)}{S} \left(\frac{1}{S} + \frac{1}{r_1} \right)
-\frac{(n_1v_2)^2}{S} \left(\frac{1}{S} + \frac{1}{r_1}\right)
-\frac{2(n_{12}v_2)(n_2v_1)}{S^2} \nonumber \\
& & \left. \qquad \qquad +\frac{2(n_1v_2)(n_2v_1)}{S^2}
-\frac{3(n_1v_1)(n_2v_2)}{2 S^2} \right\}+\frac{G^3 m_1^3}{12 r_1^3}
\nonumber \\
&+& G^3 m_1^2 m_2 \left(\frac{1}{2 r_1^3} + \frac{1}{16 r_2^3} +
\frac{1}{16 r_1^2 r_2} - \frac{r_2^2}{2 r_1^2 r_{12}^3} +
\frac{r_2^3}{2 r_1^3 r_{12}^3} - \frac{r_1^2}{32 r_2^3 r_{12}^2}
- \frac{3}{16 r_2 r_{12}^2} +
\frac{15 r_2}{32 r_1^2 r_{12}^2} \right. \nonumber \\
& & \left. \qquad \qquad - \frac{r_2^2}{2 r_1^3 r_{12}^2}
- \frac{r_2}{2 r_1^3 r_{12}} - \frac{r_{12}^2}{32 r_1^2 r_2^3} \right)
+G^3 m_1 m_2^2 \left( -\frac{1}{2 r_{12}^3} + \frac{r_2}{2 r_1 r_{12}^3}
-\frac{1}{2 r_1 r_{12}^2} \right) \nonumber \\
&+& \frac{G^2 m_1 m_2}{cr_{12}^2} \left( -\frac{3(n_{12}v_1)^2
(n_{12}v_{12})}{4}+
  \frac{3 (n_{12}v_1) (n_{12}v_{12})^2}{4} -
  \frac{3 (n_{12}v_{12})^3}{2} - \frac{(n_{12}v_1) v_{12}^2}{4}
\right. \nonumber \\ & &\qquad \qquad  \left.
+ \frac{3 (n_{12}v_{12}) v_{12}^2}{2} +
  \frac{(n_{12}v_1) v_1^2}{2} - \frac{(n_{12}v_{12}) v_1^2}{4}
- \frac{(n_{12}v_1) (v_1v_2)}{2} +
  \frac{(n_{12}v_{12}) (v_1v_2)}{2} \right) \nonumber \\
&+& \frac{G^3 m_1^2 m_2}{c} \left\{
(n_{12}v_{12}) \left( \frac{3}{8 r_1^3} + \frac{3 r_1}{8 r_{12}^4}
-\frac{3 r_2^2}{4 r_1 r_{12}^4} +
\frac{3 r_2^4}{8 r_1^3 r_{12}^4} - \frac{3}{2 r_{12}^3}
+ \frac{1}{4 r_1 r_{12}^2} -
  \frac{3 r_2^2}{4 r_1^3 r_{12}^2}\right) \right.
\nonumber \\ & & \left. \qquad \qquad
+(n_1v_{12}) \left( \frac{1}{r_{12}^3} - \frac{r_2^2}{r_1^2 r_{12}^3}
+ \frac{1}{r_1^2 r_{12}} \right) \right\}
+O(2)+1\leftrightarrow 2 \ . \label{B.1e}
\end{eqnarray}
\end{mathletters}
```

## Appendix B, Eq. (B2): the potentials regularized at body 1

"The values of the potentials at the location of body 1, following from (3.6), are":

```latex
\begin{mathletters}
\label{B.2}
\begin{eqnarray}
(V)_1 &=& {Gm_2\over r_{12}} \left\{1+{1\over c^2}\left[-{3\over 2}
{Gm_1\over r_{12}}+2v_2^2-{1\over 2}(n_{12}v_2)^2\right]
+{4\over 3}{Gm_1\over r_{12}c^3}(n_{12}v_{12}) \right.\nonumber\\
& &\qquad  +{Gm_1\over r_{12}c^4}\left[ {11\over 2}{Gm_1\over r_{12}}+
{5\over 4}{Gm_2\over r_{12}}+{15\over
8}v_1^2-{7\over 4}(v_1v_2)-{25\over 8}v_2^2 \right. \nonumber\\
& &\qquad \qquad \quad \left.+ {1\over 8}(n_{12}v_1)^2-
{25\over 4}(n_{12}v_1)(n_{12}v_2)+{33\over 8}(n_{12}v_2)^2\right] \nonumber\\
& &\left. \qquad+{1\over c^4}\left[2v_2^4-{3\over 2}(n_{12}v_2)^2v_2^2+
{3\over 8}(n_{12}v_2)^4\right]\right\}\nonumber\\
&+&{G^2m_1m_2\over c^5r_{12}^2}\left\{ -{2\over 5}{Gm_1\over r_{12}}
(n_{12}v_1) -{26\over 15}{Gm_2\over r_{12}}(n_{12}v_1)
+{22\over 5}{Gm_1\over r_{12}}(n_{12}v_2)\right.\nonumber\\
& &\qquad \qquad+{2\over 5}{Gm_2\over r_{12}}(n_{12}v_2)
-{32\over 15}(n_{12}v_1)v_1^2
+{48\over 5}(n_{12}v_1)(v_1v_2)-{122\over 15}(n_{12}v_1)v_2^2\nonumber\\
& &\qquad \qquad+{92\over 15}(n_{12}v_2)v_1^2-{184\over 15}
(n_{12}v_2)(v_1v_2)+{34\over 5}(n_{12}v_2)v_2^2
+2(n_{12}v_1)^3\nonumber\\
& &\qquad \qquad \left.- 10(n_{12}v_1)^2(n_{12}v_2)
+12(n_{12}v_1)(n_{12}v_2)^2
-4(n_{12}v_2)^3\right\}+O(6)   \ , \label{B.2a}
\\
(V_i)_1 &=& {Gm_2\over r_{12}}\biggl\{v_2^i
+{v_2^i\over c^2}\left[-2{Gm_1\over r_{12}}+v_2^2-
{1\over 2}(n_{12}v_2)^2\right] +{1\over 2}{Gm_1\over
r_{12}c^2}v_1^i \nonumber\\
& &\qquad+{Gm_1\over r_{12}c^2}n_{12}^i\left[-{3\over 2}(n_{12}v_1)+
{1\over 2}(n_{12}v_2)\right]\nonumber\\
& &\qquad+{Gm_1\over c^3r_{12}}v_1^i(n_{12}v_1)
-{Gm_1\over c^3r_{12}}v_2^i\left[{5\over 3}(n_{12}v_1)-
{2\over 3}(n_{12}v_2)\right]\nonumber\\
& &\qquad+{Gm_1\over c^3r_{12}}n_{12}^i\left[{2\over 3}{Gm_1\over r_{12}}
-{1\over 3}{Gm_2\over r_{12}}+v_1^2-(v_1v_2)-(n_{12}v_{12})^2 \right]
\biggr\}+O(4) \ , \label{B.2b}
\\
({\hat W}_{ij})_1 &=& {Gm_2\over r_{12}}\biggl\{v_2^{ij}-\delta^{ij}v_2^2+
{Gm_1\over r_{12}}\left[-2n_{12}^{ij}+\delta^{ij}\right]\nonumber\\
&&\qquad+{Gm_2\over 4r_{12}}\left[n_{12}^{ij}-\delta^{ij}\right]+4{Gm_1\over
cr_{12}}n_{12}^{(i}v_{12}^{j)}\nonumber\\
&&\qquad-{Gm_1\over cr_{12}}(n_{12}v_{12})\left(3n_{12}^{ij}+\delta^{ij}\right)
\biggr\}+O(2) \ , \label{B.2c} \\
({\hat R}_i)_1 &=& {G^2m_1m_2\over r_{12}^2}\left[-{3\over 4}v_1^i+
{5\over 4}v_2^i-{1\over
2}(n_{12}v_1)n_{12}^i-{1\over 2}(n_{12}v_2)n_{12}^i\right] \nonumber\\
&+&{G^2m_2^2\over r_{12}^2}\left[-{1\over 8}v_2^i+
{1\over 8}(n_{12}v_2)n_{12}^i\right] \nonumber\\
&+&{G^2m_1m_2\over cr_{12}^2}\biggl\{n_{12}^i\left[ {1\over 4}v_1^2-
{1\over 4}v_2^2-{3\over 4}(n_{12}v_1)^2+{3\over
4}(n_{12}v_2)^2\right]\nonumber\\
&&\qquad +{1\over 2}(n_{12}v_1)v_1^i-{1\over 2}(n_{12}v_2)v_2^i\biggr\}+O(2)
\ , \label{B.2d} \\
({\hat X})_1 &=& {G^2m_1m_2\over r_{12}^2}\left[-{3\over 2}{Gm_1\over r_{12}}+
{1\over 4}v_1^2-
2(v_1v_2)+{9\over 4}v_2^2 \right.\nonumber\\
&&\quad-\left.{11\over 4}(n_{12}v_1)^2+{9\over 2}(n_{12}v_1)(n_{12}v_2)
-{11\over 4}(n_{12}v_2)^2\right] \nonumber\\
&+&{G^2m_2^2\over r_{12}^2}\left[{1\over 12}{Gm_2\over r_{12}}-
{1\over 8}v_2^2+{1\over 8}(n_{12}v_2)^2\right]\nonumber\\
&+&\frac{G^2 m_1 m_2}{cr_{12}^2} \biggl\{ -3 (n_{12}v_1)^3 +
\frac{15}{2}(n_{12}v_1)^2 (n_{12}v_2) -
\frac{15}{2} (n_{12}v_1) (n_{12}v_2)^2 + 3 (n_{12}v_2)^3 \nonumber\\
&&\quad+ 3 (n_{12}v_1) v_1^2 - \frac{5}{2} (n_{12}v_2) v_1^2
- 5 (n_{12}v_{12}) (v_1v_2) + \frac{5}{2} (n_{12}v_1) v_2^2 -
3 (n_{12}v_2) v_2^2 \nonumber \\
&& \quad-\frac{3}{2} \frac{G m_1}{r_{12}} (n_{12}v_{12})
-\frac{5}{2} \frac{G m_2}{ r_{12}}(n_{12}v_{12}) \bigg\}
+O(2) \ . \label{B.2e}
\end{eqnarray}
\end{mathletters}
```

## Appendix B, Eq. (B3): gradients of the potentials at body 1

"The gradients of the potentials computed at body 1 (needed for the equations of motion) are":

```latex
\begin{mathletters}
\label{B.3}
\begin{eqnarray}
(\partial_iV)_1 &=& - {Gm_2\over r_{12}^2} n_{12}^i
 +{Gm_2\over r_{12}^2c^2} n_{12}^i \left[ - 2v^2_2
  + {3\over 2} (n_{12}v_2)^2 + {Gm_1\over r_{12}}\right]-
{Gm_2\over r_{12}^2c^2}(n_{12}v_2)v_2^i \nonumber\\
&+&{Gm_2\over r_{12}^2c^4} n_{12}^i \left[ -2v^4_2 +
{9\over 2} v^2_2 (n_{12}v_2)^2 -{15\over 8}(n_{12}v_2)^4 \right] \nonumber\\
&+&{G^2m_1m_2\over r_{12}^3c^4}n_{12}^i \left[ -{3\over 4} v^2_1 -
{3\over 4} v^2_2+{7\over 2} (v_1v_2) - {13\over 2} (n_{12}v_1)^2 \right.
\nonumber\\
&&\qquad \qquad \quad \left. +7(n_{12}v_1)(n_{12}v_2)+
{1\over 2}(n_{12}v_2)^2 -{1\over 4}{Gm_1\over r_{12}}-
{1\over 2}{Gm_2\over r_{12}} \right]\nonumber\\
 &+& {Gm_2\over r_{12}^2c^4}v_2^i\left[ -3v^2_2(n_{12}v_2)+
{3\over 2}(n_{12}v_2)^3 -{17\over 4}{Gm_1\over
r_{12}}(n_{12}v_1)+{9\over 4}{Gm_1\over r_{12}}(n_{12}v_2)\right]\nonumber\\
&+& {G^2m_1m_2\over r_{12}^3c^4} v_1^i\left[ -{9\over 4}(n_{12}v_2) +
{9\over 4}(n_{12}v_1)\right] \nonumber \\
&+& \frac{G^2 m_1 m_2}{c^5 r_{12}^3} n_{12}^i \left(-10 (n_{12}v_1)^3 +
  10 (n_{12}v_1)^2 (n_{12}v_2) + 10 (n_{12}v_1) (n_{12}v_2)^2 \right.
\nonumber \\ & & \qquad \qquad \quad  -
  10 (n_{12}v_2)^3 +\frac{42}{5} (n_{12}v_1) v_1^2
- \frac{22}{5} (n_{12}v_2) v_1^2 - \frac{24}{5} (n_{12}v_1) (v_1v_2)
\nonumber \\ & & \qquad \qquad \quad  -
 \frac{16}{5} (n_{12}v_2) (v_1v_2) -
  \frac{18}{5} (n_{12}v_1) v_2^2 + \frac{38}{5} (n_{12}v_2) v_2^2
+ \frac{88}{15} \frac{G m_1}{r_{12}} (n_{12}v_1) \nonumber \\
& & \qquad \qquad \quad \left.
- \frac{52}{15}\frac{G m_2}{r_{12}}(n_{12}v_1)  -
  \frac{68}{15} \frac{G m_1}{r_{12}} (n_{12}v_2)+
  \frac{24}{5} \frac{G m_2}{r_{12}}  (n_{12}v_2)\right) \nonumber \\
&+& \frac{G^2 m_1 m_2}{c^5 r_{12}^3} v_1^i
\left(10 (n_{12}v_1)^2 - 8 (n_{12}v_1) (n_{12}v_2) - 2 (n_{12}v_2)^2 -
  \frac{16}{15} \frac{G m_1}{r_{12}} +
  \frac{44}{15} \frac{G m_2}{r_{12}} \right.
\nonumber \\ & & \qquad \qquad \quad \left. - \frac{62}{15} v_1^2 +
\frac{44}{15} (v_1v_2) + \frac{6}{5} v_2^2 \right)
\nonumber \\
&+&\frac{G^2 m_1 m_2}{c^5 r_{12}^3} v_2^i \left(-6 (n_{12}v_1)^2 +
6 (n_{12}v_2)^2 + \frac{12}{5} \frac{G m_1}{r_{12}} -
\frac{8}{5} \frac{G m_2}{r_{12}} + \frac{14}{5} v_1^2 \right.
\nonumber \\ & & \qquad \qquad \quad \left. - \frac{4}{15} (v_1v_2)-
 \frac{38}{15} v_2^2\right) +O(6) \ , \label{B.3a}\\
(\partial_jV_i)_1 &=& -{Gm_2\over r_{12}^2}n_{12}^jv_2^i+
{Gm_2\over r_{12}^2c^2}\left\{n_{12}^{ij}{Gm_1\over r_{12}}\left[ -{3\over 2}
(n_{12}v_1) + {5\over 2}(n_{12}v_2)\right]\right.\nonumber\\
& & \qquad
+{Gm_1\over 2r_{12}}n_{12}^jv_1^i+n_{12}^jv_2^i\left[-v_2^2+{3\over 2}
(n_{12}v_2)^2+{Gm_1\over 2r_{12}}\right]\nonumber\\
& & \qquad \left.-{Gm_1\over r_{12}}n_{12}^iv_2^j-v_2^{ij}(n_{12}v_2)\right\}
\nonumber \\
&+& \frac{G^2 m_1 m_2}{c^3 r_{12}^3} \biggl\{ n_{12}^{ij} \left(
-5 (n_{12}v_{12})^2 + v_{12}^2
+ \frac{1}{3}\frac{G m_1}{r_{12}} +
\frac{1}{3} \frac{G m_2}{r_{12}} \right)\nonumber\\
&&\qquad +6 n_{12}^{(i} v_1^{j)} (n_{12}v_{12})
-6 n_{12}^{(i} v_2^{j)} (n_{12}v_{12})
-\frac{4}{3} v_{12}^{ij}\biggr\}+O(4) \ , \label{B.3b} \\
(\partial_k{\hat W}_{ij})_1 &=& {Gm_2\over r_{12}^2}\biggl\{-n_{12}^kv_2^{ij}+
\delta^{ij}n_{12}^k\left[-{Gm_1\over 2r_{12}}+{Gm_2\over 2r_{12}}+v_2^2\right]
\nonumber\\
&&\quad\quad+\delta^{k(i}n_{12}^{j)}\left[-{3Gm_1\over 2r_{12}}+{Gm_2\over 2r_{
12}}\right]
+n_{12}^{ijk}\left[{2Gm_1\over r_{12}}-{Gm_2\over r_{12}}\right]\biggr\} +O(2)
\ , \label{B.3c}\\
(\partial_j{\hat R}_i)_1 &=&
{Gm_2\over r_{12}^2}\biggl\{\delta^{ij}\left[-{5\over 8}{Gm_1\over r_{12}}(n_{1
2}v_1)+{1\over 4}{Gm_1\over r_{12}}(n_{12}v_2)+
{1\over 8}{Gm_2\over r_{12}}(n_{12}v_2)\right]\nonumber\\
&&\quad\quad +n_{12}^{ij}\left[{1\over 2}{Gm_1\over r_{12}}(n_{12}v_1)+
{1\over 2}{Gm_1\over r_{12}}(n_{12}v_2)-
{1\over 2}{Gm_2\over r_{12}}(n_{12}v_2)\right]\nonumber\\
&&\quad\quad +n_{12}^i\left[{1\over 4}{Gm_1\over r_{12}}v_1^j-
{5\over 8}{Gm_1\over r_{12}}v_2^j+
{1\over 8}{Gm_2\over r_{12}}v_2^j\right]\nonumber\\
&&\quad\quad +n_{12}^j\left[{7\over 8}{Gm_1\over r_{12}}v_1^i-
{9\over 8}{Gm_1\over r_{12}}v_2^i+
{1\over 4}{Gm_2\over r_{12}}v_2^i\right]\biggr\}+O(2) \ , \label{B.3d}
\\
(\partial_i{\hat X})_1 &=&
{G^2m_1m_2\over r_{12}^3}\biggl\{n_{12}^{i}\left[{1\over 2}{Gm_1\over r_{12}}-
{3\over 2}{Gm_2\over r_{12}}+
{9\over 2}(n_{12}v_1)^2-8(n_{12}v_1)(n_{12}v_2)\right.\nonumber\\
& &\qquad \qquad \quad +\left.{9\over 2}(n_{12}v_2)^2
+{7\over 4}(v_1v_2)-2v_2^2\right]\nonumber\\
&&\quad\quad -{5\over 2}(n_{12}v_1)v_1^i
+3(n_{12}v_2)v_1^i+{5\over 4}(n_{12}v_1)v_2^i-
{5\over 2}(n_{12}v_2)v_2^i\biggr\}\nonumber\\
&+&{G^2m_2^2\over r_{12}^3}\biggl\{n_{12}^{i}\left[-{1\over 4}
{Gm_2\over r_{12}}- {1\over 2}(n_{12}v_2)^2+{1\over 4}v_2^2\right]
+{1\over 4}(n_{12}v_2)v_2^i\biggr\}\nonumber\\
&+&\frac{G^3 m_1 m_2^2}{cr_{12}^4}\biggl\{ 3n_{12}^i (n_{12}v_{12})
- 2v_{12}^i\biggr\} +O(2) \ . \label{B.3e}
\end{eqnarray}
\end{mathletters}
```

## Section VII, Eq. (7.2): the explicit metric of two bodies at any field point

This is Eq. (2.5) with the potentials (B1) inserted and like terms combined. $g_{00}+1$ is (7.2a), $g_{0i}$ is (7.2b), $g_{ij}-\delta_{ij}$ is (7.2c). The sentence after the block is the paper's.

```latex
\begin{mathletters}
\label{7.2}
\begin{eqnarray}
g_{00}+1 & = &
\frac{2 G m_1}{c^2 r_1}+
\frac{1}{c^4}
	\left[\frac{G m_1}{r_1} \left(-(n_1v_1)^2+4 v_1^2 \right)
       -2 \frac{G^2 m_1^2}{r_1^2}\right. \nonumber \\ & &
	\qquad \qquad \quad+
	\left. G^2 m_1 m_2 \left(-\frac{2}{r_1 r_2}-
				\frac{r_1}{2 r_{12}^3}+
				\frac{r_1^2}{2 r_2 r_{12}^3}-
				\frac{5}{2 r_2 r_{12}}\right)\right]
+\frac{4 G^2 m_1 m_2}{3 c^5 r_{12}^2} (n_{12}v_{12})
\nonumber \\&+&
\frac{1}{c^6} \left[
\frac{G m_1}{r_1} \left(\frac{3}{4} (n_1v_1)^4-3 (n_1v_1)^2 v_1^2+
							4 v_1^4\right)+
\frac{G^2 m_1^2}{r_1^2} \left(3 (n_1v_1)^2-v_1^2\right)
				+2 \frac{G^3 m_1^3}{r_1^3}\right.
				\nonumber \\ &+&
G^2 m_1 m_2 \left(v_1^2 \left(\frac{3 r_1^3}{8 r_{12}^5}-
	         \frac{3 r_1^2 r_2}{8 r_{12}^5}-
		 \frac{3 r_1 r_2^2}{8 r_{12}^5}+
   		 \frac{3 r_2^3}{8 r_{12}^5}-\frac{37 r_1}{8 r_{12}^3}+
		 \frac{r_1^2}{r_2 r_{12}^3}+\frac{3 r_2}{8 r_{12}^3}
\right. \right. \nonumber
		 \\ & & \left.
		 \qquad \qquad \quad +
		 \frac{2 r_2^2}{r_1 r_{12}^3}+\frac{6}{r_1 r_{12}}-
   		 \frac{5}{r_2 r_{12}}-\frac{8r_{12}}{r_1r_2S}
		+\frac{16}{r_{12} S}\right)\nonumber \\  & &
		 		\qquad \quad+
	         (v_1v_2) \left(\frac{8}{r_1 r_2}-\frac{3 r_1^3}{4 r_{12}^5}+
		 \frac{3 r_1^2 r_2}{4 r_{12}^5}+
  		 \frac{13 r_1}{4 r_{12}^3}-\frac{2 r_1^2}{r_2 r_{12}^3}-
		 \frac{6}{r_1 r_{12}}-
		 \frac{16}{r_1 S}- \frac{12}{r_{12} S}\right) \nonumber \\ & &
		    		\qquad \quad+
	   (n_{12}v_1)^2 \left(-\frac{15 r_1^3}{8 r_{12}^5}+
			       \frac{15 r_1^2 r_2}{8 r_{12}^5}+
  		    	       \frac{15 r_1 r_2^2}{8 r_{12}^5}-
			       \frac{15 r_2^3}{8 r_{12}^5}+
  		    	       \frac{57 r_1}{8 r_{12}^3}-
			       \frac{3 r_1^2}{4 r_2 r_{12}^3}-
  		    	       \frac{33 r_2}{8 r_{12}^3}
			       \right. \nonumber \\ & & \left.
			       \qquad \qquad \qquad \quad+
			       \frac{7}{4 r_2 r_{12}} -
			       \frac{16}{S^2}-\frac{16}{r_{12} S}\right)
			\nonumber \\ & & \quad \qquad +
	   (n_{12}v_1) (n_{12}v_2) \left(\frac{15 r_1^3}{4 r_{12}^5}-
			\frac{15 r_1^2 r_2}{4 r_{12}^5}-
  			\frac{9 r_1}{4 r_{12}^3}+\frac{12}{S^2}+
			\frac{12}{r_{12} S}\right) \nonumber
			\\ & & \quad \qquad+
	   (n_1v_1)^2 \left(\frac{2}{r_1 r_2}-\frac{r_1}{4 r_{12}^3}-
			    \frac{3 r_2^2}{4 r_1 r_{12}^3}+
  		   	    \frac{7}{4 r_1 r_{12}}-\frac{8}{S^2}-
			    \frac{8}{r_1 S}\right) \nonumber
			\\ & & \quad \qquad+
	   (n_1v_1) (n_1v_2) \left(\frac{r_1}{r_{12}^3}+\frac{16}{S^2}+
				   \frac{16}{r_1 S}\right)-
	   (n_1v_2)^2 \left(\frac{8}{S^2}+\frac{8}{r_1 S}\right) \nonumber
	   \\ & & \qquad \quad+
	   (n_{12}v_1) (n_1v_1) \left(-\frac{3 r_1^2}{r_{12}^4}+
				\frac{3 r_2^2}{2 r_{12}^4}+
				\frac{3}{2 r_{12}^2}+\frac{16}{S^2}\right)+
	   \frac{16 (n_1v_2) (n_2v_1)}{S^2}
			\nonumber \\ & & \qquad \quad +
	   (n_{12}v_2) (n_1v_1) \left(\frac{3 r_1^2}{r_{12}^4}-
				      \frac{3 r_2^2}{2 r_{12}^4}+
				      \frac{13}{2 r_{12}^2}-
				      \frac{40}{S^2}\right)-
		\frac{12(n_1v_1) (n_2v_2)}{S^2} \nonumber
			\\ & & \left. \qquad \quad +
	   (n_{12}v_1) (n_1v_2) \left(\frac{3 r_1^2}{2 r_{12}^4}+
				      \frac{4}{r_{12}^2}+
				      \frac{16}{S^2} \right)+
	   (n_{12}v_2) (n_1v_2) \left(\frac{-3 r_1^2}{2 r_{12}^4}-
				      \frac{3}{r_{12}^2}+\frac{16}{S^2}\right)
					\right)
			\nonumber \\ &+&
G^3 m_1^2 m_2 \left(\frac{4}{r_1^3}+\frac{1}{2 r_2^3}+\frac{9}{2 r_1^2 r_2}-
	       \frac{r_1^3}{4 r_{12}^6}+\frac{3 r_1^4}{16 r_2 r_{12}^6}-
	       \frac{r_1^2 r_2}{8 r_{12}^6}+
  	       \frac{r_1 r_2^2}{4 r_{12}^6}-
	       \frac{r_2^3}{16 r_{12}^6}+
	       \frac{5 r_1}{4 r_{12}^4}
			\right. \nonumber \\ & & \qquad \
		-\frac{23 r_1^2}{8 r_2 r_{12}^4}+
	       \frac{43 r_2}{8 r_{12}^4}-
  	       \frac{5 r_2^2}{2 r_1 r_{12}^4}-\frac{3}{r_{12}^3}+
	       \frac{3 r_1}{r_2 r_{12}^3}+\frac{r_2}{r_1 r_{12}^3}-
	       \frac{5 r_2^2}{r_1^2 r_{12}^3}+
  	       \frac{4 r_2^3}{r_1^3 r_{12}^3}+\frac{3}{2 r_1 r_{12}^2}
	\nonumber \\ & & \left. \left. \qquad \ -
  	       \frac{r_1^2}{4 r_2^3 r_{12}^2} +
	       \frac{3}{16 r_2 r_{12}^2}+
  	       \frac{15 r_2}{4 r_1^2 r_{12}^2}-\frac{4 r_2^2}{r_1^3 r_{12}^2}+
  	       \frac{5}{r_1^2 r_{12}}+\frac{5}{r_1 r_2 r_{12}}-
	       \frac{4 r_2}{r_1^3 r_{12}}-
	       \frac{r_{12}^2}{4 r_1^2 r_2^3}\right)\right] \nonumber \\ &+&
\frac{1}{c^7} \left[
G^2 m_1 m_2
   \left((n_{12}v_{12})^2 (n_1v_1) \left(-\frac{8 r_1^3}{r_{12}^5}-
	 \frac{16 r_1}{r_{12}^3}\right)+
     	(n_{12}v_{12})^2 (n_1v_2) \left(\frac{8 r_1^3}{r_{12}^5}+
	\frac{5 r_1}{r_{12}^3}\right) \right. \right.
\nonumber \\ & & \qquad \qquad \quad +
     	(n_{12}v_{12})^3 \left(-\frac{7 r_1^4}{2 r_{12}^6}+
	\frac{7 r_1^2 r_2^2}{2 r_{12}^6}-
        \frac{11 r_1^2}{r_{12}^4}-\frac{37}{4 r_{12}^2}\right)
\nonumber \\ & & \qquad \qquad \quad +
     	(n_{12}v_1) (n_{12}v_{12})^2 \left(\frac{20 r_1^2}{r_{12}^4}-
	\frac{11}{2 r_{12}^2}\right)-
     	4 (n_{12}v_{12}) (n_1v_1)^2 \frac{r_1^2}{r_{12}^4}
\nonumber \\ & & \qquad \qquad \quad
     	+4 (n_{12}v_{12}) (n_1v_1) (n_1v_2)\frac{r_1^2}{r_{12}^4}+
     	(n_{12}v_1)^2 (n_1v_1)\frac{r_1}{r_{12}^3}
\nonumber \\ & & \qquad \qquad \quad+
     	22 (n_{12}v_1) (n_{12}v_{12}) (n_1v_1)\frac{r_1}{r_{12}^3}-
     	(n_{12}v_1)^2 (n_1v_2)\frac{r_1}{r_{12}^3}
\nonumber \\ & & \qquad \qquad \quad+
	4 (n_{12}v_1) (n_{12}v_{12}) (n_1v_2)\frac{r_1}{r_{12}^3}+
     	11 (n_{12}v_1)^2 (n_{12}v_{12})\frac{1}{2 r_{12}^2}
\nonumber \\ & & \qquad \qquad \quad+
     	(n_1v_2) v_{12}^2 \left(-\frac{8 r_1^3}{5 r_{12}^5}-
		\frac{2 r_1}{3 r_{12}^3}\right)+
     	(n_1v_1) v_{12}^2 \left(\frac{8 r_1^3}{5 r_{12}^5}+
		\frac{11 r_1}{3 r_{12}^3}\right)
\nonumber \\ & & \qquad \qquad \quad -
     	(n_{12}v_1) v_{12}^2 \left(\frac{4 r_1^2}{r_{12}^4}+
	\frac{5}{2 r_{12}^2}\right)-
     	(n_{12}v_{12}) v_1^2 \left(\frac{12 r_1^2}{r_{12}^4}+
		\frac{5}{2 r_{12}^2}\right)
\nonumber \\ & & \qquad \qquad \quad+
     	(n_{12}v_{12}) v_{12}^2 \left(\frac{3 r_1^4}{2 r_{12}^6}-
	\frac{3 r_1^2 r_2^2}{2 r_{12}^6}+
       \frac{7 r_1^2}{r_{12}^4}+\frac{27}{4 r_{12}^2}\right)
\nonumber \\ & & \qquad \qquad \quad-
     	29 (n_1v_1) v_1^2\frac{r_1}{3 r_{12}^3}+
	(n_1v_2) v_1^2\frac{r_1}{r_{12}^3}+
     	5 (n_{12}v_1) v_1^2\frac{1}{r_{12}^2}
\nonumber \\ & & \qquad \qquad \quad+
     	(n_{12}v_{12}) (v_1v_2) \left(\frac{12 r_1^2}{r_{12}^4}+
		\frac{3}{r_{12}^2}\right)+
     	8 (n_1v_1) (v_1v_2)\frac{r_1}{r_{12}^3}
\nonumber \\ & & \left. \qquad \quad \qquad+
	2 (n_1v_2) (v_1v_2)\frac{r_1}{3 r_{12}^3}-
     	5 (n_{12}v_1) (v_1v_2)\frac{1}{r_{12}^2}\right)
\nonumber \\ &+ &
G^3 m_1^2 m_2 \left((n_1v_{12}) \left(-\frac{8 r_1^3}{15 r_{12}^6}+
		\frac{8 r_1 r_2^2}{15 r_{12}^6}-
        \frac{16 r_1}{3 r_{12}^4}+\frac{8}{r_{12}^3}-
	\frac{8 r_2^2}{r_1^2 r_{12}^3}+
  	\frac{8}{r_1^2 r_{12}}\right) \right.
\nonumber \\ & & \qquad \quad +(n_{12}v_1)
      	\left(-\frac{4 r_1^2}{3 r_{12}^5}+
		\frac{4 r_2^2}{3 r_{12}^5}+\frac{20}{3 r_{12}^3}\right)+
     	8 (n_1v_1)\frac{r_1}{3 r_{12}^4}
\nonumber \\ & & \qquad \quad +
     (n_{12}v_{12}) \left(\frac{3}{r_1^3}-\frac{4 r_1^2}{3 r_{12}^5}+
		\frac{68 r_2^2}{15 r_{12}^5}+
       \frac{3 r_1}{r_{12}^4}-\frac{6 r_2^2}{r_1 r_{12}^4}+
	\frac{3 r_2^4}{r_1^3 r_{12}^4}-
       \frac{76}{3 r_{12}^3} \right.
\nonumber \\ & &  \left. \left. \left. \quad \qquad
\qquad \qquad +\frac{2}{r_1 r_{12}^2}-
	\frac{6 r_2^2}{r_1^3 r_{12}^2}\right)\right)
\right]+O(8)+1 \leftrightarrow 2     \ , \label{7.2a}
\\ \nonumber\\
g_{0i} & = &
-4 \frac{G m_1}{c^3 r_1} v_1^i+
\frac{1}{c^5} \left[
n_1^i \left(-\frac{G^2 m_1^2}{r_1^2} (n_1v_1)+
		\frac{G^2 m_1 m_2}{S^2} \left(-16 (n_{12}v_1)+
	12 (n_{12}v_2)
			\right. \right. \right. \nonumber \\ & & \left.
			\frac{}{}
			\left. \qquad \qquad \qquad
	-16 (n_2v_1) +12 (n_2v_2) \right) \right) \nonumber \\
			& & \qquad \qquad +
n_{12}^i G^2 m_1 m_2 \left(-6 (n_{12}v_{12}) \frac{r_1}{r_{12}^3}-
	4 (n_1v_1)\frac{1}{r_{12}^2}+12 (n_1v_1)\frac{1}{S^2} \right.
			\nonumber \\ & & \left. \qquad \qquad \qquad \qquad
			\qquad \qquad \qquad-
	16 (n_1v_2)\frac{1}{S^2}+
	4 (n_{12}v_1)\frac{1}{S} \left(\frac{1}{S}+
		\frac{1}{r_{12}}\right) \right) \nonumber \\ & &
			\qquad \qquad +
v_1^i \left(\frac{G m_1}{r_1} \left(2 (n_1v_1)^2-4 v_1^2\right)+
		\frac{G^2 m_1^2}{r_1^2}+
	G^2 m_1 m_2 \left(\frac{3 r_1}{r_{12}^3}-\frac{2 r_2}{r_{12}^3}\right)
			\right. \nonumber \\ & & \left. \left.
			\qquad \qquad \qquad \qquad  +
	G^2 m_1 m_2  \left(-\frac{r_2^2}{r_1 r_{12}^3}  -\frac{3}{r_1 r_{12}}+
  		\frac{8}{r_2 r_{12}}-\frac{4}{r_{12} S}\right)\right)\right]
			\nonumber \\ &+&
\frac{1}{c^6} \left[
n_{12}^i \left(G^2 m_1 m_2 \left(-10 (n_{12}v_{12})^2 \frac{r_1^2}{r_{12}^4}-
	12 (n_{12}v_{12}) (n_1v_1) \frac{r_1}{r_{12}^3}+
	2 v_{12}^2 \frac{r_1^2}{r_{12}^4}-4 \frac{v_1^2}{r_{12}^2}\right)
			\right. \right. \nonumber \\ & & \left. \qquad \quad+
        G^3 m_1^2 m_2 \left(\frac{2 r_1^2}{3 r_{12}^5}-
		\frac{2 r_2^2}{3 r_{12}^5}-
		\frac{2}{r_{12}^3}\right)\right)
			\nonumber \\
			& & \qquad +
v_1^i \frac{G^2 m_1 m_2}{r_{12}^3} \left(\frac{16 (n_1v_{12}) r_1}{3}-
4 (n_{12}v_2) r_{12} \right)
			\nonumber \\ & & \qquad \left. +
v_{12}^i \frac{G^2 m_1 m_2}{r_{12}^2} \left(-2 (n_{12}v_1)+
	6 (n_{12}v_{12}) \frac{r_1^2}{r_{12}^2}\right) \right]+O(7)+
1 \leftrightarrow 2             \ , \label{7.2b}
\\ \nonumber\\
g_{ij} -\delta_{ij} & = &
2 \frac{G m_1}{c^2 r_1} \delta^{ij}+
\frac{1}{c^4} \left[
\delta^{ij} \left(
		-\frac{G m_1}{r_1} (n_1v_1)^2+\frac{G^2 m_1^2}{r_1^2}
			\right. \right.
			\nonumber \\ & & \left. \qquad \qquad
			\qquad \qquad   +
		G^2 m_1 m_2 \left(\frac{2}{r_1 r_2}-\frac{r_1}{2 r_{12}^3}+
		\frac{r_1^2}{2 r_2 r_{12}^3}-\frac{5}{2 r_1 r_{12}}+
  		\frac{4}{r_{12} S}\right)\right)
			\nonumber \\ & & \qquad \qquad \qquad +
4 \frac{G m_1}{r_1} v_1^i v_1^j+
\frac{G^2 m_1^2}{r_1^2} n_1^i n_1^j -
4 G^2 m_1 m_2 n_{12}^i n_{12}^j \left(\frac{1}{S^2}+
		\frac{1}{r_{12} S}\right)
			\nonumber \\ & & \left.
			\qquad \qquad \qquad+
\frac{4 G^2 m_1 m_2}{S^2} \left(n_1^{(i} n_2^{j)}+
	2n_1^{(i} n_{12}^{j)} \right)\right]
			\nonumber \\ &+&
\frac{G^2 m_1 m_2}{c^5r_{12}^2} \left(-\frac{2}{3}(n_{12}v_{12})\delta^{ij}-
6 (n_{12}v_{12}) n_{12}^i n_{12}^j+ 8 n_{12}^{(i} v_{12}^{j)}\right)+O(6)+
1 \leftrightarrow 2    \ , \label{7.2c}
\end{eqnarray}
\end{mathletters}
(where we recall $S=r_1+r_2+r_{12}$). The harmonic-coordinate conditions
satisfied by (\ref{7.2}) read, with this order of approximation:
$\partial_\nu [(-g)^{1/2} g^{0\nu}]=O(7)$ and
$\partial_\nu [(-g)^{1/2} g^{i\nu}]=O(6)$.
```

The paper's footnote N12, attached to Eq. (7.2):

```latex
\bibitem{N12}Note that many of the ``odd'' terms in the metric, which are
associated with radiation-reaction, are only functions of time $t$ or depend
smoothly on the field point ${\bf x}$ (through for instance the product
$r_1{\bf n}_1={\bf x}-{\bf y}_1$). However, when presented into a form which
shows this fact more explicitly (expressing all the terms with the help of
the vectors ${\bf r}_1$ and ${\bf r}_{12}$) the metric is not simpler.
```

## Section VII, Eq. (7.6): the metric regularized at body 1

```latex
\begin{equation}
(g_{\mu\nu})_1(t) = g_{\mu\nu}[{\bf y}_1(t); {\bf y}_1(t),{\bf y}_2(t);
{\bf v}_1(t),{\bf v}_2(t)] \ , \label{7.5}
\end{equation}
```

```latex
\begin{mathletters}
\label{7.6}
\begin{eqnarray}
(g_{00})_1 &=&-1+2 \frac{G m_2}{c^2r_{12}}
+\frac{G m_2}{c^4 r_{12}} \left(4 v_2^2-(n_{12}v_2)^2-3\frac{G m_1}{r_{12}}
- 2 \frac{G m_2}{r_{12}} \right) \nonumber \\
&+&\frac{8G^2 m_1 m_2}{3c^5 r_{12}^2}(n_{12}v_{12})
+ \frac{G m_2}{c^6 r_{12}} \left( \frac{3}{4} (n_{12}v_2)^4-
3 (n_{12}v_2)^2 v_2^2 + 4 v_2^4 \right) \nonumber \\ &+&
\frac{G^2 m_1 m_2}{c^6 r_{12}^2} \left(-\frac{87}{4} (n_{12}v_1)^2 +
\frac{47}{2}(n_{12}v_1) (n_{12}v_2)-
  \frac{55}{4} (n_{12}v_2)^2 + \frac{23}{4} v_1^2 -
  \frac{39}{2} (v_1v_2) \right) \nonumber \\
&+& \frac{47}{4} \frac{G^2 m_1 m_2}{c^6 r_{12}^2} v_2^2+
\frac{G m_2}{c^6 r_{12}} \left\{ \frac{G m_2}{r_{12}}
\left[3 (n_{12}v_2)^2 - v_2^2 \right]
-\frac{G^2 m_1^2}{r_{12}^2} + \frac{17}{2} \frac{G^2 m_1 m_2}{r_{12}^2} +
  2 \frac{G^2 m_2^2}{r_{12}^2}  \right\} \nonumber \\
&+& \frac{G^2 m_1 m_2}{c^7 r_{12}^2} \left\{ -20 (n_{12}v_1)^3 +
40 (n_{12}v_1)^2 (n_{12}v_2) - 36 (n_{12}v_1) (n_{12}v_2)^2 +
  16 (n_{12}v_2)^3 \right. \nonumber
\\ & & \qquad \qquad+\frac{296}{15} (n_{12}v_1) v_1^2  -
\frac{116}{15} (n_{12}v_2) v_1^2 -
  \frac{104}{5} (n_{12}v_1) (v_1v_2) + \frac{232}{15} (n_{12}v_2) (v_1v_2)
\nonumber \\ & & \qquad \qquad+
 \frac{56}{15} (n_{12}v_1) v_2^2 - \frac{52}{5} (n_{12}v_2) v_2^2
+\frac{G m_1}{r_{12}} \left(-\frac{64}{5} (n_{12}v_1)+
\frac{104}{5} (n_{12}v_2) \right)
\nonumber \\ & & \qquad \qquad \left.
+\frac{G m_2}{r_{12}} \left(-\frac{144}{5}(n_{12}v_1)+
\frac{392}{15}(n_{12}v_2)\right)
\right\}+O(8)   \ , \label{7.6a}\\
(g_{0i})_1 &=&-4 \frac{G m_2}{c^3 r_{12}} v_2^i
+\frac{G m_2}{c^5 r_{12}}\left\{n_{12}^i
\left[ \frac{G m_1}{r_{12}} \left(10 (n_{12}v_1) + 2 (n_{12}v_2)\right)
-\frac{G m_2}{r_{12}} (n_{12}v_2)\right]
\right\} \nonumber \\
&+& \frac{G m_2}{c^5 r_{12}} \left\{ 4 \frac{G m_1}{r_{12}} v_1^i
+v_2^i \left(2 (n_{12}v_2)^2 -
4 v_2^2- 2 \frac{G m_1}{r_{12}} +
\frac{G m_2}{r_{12}} \right) \right\}  \nonumber \\
&+&\frac{G^2 m_1 m_2}{c^6 r_{12}^2} \left\{n_{12}^i \left(10 (n_{12}v_1)^2 -
8 (n_{12}v_1) (n_{12}v_2)
- 2 (n_{12}v_2)^2 - 6 v_1^2  \right. \right.
\nonumber \\ & & \left. \qquad \qquad \qquad + 4 (v_1v_2) + 2 v_2^2-
\frac{8}{3} \frac{G m_1}{r_{12}} + \frac{4}{3} \frac{G m_2}{r_{12}} \right)
\nonumber \\ & &   \left. \qquad \qquad -8 (n_{12}v_1) v_1^i
+v_2^i \left(\frac{20}{3} (n_{12}v_1) +
\frac{4}{3} (n_{12}v_2)\right) \right\} +O(7) \ , \label{7.6b} \\
(g_{ij})_1 &=& \delta^{ij}
+2 \frac{G m_2}{c^2 r_{12}} \delta^{ij}
+\frac{G m_2}{c^4 r_{12}} \delta^{ij} \left(-(n_{12}v_2)^2 +
 \frac{G m_1}{r_{12}} + \frac{G m_2}{r_{12}} \right) \nonumber \\
&+&\frac{G m_2}{c^4 r_{12}} \left\{n_{12}^{ij}
\left(-8 \frac{G m_1}{r_{12}} +
\frac{G m_2}{r_{12}}\right) + 4 v_2^{ij}\right\}
-\frac{4 G^2 m_1 m_2}{3c^5 r_{12}^2}
\delta^{ij} (n_{12}v_{12}) \nonumber \\
&+&\frac{G^2 m_1 m_2}{c^5 r_{12}^2} \left\{
-12 n_{12}^{ij}(n_{12}v_{12})+
16 n_{12}^{(i} v_{12}^{j)}\right\} +O(6) \ . \label{7.6c}
\end{eqnarray}
\end{mathletters}
```

## Living Review (Blanchet 2014), matching equations

Section 5.3 of the Living Review ("The 3.5PN metric for general matter systems") gives the same
potentials. Differences from the 1998 paper:

- The retarded integral is written $\Box^{-1}_\mathrm{ret}$ in place of $\Box^{-1}_R$, and
  remainders as $\mathcal{O}(1/c^n)$ in place of $O(n)$.
- The metric, Eq. (144), goes to 3.5PN and so carries the extra potentials $\hat T$,
  $\hat Y_i$, $\hat Z_{ij}$ at order $1/c^8$, $1/c^7$, $1/c^6$. Cut at the 1998 paper's order,
  it is Eq. (2.5) above.
- The definitions of $\sigma$, $V$, $V_i$, $\hat W_{ij}$, $\hat R_i$, $\hat X$ are the same as
  (2.1), (2.3), (2.4); only the order of the terms in the source of $\hat X$ differs.
- The Living Review does **not** list the two-body potentials (B1) or their values at body 1
  (B2). It gives the explicit two-body metric instead, Eq. (241), written as
  $g_{\alpha\beta} = \eta_{\alpha\beta} + k_{\alpha\beta}$ and credited to the 1998 paper.

Equation numbers below are those printed in the arXiv PDF of version 2 (2014). The 2024 revision
of the review (arXiv version 6) has the same definitions under the same labels (`gmunu3PN`,
`sigmadef`, `pot1PN`, `pot2PN`); its equation numbers are different and were not looked up.

The 3.5PN metric, Eqs. (144a)-(144c):

```latex
\begin{subequations}\label{gmunu3PN}
  \begin{align}
  g_{00} &= - 1 + \frac{2}{c^2} V - \frac{2}{c^4} V^2 + \frac{8}{c^6}
  \left( \hat{X} + V_i V_i + \frac{V^3}{6} \right) \nonumber\\ &+
  \frac{32}{c^8} \left( \hat{T} - \frac{1}{2} V \hat{X} + \hat{R}_i
  V_i - \frac{1}{2} V V_i V_i - \frac{1}{48} V^4 \right) + \mathcal{O}
  \left( \frac{1}{c^{10}} \right)\,, \\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g_{0i} &= - \frac{4}{c^3} V_i - \frac{8}{c^5} \hat{R}_i -
\frac{16}{c^7} \left( \hat{Y}_i + \frac{1}{2} \hat{W}_{ij} V_j +
\frac{1}{2} V^2 V_i \right) + \mathcal{O}\left(\frac{1}{c^9}\right)\,,
\\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g_{ij} &= \delta_{ij} \left[ 1 + \frac{2}{c^2} V + \frac{2}{c^4} V^2 +
  \frac{8}{c^6} \left(\hat{X} + V_k V_k + \frac{V^3}{6} \right)
  \right] + \frac{4}{c^4} \hat{W}_{ij} \nonumber\\ &+ \frac{16}{c^6}
\left( \hat{Z}_{ij} + \frac{1}{2} V \hat{W}_{ij} - V_i V_j \right) +
     \mathcal{O}\left(\frac{1}{c^8} \right)\,.
  \end{align}
\end{subequations}
```

Source densities, Eqs. (145a)-(145c):

```latex
\begin{subequations}\label{sigmadef}
  \begin{align}
    \sigma &= \frac{T^{00}+T^{ii}}{c^2}\,, \\ \sigma_{i} &=
    \frac{T^{0i}}{c}\,, \\ \sigma_{ij} &= T^{ij}\,.
  \end{align}
\end{subequations}
```

$V$ and $V_i$, Eqs. (146a), (146b):

```latex
\begin{subequations}\label{pot1PN}
  \begin{align}
    V &= \Box^{-1}_\mathrm{ret} \left[ -4 \pi G \sigma \right]\,,\\
    V_i &= \Box^{-1}_\mathrm{ret} \left[ -4 \pi G \sigma_i \right]\,.
  \end{align}
\end{subequations}
```

$\hat X$, $\hat R_i$, $\hat W_{ij}$, Eqs. (147a)-(147c):

```latex
\begin{subequations}\label{pot2PN}
  \begin{align}
    \hat{X} &= \Box^{-1}_\mathrm{ret} \left[ - 4 \pi G
      V \sigma_{ii}+\hat{W}_{ij} \partial_{ij} V+2 V_i \partial_t
      \partial_i V+ V \partial_t^2 V + \frac{3}{2} (\partial_t V)^2 -
      2 \partial_i V_j \partial_j V_i \right]\,, \\
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\hat{R}_i &= \Box^{-1}_\mathrm{ret} \left[ -4 \pi G (V \sigma_i-V_i
  \sigma)-2 \partial_k V \partial_i V_k-\frac{3}{2} \partial_t V
  \partial_i V \right], \\ \hat{W}_{ij} &= \Box^{-1}_\mathrm{ret}
\left[ -4 \pi G (\sigma_{ij}-\delta_{ij} \sigma_{kk})- \partial_i V
  \partial_j V \right]\,.
  \end{align}
\end{subequations}
```

Section 7.5 ("The 2.5PN metric in the near zone") also records the elementary functions the two-body potentials are built from. The Poisson equation and its solution $g = \ln S$, Eqs. (237), (238):

```latex
\begin{equation}
\Delta g = \frac{1}{r_1 r_2} \,,  \label{Laplaceg}
\end{equation}
```

```latex
\begin{subequations}\label{lnS}
\begin{align}
g &= \ln S \,, \\ S &\equiv r_1+r_2+r_{12} \,,
\end{align}
\end{subequations}
```

The harder Poisson equations and their closed-form solutions $K_1$, $H_1$, Eqs. (239), (240). Here ${}_\mathrm{a}\partial_i$ is the derivative with respect to the source point $y_\mathrm{a}^i$ and $\Delta_\mathrm{a}$ the Laplacian with respect to it:

```latex
\begin{subequations}\label{eqK1H1}
\begin{align}
\Delta K_1 &= 2\,\partial_{i}\partial_{j}\left(\frac{1}{r_2}\right)
\partial_{i}\partial_{j}\ln r_1 \,, \label{eqK1} \\ \Delta H_1 &=
2\,\partial_{i}\partial_{j}\left(\frac{1}{r_1}\right)
\mathop{\partial}_{1}{}_{\!\!i}\!\mathop{\partial}_{2}{}_{\!\!j}g \,,
\label{eqH1}
\end{align}
\end{subequations}
```

```latex
\begin{subequations}\label{K1H1}
\begin{align}
K_1 &= \left(\frac{1}{2}\Delta-\Delta_1\right) \left[{\ln r_1\over
    r_2}\right] + \frac{1}{2}\Delta_2\left[{\ln r_{12}\over
    r_2}\right] + {r_2\over 2r_{12}^2r_1^2}+{1\over r_{12}^2r_2}
\,, \label{K1} \\ H_1 &= \frac{1}{2}\Delta_1\left[{g\over r_1}+{\ln
    r_1\over r_{12}}- \Delta_1\left({r_1+r_{12}\over 2}g
  \right)\right] - {r_2\over 2r_1^2r_{12}^2}\nonumber\\ &+
\partial_i\!\mathop{\partial}_{2}{}_{\!\!i}\left[{\ln r_{12} \over
    r_1}+ {\ln r_1\over 2r_{12}}\right] +
\frac{1}{2}\Delta_2\left[{\ln r_{12}\over r_{1}}\right] \,, \label{H1}
\end{align}
\end{subequations}
```

The explicit near-zone metric, Eq. (241) of the Living Review, is not copied here because it
repeats Eq. (7.2) above. This was checked by script: with spacing and alignment commands
removed, the two LaTeX sources agree term for term. The only differences are that the Living
Review writes $k_{00}$, $k_{0i}$, $k_{ij}$ for $g_{00}+1$, $g_{0i}$, $g_{ij}-\delta_{ij}$,
typesets the leading $g_{0i}$ term as $-\frac{4Gm_1}{c^3 r_1}v_1^i$, and writes the remainders
as $\mathcal{O}(1/c^8)$, $\mathcal{O}(1/c^7)$, $\mathcal{O}(1/c^6)$.

## Known vs new

Everything here is known, published mathematics, copied from the two sources cited at the top.
Nothing is new and nothing is ours.

## Assumptions (of the source, as it states them)

- General relativity, harmonic (De Donder) coordinates, post-Newtonian expansion valid in the
  near zone of a slowly moving, weakly stressed, self-gravitating source.
- The two bodies are point masses, with their infinite self-field removed by the Hadamard
  partie finie, Eqs. (3.5)-(3.7). The paper notes this choice of regularization is part of the
  physical model.
- Accelerations are replaced using the lower-order equations of motion (order reduction),
  Eq. (3.11).
- The bodies have no spin.
