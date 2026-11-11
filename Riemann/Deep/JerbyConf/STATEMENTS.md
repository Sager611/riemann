# Jerby, arXiv 2610.21952 ("On the Sampling of Certain Dirichlet and Prime Sums at Gram Points and their Viscosity"): finite statements, exact quotes

Source: `https://arxiv.org/e-print/2610.21952` (gzip tar of 15 Oct 2026; `Clouds.tex`, 2,359 lines). Line numbers below are `Clouds.tex` lines; quotes are verbatim TeX (`nnnn| ` prefix = line number). Data: Zenodo record 23390338 (`Zonedo.zip`, 166.6 MB): `GramPoints_CORE_1_250000.csv` (n, g, Z, Z', mu_exact, class) and `DH_GramPoints_12500_Zenodo_Database.xlsx` (2 ≤ n ≤ 12501).

## 0. Definitions

Theta, Gram points, Gram's law, good/bad/isolated:

```
  92| \[
  93| \theta(t)
  94| :=
  95| \Im\log\Gamma\!\left(\frac14+\frac{it}{2}\right)
  96| -\frac{t}{2}\log\pi .
  97| \]
  98| By Stirling's formula,
  99| \[
 100| \theta(t)
 101| =
 102| \frac{t}{2}\log\!\left(\frac{t}{2\pi}\right)
 103| -\frac{t}{2}
 104| -\frac{\pi}{8}
 105| +O\!\left(\frac1t\right),
 106| \qquad t\to\infty.
 107| \]
 108| The \(n\)-th Gram point \(g_n\) is defined by
 109| \[

 116| usually obey the rule
 117| \begin{equation}
 118| \label{eq:Gram}
 119|         (-1)^n Z(g_n)>0.

 128| However, Gram's law is not always true.  The first violations were discovered
 129| by Hutchinson in 1925, who found that the Gram points \(g_{126}\) and \(g_{137}\)
 130| fail the expected sign rule \cite{Hutchinson}.  We hence refer to \(g_n\) as a
 131| \emph{good Gram point} if it satisfies Gram's law \eqref{eq:Gram}, and as a
 132| \emph{bad Gram point} if it fails. We refer to a bad Gram point $g_n$ as \emph{isolated} if its neighbours $g_{n \pm 1}$ are good. A modern experimental study of the relative frequency of good and bad Gram
 133| points was carried out by Shanker~\cite{ShankerGramRatio}.
```
AFE length and viscosity:
```
 236| where \(        N(t)
 237|         =
 238|         \left\lfloor
 239|         \sqrt{\frac{t}{2\pi}}
 240|         \right\rfloor

 312| collision of neighboring zeros, rather than an accidental numerical pattern. To quantify this mechanism, we introduced the local
 313| invariant
 314| \begin{equation}
 315| \label{eq:visc}
 316|         \mu(g_n)
 317|         :=
 318|         \left|
 319|         \frac{Z'(g_n)}{Z(g_n)}
 320|         \right|,
 321| \end{equation}
 322| which we refer to as the \emph{viscosity} of the Gram point \(g_n\).
```
Variation space, core zeros, confined pair, discriminant:
```
 448|  provides a setting in
 449| which such a dynamical interpretation becomes possible. Consider the
 450| parameter space
 451| \[
 452|     \mathcal Z_N
 453|     =
 454|     \left\{
 455|     Z_N(t;\mathbf a)
 456|     =
 457|     Z_0(t)
 458|     +
 459|     \sum_{k=1}^{N}
 460|     \frac{a_k}{\sqrt{k+1}}
 461|     \cos\!\bigl(\theta(t)-t\log(k+1)\bigr)
 462|     :
 463|     \mathbf a=(a_1,\ldots,a_N)\in\mathbb R^N
 464|     \right\},
 465| \]
 466| where
 467| \[

 470| is the core function. The origin \( \mathbf a=0 \) therefore corresponds to
 471| the core, whereas \( \mathbf a=\mathbf 1=(1,\ldots,1) \) gives the arithmetic
 472| Hardy section, which satisfies the approximation
 473| \[ Z(t)= Z_N(t) + O \left ( t^{-\frac{1}{4}} \right) \] in the range $2N \leq t \leq 2N+2 $. We consider a continuous path
 474| \[
 475|         \gamma:[0,1]\longrightarrow\mathbb R^N,
 476|         \qquad
 477|         \gamma(0)=0,
 478|         \qquad
 479|         \gamma(1)=\mathbf 1,
 480| \]
 481| and follow the zeros of
 482| \(
 483|         Z_N(t;\gamma(a))
 484| \)
 485| as \( a \) varies. Associated with a Gram point \( g_n \), the core function
 486| has two neighboring zeros

 503| \).
 504| 
 505| \begin{definition}
 506| We will refer to the pair \( t_-(a),t_+(a) \) associated with an isolated Gram point \( g_n \) as a \emph{confined} pair of zeros if, throughout the deformation \( 0\leq a\leq1 \), these
 507| two branches remain the only zeros of \( Z(t; \gamma(a) ) \), counted with multiplicity,
 508| in the vertical Gram strip
 509| \[
 510|        g_{n-1} < Re(t) \leq g_{n+1}.
 511| \]

 527|         :=
 528|         Z_N\!\left(g_n(\mathbf a);\mathbf a\right).
 529| \]
 530| A collision of \( t_-(a) \) and \( t_+(a) \) occurs precisely when
 531| \(
 532|         \Delta_n(\gamma(a))=0.
 533| \)
 534| Consequently, for a confined pair the question of whether the two zeros
 535| remain real is reduced locally to the problem of avoiding the discriminant
 536| hypersurface
```
Notes on the definitions (ours):
- *Confined* (l. 505–511) is path-dependent: any continuous path γ from 0 to 1, and the count is over complex zeros in the strip, for every a. At a = 1 it forces exactly two zeros (with multiplicity) in the block. The paper's only operational use (l. 630–636 below) identifies "confined" with "Rosser's rule holds for the length-2 block". We use that reading: **confined ⟺ the block (g_{n−1}, g_{n+1}] holds exactly two zeros of Z**. Because Z(g_{n−1}), Z(g_n), Z(g_{n+1}) share one sign at an isolated bad point, each Gram interval holds an even number of zeros, so two zeros means a pair inside one Gram interval and none in the other. RH is verified far beyond our heights (Platt–Trudgian, to 3·10^12), so every zero in the strip is real.
- l. 1647–1652 call the same-side placement of the two zeros "Empirically" observed; it is forced by the sign argument above.
- l. 473: "Z(t) = Z_N(t) + O(t^{-1/4}) in the range 2N ≤ t ≤ 2N+2" cannot hold as printed: at a = 1, Z_N(t;1) = Σ_{k≤N+1} k^{-1/2} cos(θ(t) − t log k) is half the Riemann–Siegel main sum, and the Riemann–Siegel range is 2πM² ≤ t < 2π(M+1)² with M = N+1. Zero sets are unaffected.
- Normalisation used below: s_n μ(g_n) with s_n = 2π/log(g_n/2π) (the paper's Δ_n, l. 2016–2021) equals π μ(g_n)/θ′(g_n) up to O(g_n^{-2}); our scanner prints μ/θ′ ("norm").

## 1. Finite statements (6)

### F1. Repulsion bound (l. 323–326); data (l. 338–353, 637–648)
```
 323| Initial experiments in \cite{JerbyDiscriminant} suggest that isolated bad Gram points satisfy the following repulsion property
 324| \begin{equation}
 325| \mu(g_n)>4,
 326| \end{equation} to height of up to $g_n$ of index $n \sim 2.5 \cdot 10^{7}$.

 338| \end{figure}
 339| 
 340| 
 341| The figure is based on the database of the first \(250{,}000\) Gram points,
 342| containing \(220{,}373\) good Gram points and \(24{,}269\) isolated bad Gram
 343| points. The distribution associated with good Gram points is
 344| broad and approximately symmetric on the logarithmic scale, with substantial
 345| tails extending in both directions.  By contrast, the distribution associated
 346| with isolated bad Gram points is strongly asymmetric. It has a sharp lower edge
 347| immediately above the repulsion bound \(\mu=4\), rises rapidly to a pronounced peak, and then
 348| decays through a long right-hand tail. Thus the numerical evidence points not only to a lower repulsion barrier,
 349| but to a more rigid statistical organization of isolated bad Gram points
 350| according to the size of \(\mu(g_n)\). In the present
 351| data set, the smallest observed value at an isolated bad Gram point is
 352| \[
 353|         \mu(g_n)\approx 4.26178.

 637| 
 638| It should be noted that in our previous work~\cite{JerbyDiscriminant},
 639| the viscosity bound
 640| \(
 641|         \mu(g_n)>4
 642| \)
 643| was investigated at substantially larger indices, including a systematic
 644| computation of all bad Gram points up to the range
 645| \( n\leq 2.5 \times10^{7}\).
 646| The numerical evidence revealed a sharp distinction between isolated and
 647| non-isolated bad Gram points. The bound was seen to persist for the isolated
 648| configurations, whereas bad Gram points belonging to longer Gram blocks
 649| could have viscosity far below \(4\).  In particular, the low-viscosity
 650| ``corrupt'' points observed there were non-isolated. An extreme example is the non-isolated
```
Range tested by the author: n ≤ 250,000 (t < 170,572) in this paper; "all bad Gram points up to ... n ≤ 2.5·10^7" in 2406.06548. **Status: false** (jerbyfals, 11 Nov 2026: three Arb-certified isolated bad points with μ < 4, all zero-free blocks; `Riemann/Deep/JerbyFals/REPORT.md`).

### F2. Confined-pair law (l. 616–663, Remark "Confined zeros and Rosser's law")
```
 616| \begin{remark}[Confined zeros and Rosser's law]
 617| The confinement condition is stronger than requiring \( g_n \) to be an
 618| isolated bad Gram point. Recall that a \emph{Gram block} of length \( N \)
 619| is a sequence
 620| \[
 621| g_m,g_{m+1},\ldots,g_{m+N},
 622| \]
 623| whose endpoints are good and whose intermediate Gram points are bad
 624| \cite{Edwards,RosserYoheSchoenfeld}. Thus an isolated bad Gram point
 625| \( g_n \) is the unique interior point of a Gram block of length two,
 626| \(
 627| (g_{n-1},g_{n+1}]
 628| \). Rosser's rule asserts that a Gram block of length \( N \) contains exactly
 629| \( N \) zeros of \( Z(t) \) on the real axis. Hence a length-two block
 630| normally contains two zeros. However, Rosser's rule is not universal, and
 631| Trudgian proved that it fails a positive proportion of the time
 632| \cite{Trudgian}. Thus isolated badness does not, in general, guarantee a
 633| two-zero local configuration. It should be noted, however, that among the first \(250{,}000\) Gram points
 634| in the data set considered, no violations of Rosser's rule occur.
 635| Consequently, every isolated bad Gram point encountered in this range is
 636| associated with a confined pair.
 637| 
 638| It should be noted that in our previous work~\cite{JerbyDiscriminant},
 639| the viscosity bound
 640| \(
 641|         \mu(g_n)>4
 642| \)
 643| was investigated at substantially larger indices, including a systematic
 644| computation of all bad Gram points up to the range
 645| \( n\leq 2.5 \times10^{7}\).
 646| The numerical evidence revealed a sharp distinction between isolated and
 647| non-isolated bad Gram points. The bound was seen to persist for the isolated
 648| configurations, whereas bad Gram points belonging to longer Gram blocks
 649| could have viscosity far below \(4\).  In particular, the low-viscosity
 650| ``corrupt'' points observed there were non-isolated. An extreme example is the non-isolated
 651| \(g_{9807962}\), for which
 652| \(\mu(g_{9807962})=0.0750883\).
 653| Finally, we note that the data suggest more than the lower bound
 654| \(\mu(g_n)>4\) for isolated bad Gram points associated with a confined pair
 655| of zeros.  They also suggest the asymptotic relation
 656| \begin{equation}
 657| \label{eq:liminf_scaled_mu}
 658| \liminf_{\substack{n\to\infty\\
 659| g_n\ \mathrm{isolated\ bad,\ confined}}}
 660| \frac{2\pi}{\log(g_n/2\pi)}\,\mu(g_n)
 661| =
 662| 2.
 663| \end{equation}
```
Range: the 250,000-point database (n ≤ 250,000, t < 170,572); the author gives no table for the liminf. A liminf law; tested by its trend.

### F3. Trunk-slope law (l. 1762–1770, 1822–1864, Remark "Partial-sum trunks near extremal configurations")
```
1762| \begin{remark}[Partial-sum trunks near extremal configurations] For \(t\in\mathbb R\), define the partial AFE sums
1763| \[
1764|         S_j(t)
1765|         :=
1766|         \sum_{k=1}^{j}
1767|         \frac{e^{it\log k}}{\sqrt{k}},
1768|         \qquad
1769|         1\le j\le N(t).
1770| \]

1822|     points \(n=8569,45924,64256,39324,96237\).  The colors are as in
1823|     Figure~\ref{fig:generic_trunks}.}
1824|     \label{fig:low_visc_trunks}
1825| \end{figure}
1826| 
1827| In these examples, one neighboring trunk terminates near a comparatively
1828| small endpoint, while the opposite neighboring trunk reaches a much larger
1829| endpoint.  Writing
1830| \[
1831|         A_{\mathrm{act}}
1832|         :=
1833|         \max\{|Z(g_{n-1})|,|Z(g_{n+1})|\},
1834| \]
1835| the observed coherent configurations suggest the rough slope estimate
1836| \[
1837|         |Z'(g_n)|
1838|         \asymp
1839|         \frac{A_{\mathrm{act}}}{g_{n+1}-g_n}
1840|         \sim
1841|         \frac{\theta'(g_n)}{\pi}A_{\mathrm{act}}.
1842| \]
1843| Consequently,
1844| \[
1845|         \mu(g_n)
1846|         =
1847|         \left|\frac{Z'(g_n)}{Z(g_n)}\right|
1848|         \asymp
1849|         \frac{\theta'(g_n)}{\pi}
1850|         \frac{A_{\mathrm{act}}}{|Z(g_n)|}.
1851| \]
1852| This gives an intuitive explanation for why very small viscosity is difficult
1853| to achieve and how the bound arises in such coherent configurations. The same phase organization that
1854| produces cancellation at the central bad point tends, after one Gram step, to
1855| produce a compensating neighboring excursion, and this excursion is measured
1856| on the derivative scale by division by the Gram spacing.
1857| 
1858| We emphasize that this is only a qualitative mechanism.  We do not expect
1859| every extremal point to arise from such a coherent three-trunk
1860| configuration.  The figures of such coherent configurations should therefore be read as evidence for one
1861| visible source of the viscosity barrier, not as a classification of all
1862| near-extremal isolated bad Gram points.
1863| \end{remark}
```
Range: the five lowest-viscosity points n = 8569, 45924, 64256, 39324, 96237 (Fig. low-visc-trunks) and five generic points; "≍" has no stated constants and the author calls it qualitative (l. 1858–1863). Testable form: the ratio ρ_n = |Z′(g_n)|(g_{n+1} − g_n)/A_act stays within fixed positive bounds at low-viscosity isolated bad points.

### F4. Lower-scale law (l. 2010–2051)
```
2010| denotes the local mean zero spacing, and if the two distances from \(g_n\) are
2011| written in normalized form as
2012| \[
2013|         g_n-\gamma_1 = a_1\Delta_n,
2014|         \qquad
2015|         g_n-\gamma_2 = a_2\Delta_n,
2016| \]
2017| then the contribution of this pair to the normalized viscosity is of the form
2018| \[
2019|         \frac{2\pi}{\log(g_n/2\pi)}
2020|         \left|
2021|         \frac{1}{g_n-\gamma_1}
2022|         +
2023|         \frac{1}{g_n-\gamma_2}
2024|         \right|
2025|         =
2026|         \left|
2027|         \frac{1}{a_1}
2028|         +
2029|         \frac{1}{a_2}
2030|         \right|.
2031| \]
2032| For a confined pair in a single neighboring Gram interval one expects
2033| \(0<a_1,a_2\leq 1\), and therefore this contribution is at least \(2\), with
2034| equality corresponding to the extremal case in which both zeros are pushed
2035| toward the far endpoint of that interval.  This gives the basic intuition
2036| behind the normalized lower scale
2037| \[
2038|         \frac{2\pi}{\log(g_n/2\pi)}\mu(g_n)
2039|         \gtrsim 2
2040| \]
2041| for isolated bad Gram points of confined type.
2042| 
2043| 
2044| It is important to note that this argument should not be read as a proof of the observed lower envelope.
2045| It keeps only the contribution of the confined zero pair and does not yet
2046| account for the nearby zeros outside the Gram block
2047| \([g_{n-1},g_{n+1}]\), nor for the more distant zeros and slowly varying
2048| background terms in the logarithmic derivative.  Nevertheless, it explains
2049| why the viscosity is a natural quantity for detecting the local zero
2050| configuration, and why extremal low-viscosity bad points should correspond to
2051| highly coherent arrangements of the confined zero pair.
```
Range: heuristic, no data stated. Two testable parts: (a) the pair term |1/a_1 + 1/a_2| ≥ 2 given 0 < a_1, a_2 ≤ 1; (b) s_n μ(g_n) ≳ 2 at confined isolated bad points, i.e. the other zeros and background do not pull the value far below the pair term.

### F5. Davenport–Heilbronn collision parameter (l. 678–701, 744–800, 816–831)
```
 678| \[
 679|         D(s)
 680|         =
 681|         \frac{1-i\kappa}{2}L(s,\chi_{5,2})
 682|         +
 683|         \frac{1+i\kappa}{2}L(s,\overline{\chi}_{5,2}),
 684| \]
 685| where \( \chi_{5,2} \) denotes the non-real primitive Dirichlet character modulo \(5\), chosen so that \( \chi_{5,2}(2)=i \) and
 686| \[
 687|         \kappa
 688|         =
 689|         \frac{\sqrt{10-2\sqrt5}-2}{\sqrt5-1}.
 690| \]
 691| Its completed form
 692| \[
 693|         \xi_{DH}(s)
 694|         =
 695|         \left(\frac{\pi}{5}\right)^{-s/2}
 696|         \Gamma\!\left(\frac{1+s}{2}\right)D(s)
 697| \]
 698| satisfies the symmetric functional equation
 699| \[
 700|         \xi_{DH}(s)=\xi_{DH}(1-s).
 701| \]

 716|         =
 717|         \frac12
 718|         \left[
 719|         \left(\frac{\pi}{5}\right)^{-s/2}
 720|         \Gamma\!\left(\frac{1+s}{2}\right)
 721|         +
 722|         \left(\frac{\pi}{5}\right)^{(s-1)/2}
 723|         \Gamma\!\left(\frac{2-s}{2}\right)
 724|         \right],
 725| \]
 726| and, after normalization on the critical line, gives a real core function
 727| \( Z_0^{DH}(t) \). As in the \( Z \)-function case, the arithmetic terms are
 728| then activated through a parameter space
 729| \[
 730|         \mathcal Z_N^{DH}
 731|         =
 732|         \bigl\{
 733|         Z_N^{DH}(t;\mathbf a)
 734|         \bigr\},
 735| \]
 736| with the core at \( \mathbf a=0 \) and the arithmetic Davenport--Heilbronn
 737| section at \( \mathbf a=\mathbf 1 \). Along the diagonal deformation
 738| \(\mathbf a=r\mathbf 1 \), for
 739| $0\leq r\leq1$,
 740| one obtains a family
 741| \(
 742|         Z_N^{DH}(t;r)
 743| \)
 744| and a continued extremum \( g_n^{DH}(r) \), exactly as for the Hardy

 762| To see this, define the Davenport--Heilbronn Gram discriminant by
 763| \[
 764|         \Delta_{44}^{DH}(r)
 765|         =
 766|         Z_N^{DH}
 767|         \bigl(
 768|         g_{44}^{DH}(r);r
 769|         \bigr).
 770| \]
 771| In striking contrast with the isolated bad Gram point \( g_{126} \) of the
 772| Hardy \( Z \)-function, the first-order approximation
 773| \(
 774|        Z_N^{DH}(g_{44}^{DH};r)
 775| \)
 776| tracks the corrected discriminant very closely. There is essentially no
 777| large nonlinear correction which restores the correct sign. Instead,
 778| \(\Delta_{44}^{DH}(r)\) itself reaches zero. In the numerical deformation
 779| shown below this occurs at
 780| \[
 781|         r_*
 782|         \approx
 783|         0.72346.
 784| \]
 785| At this parameter value,
 786| \[
 787|         \Delta_{44}^{DH}(r_*)=0,
 788| \]
 789| and equivalently
 790| \[
 791|         t_-^{DH}(r_*)
 792|         =
 793|         t_+^{DH}(r_*)
 794|         =
 795|         g_{44}^{DH}(r_*).
 796| \]
 797| The two real zero branches therefore collide at a multiple zero. Beyond the
 798| collision they cease to form two distinct real zeros and continue as the
 799| non-real pair responsible for the failure of the Davenport--Heilbronn
 800| analogue of the Riemann Hypothesis.
```
### F6. Davenport–Heilbronn viscosity at g_44 (l. 838–847)
```
 838| zeros into the complex plane.
 839| 
 840| It should be noted that the analogous viscosity in this
 841| Davenport--Heilbronn example is small: for the section used in
 842| Figure~\ref{fig:DH_collision} one obtains numerically
 843| \[
 844|         \mu_{DH}(g_{44}^{DH})
 845|         \approx 1.45,
 846| \]
 847| far below the conjectural threshold \(4\). Correspondingly, the transverse
```
F5 and F6 depend on the DH section Z_N^DH(t; a), which the paper does not write out: no N and no arithmetic terms are given (l. 727–744: "the arithmetic terms are then activated through a parameter space"). The Zenodo workbook gives g_44^DH, its unnormalised values ZDH(g_n) (≈ 10^{-29}, i.e. the completed ξ_DH(1/2 + it) including the Γ-factor decay) and mu_exact.

### Side numerics (data claims)
```
1866| \subsection{The distribution of $\mu(g_n^{DH})$ for isolated bad DH Gram points}
1867| 
1868| Recall that Fig.~\ref{fig:viscosity_distributions} revealed a striking
1869| separation between the viscosity distributions of good and isolated bad Gram
1870| points for $Z(t)$.  While the good-point distribution is broad and extends
1871| across a wide range of values, the isolated-bad distribution develops a sharp
1872| lower edge immediately above
1873| \(
1874| \mu=4.
1875| \)
1876| In the data considered there, no isolated bad Gram point associated with a
1877| confined pair was observed below this threshold, motivating the repulsion Conjecture.
1878| 

1884| \begin{figure}[H]
1885| \centering
1886| \includegraphics[width=0.58\textwidth]{DH-visc-dist.png}
1887| \caption{Empirical distributions of
1888| $\ln\mu_{DH}(g_n^{DH})$ for the $11\,525$ good DH Gram points (blue) and
1889| $751$ isolated bad DH Gram points (red) in the data set
1890| $2\leq n\leq12501$.  The vertical line marks $\mu=4$.  In sharp contrast
1891| with Fig.~\ref{fig:viscosity_distributions}, the isolated-bad DH distribution
```
Checked against the Zenodo data: 220,373 good and 24,269 isolated bad among n ≤ 250,000 (l. 342–343; the Zenodo README says GBG = 24,270); smallest isolated-bad μ ≈ 4.26178 (l. 353); 11,525 good and 751 isolated bad DH Gram points for 2 ≤ n ≤ 12501 (l. 1888–1889).

## 2. Qualitative statements (not tested here)

Two-lobe split of the isolated-bad prime cloud (l. 296–304); correlation of lobe position with μ (l. 370–375); Spearman correlation decreasing with multiplicative depth (l. 965–991); lowest-viscosity points at the lobe tips with one zero pushed toward a wall (l. 1735–1745); DH: no lower boundary at μ = 4 (l. 1891–1902); DH isolated-bad prime cloud is one component (l. 1924–1931).
