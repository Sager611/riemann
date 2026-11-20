BARRIER: certificates of Weil positivity built from short-window (link) spectral data, as in Garland's method and Oppenheim's trickling-down, need link floors ≥ Λ*(c, c′), the Kronecker-aligned prime atoms no link holds (≥ 0.5616 past c = 47, ≥ 0.1169 past λ² = 37.923), while links have dozens of eigenvalues below 10⁻³ (nodes N1, N13).

# explore64: high-dimensional expanders → N1, N13 (Orca exploration, 20 Nov 2026, from 15:20 UTC)

Code: this directory (`hdx_*.py`). Logs: [logs/](logs/). Conventions of ../ExploreFunctional (xform.py, wmat.py): I = [−a, a], a = ½ ln c,
Q(f) = Arch(f) + pole(f) − 2Σ_{n<c} w_n C_f(ln n), w_n = Λ(n)/√n, C_f(h) = ∫ f(x) f(x−h) dx, Arch(f) = D‖f‖² + ½∬ κ(|x−y|)(f(x)−f(y))²,
κ(z) = e^{z/2}/(e^z − e^{−z}), pole(f) = 2⟨f, e^{x/2}⟩⟨f, e^{−x/2}⟩. So Q(f) = ∬ f(x) f(y) K(x−y) with K = Dδ − κ(|u|) + 2cosh(u/2) − Σ w_n(δ_{ln n} + δ_{−ln n}) (finite part at 0).
e_k(c) = k-th eigenvalue of Q on L²[−a, a] (both parities), λ² = c in CCM's notation.

## 1. Transfer
Tool: Garland's method / Oppenheim's trickling-down (and spectral-independence local-to-global): a global quadratic form is an average of
its restrictions to links plus a remainder; link spectral gaps (scalar floor, or floor above the link's trivial eigenvector) bound the
global spectrum. No group: the decomposition is of Weil's form itself, links = short windows (translates; K is translation invariant).
Nodes considered: N1 and N13 (chosen; one barrier covers both); N6 (trickling-down for Lorentzian polynomials gives log-concavity of
Jensen coefficients = Turán inequalities, known: Csordas–Norfolk–Varga 1986, Dimitrov–Lucas; weaker than hyperbolicity); N5, N8 (no
link structure in the quantity). R474's Zuk/property (T) route through a group action is not repeated.

Localization identity (exact; IMS form of Garland's averaging). Let ν be a measure on links W_b ⊂ I (each of diameter ≤ ln c′, c′ < c) and
χ_b cutoffs supported in W_b with ∫ |χ_b(x)|² ν(db) = 1 on I and (1 − Re M)κ(|x−y|) ∈ L¹(I²) ("admissible": Lipschitz cutoffs, translates of
hard windows, tilings by indicators). Put M(x, y) = ∫ χ_b(x) χ_b(y) ν(db)
(PSD, M(x, x) = 1, |M| ≤ 1, M = 0 for |x−y| > ln c′). Then for f ∈ L²(I) in the form domain
  Q_c(f) = ∫ Q(χ_b f) ν(db) + R_M(f),   R_M(f) = ∬ f(x) f(y) K(x−y)(1 − M(x, y)) dx dy,
where the δ-part of K is killed (1 − M = 0 on the diagonal) and −κ(1 − M) is locally integrable. Proof: expand Arch, pole and atoms
term by term (Arch: ∫ν (χ_b f(x) − χ_b f(y))² = (f(x) − f(y))² + 2(1 − M) f(x) f(y)).
Translation links (χ_b = χ(· − b), ν = db, ∫χ² = 1) give M = m(x − y), m = χ ⋆ χ̃; hard links χ = 1_{[−a′,a′]}/√(2a′): 1 − m(z) = min(|z|/2a′, 1).
Garland certificate (scalar): Q(χ_b f) ≥ e₁(c′)‖χ_b f‖² ⟹ Q_c ≥ e₁(c′) + λ_min(R_M). Trickling-down (rank r): Q(χ_b·) ≥ μ_b − F_b, rank F_b ≤ r,
so μ_b ≤ e_{r+1}(c′); r = 1 with the link ground state deflated is Oppenheim's use of the trivial eigenvector.

## 2. Statements (proofs here; numerics in §3)
**A (localization identity)** — §1, exact (real cutoffs; complex ones: Hermitian M, Re M throughout).

**B (Kronecker ceiling).** Fix c > c′ > 1 and an admissible (ν, χ) (§1: ∫ν χ_b² = 1, (1 − M)κ ∈ L¹(I²)) with links of diameter ≤ ln c′; link bounds B_b = μ_b − F_b ≤ Q on
L²(W_b), F_b compact (ν-ess. bounded); floor function μ(x) := ∫ μ_b χ_b(x)² ν(db), μ̄ := ess sup μ_b; certified operator
G := ∫ χ_b B_b χ_b ν(db) + R_M ≤ Q_c (proving Q_c ≥ 0, or e₂(Q_c) > 0, from G ≤ Q_c alone needs e₁(G) ≥ 0, resp. e₂(G) > 0). Put
  ⟨u, P_M u⟩ := Σ_{n<c} 2w_n ∫ u(x) u(x+ln n)(1 − M(x, x+ln n)) dx,  Λ_M := sup_{u ≥ 0 bounded} ⟨u, P_M u⟩/‖u‖²,  Λ*(c, c′) := inf_M Λ_M.
Then for every bounded u ≥ 0 and every k: e_k(G) ≤ [∫ μ u² − ⟨u, P_M u⟩]/‖u‖²; in particular e_k(G) ≤ μ̄ − Λ_M ≤ μ̄ − Λ*(c, c′).
Lower bounds for Λ* (any admissible M):
 (i) Σ_n 2w_n max(0, ⟨u, T_{ln n} u⟩/‖u‖² − cos(π/(⌈ln c′/ln n⌉ + 1))), any bounded u ≥ 0 (T_s = shift truncated to I);
 (ii) [1/Σ_i 1/(2 w_{q_i} α_{i−1} α_i)] / Σ_i α_i², any α_i > 0 and any rational r = Π q_i^{±1} ∈ (c′, c) (q_i prime powers < c) whose path of
     steps ±ln q_i has distinct positions p_i and span max p_i − min p_i < ln c ("integer/rational cycles").
Proof. f_T = u(x) cos(Tx), u bounded ≥ 0, T → ∞ along T ln p → 0 (mod 2π) for all p < c (recurrence of a linear flow on the torus). Then
‖f_T‖² → ½‖u‖²; link part → ½∫ μ u² (f_T ⇀ 0 and F_b compact, dominated by C‖χ_b u‖²); smooth part of R_M → 0 (Riemann–Lebesgue: the
kernel (1 − M)(2cosh(z/2) − κ) times u⊗u is in L¹(I²)); atom part → −½ Σ 2w_n cos(T ln n) ∫ u u_{ln n}(1 − M) → −½⟨u, P_M u⟩. k such times with
pairwise gaps → ∞ make all cross terms vanish: min-max on their span. (i): ∫ u(x)u(x+s)M(x,x+s) dx = ∫ ν(db) ⟨uχ_b, T_s(uχ_b)⟩ ≤ ω_s ‖u‖²,
ω_s = cos(π/(m+1)), m = ⌈ln c′/s⌉, the numerical radius of T_s on an interval of length ln c′ (a direct integral of path graphs with
≤ m vertices; longer chains only at endpoints, a null set); 1 − M ≥ 0 and u ≥ 0 allow max(0, ·). (ii): v_x := (χ_b(x))_b ∈ L²(ν) are unit
vectors with ⟨v_x, v_y⟩ = M(x, y); the end point is at distance ln r > ln c′, so no link holds both and v_{x₀} ⊥ v_{x_k}; the triangle
inequality gives Σ_i √(1 − M(x_{i−1}, x_i)) ≥ 1 for every start x₀ of a set I₀ of positive length (span < ln c); with u = α_i on I₀ + p_i
(disjoint after shrinking I₀) and dropping all other (nonnegative) terms, min Σ A_i ρ_i² subject to Σ ρ_i ≥ 1 is 1/Σ(1/A_i). ∎
Consequences. (a) A certificate needs ∫ μ u² ≥ ⟨u, P_M u⟩ for all bounded u ≥ 0, hence μ̄ ≥ Λ*(c, c′); on the packets of a cycle (ii) the
α²-weighted average of μ is ≥ the cycle bound. For N13 (gap e₂ > e₁ ≥ 0) a comparison certificate needs e₂(G) > (an upper bound of
e₁(Q_c)) ≥ 0, so again μ̄ > Λ*. (b) If rank F_b ≤ r then μ_b ≤ e_{r+1}(W_b) (min-max); for equal-length families (translates or tilings of
windows of length ln c′, the regular setting of Garland and Oppenheim) every link with μ_b > Λ* − ε (a set of positive ν-measure) has
rank ≥ N_{c′}(Λ* − ε), N_{c′}(θ) := #{k : e_k(c′) < θ}; Garland is r = 0, Oppenheim's trickling-down r = 1. For compact F_b of any rank,
F_b needs ≥ N_{W_b}(μ_b − δ) eigenvalues above δ, every δ > 0. (c) Mixed lengths: (a) holds verbatim but forces large floors only in
this weighted sense (a high-floor set met by every cycle's packets), not on every link.
Corollary (record windows). Λ*(c, 47) ≥ w₄₇ = ln 47/√47 = 0.56160 for every c > 47: pairs at distance exactly ln 47 lie in a common link
only on a null set (Fubini: one endpoint pair per link), and u = 1_{[−a, −a+δ]} + 1_{[−a+ln 47, −a+ln 47+δ]}, 0 < δ ≤ min(ln(c/47), ln 47), gives
2w₄₇⟨u, T_{ln 47} u⟩/‖u‖² = w₄₇.
Monotonicity (used in 3.5). Λ*(c, c′) is non-increasing in c′ (larger c′, larger admissible class) and non-decreasing in c (restrict a
localization of the larger window to a subwindow of length ln c₁; extra atoms add nonnegative terms); e_k(c′) is non-increasing in c′
(nested windows, min-max).

**C (smallness is inherited).** G ≤ Q_c, so the min-max values satisfy e_k(G) ≤ e_k(Q_c): at c = 47 the first 89 are < 10⁻³ (§3.4),
and e₁(G) ≤ e₁(Q_c). (This alone does not make G harder than Q_c — G = εI would be easy — but here G = floors − F + R_M keeps each atom
with weight w_n(1 − M), all of it when the links are atom-free.)

**D (short links are ζ-blind).** Links with c′ ≤ 2 contain no atom (an atom at distance exactly ln 2 is a null set): their spectra depend
only on the Γ-factor and the pole (the same for every Dirichlet series with ζ's Γ-factor and pole); all arithmetic sits in R_M.

## 3. Tests (numbers; logs in logs/)
**3.1 Smallest case: the natural Garland bound (float64, Legendre parity basis, 24 modes/sector; hdx_garland.py, logs/garland_n24.log).**
Translation links of half-width a′ = ½ ln c′, hard (1 − m = min(z/2a′, 1)) or cosine cutoffs. Garland: Q_c ≥ e₁(c′) + λ_min(R).
| c | c′ | e₁(c′), e₂(c′) even / odd | λ_min R hard (even/odd) | λ_min R cosine |
|---|---|---|---|---|
| 10 | 3.162 = c^½ | 9.3e−9, 4.1e−4 / 3.1e−6, 2.5e−2 | −1.493 / −1.512 | −1.580 / −1.597 |
| 10 | 5.623 = c^¾ | < 1e−13 (float floor), e₂ odd 3e−12 | −1.288 / −1.297 | −1.429 / −1.447 |
| 10 | 7.943 = c^0.9 | < 1e−13 | −1.123 / −1.125 | −1.334 / −1.349 |
| 20 | 4.47, 9.46, 14.8 | ≤ 1.3e−9 | −1.51, −1.49, −1.37 | −1.50, −1.52, −1.51 |
| 30 | 5.48, 12.8, 21.4 | ≤ 1.5e−11 | −1.64, −1.59, −1.40 | −1.63, −1.65, −1.63 |
| 47 | 6.86, 17.95, 31.98 | < 1e−13 | −1.39, −1.37, −1.26 | −1.52, −1.42, −1.41 |
Garland fails everywhere by ≥ 1.1 (margins ≤ 2.5e−2, mostly < 1e−9). Which part (hdx_decomp.py, logs/decomp_n24.log): λ_min of R_atoms
alone −1.1 … −6.3, R_pole + R_arch −0.3 … −5.0, and R restricted to modes ≥ 8 (high frequency) still −1.0 … −1.4: the negativity is the
prime atoms longer than (or partly outside) the links, not the pole. Low-mode Galerkin understates it (next items).
Nested windows, links almost the whole window (hdx_nested.py, logs/nested_n24.log): λ_min R (hard, even/odd) = −1.11/−1.12 for
(c, c′) = (40, 227/6), −1.14/−1.11 for (41.5, 227/6), −1.15/−1.20 for (49, 47), −1.18/−1.20 for (50, 47) (cosine links −1.28 … −1.44), against
link margins e₁, e₂ < 10⁻³: the translation average leaves every atom partly uncovered (weight ln n/ln c′), the smallest N13/N1 tests fail.
**3.2 Kronecker witnesses on actual functions (hdx_kron.py, logs/kron.log; f_T = (1 − x²/a²)² cos Tx, 40,001-point quadrature).**
c = 10, c′ = √10, T = 108.795: R(f_T)/‖f_T‖² atoms −0.7077 (formula of B: −0.7077), smooth part −0.0002. c = 47, c′ = √47, T = 280.833:
atoms −1.5969 (formula −1.5970), smooth −0.0000. The limit formula already holds to 4 digits at these T; the ceiling is its Kronecker limit with
the optimal envelope.
**3.3 Certified ceilings Λ*(c, c′) (arb lower bounds; hdx_ceiling.py N = 3000 cells, hdx_path.py, hdx_rat.py, hdx_beta0.py).**
| c | c′ | Λ* ≥ | source |
|---|---|---|---|
| 10 | √10, 10^¾, 10^0.9 | 1.6245, 1.0579, 0.5104 | atoms 4–9 long-range; 2, 3 via (i) |
| 20, 30 | √c | 2.9206, 3.6234 | long-range atoms |
| 47 | √47, 47^¾, 47^0.9 | 5.1538, 3.1687, 1.3715 | long-range atoms |
| 47 | 44 / 46.5 / 46.9 | 0.1988 (46 = 23·2) / 0.1084 (29·37/23) / 0.0560 (43·4·3/11) | cycles (ii) |
| c > 47 | 47 | 0.56160 = w₄₇ (exact); c = 50: 0.7144, 53: 0.7724, 60: 1.2744, 71: 2.0556, 100: 3.9659 | atom 47, then 49, 53, … |
| 48.5 | 47 | 0.1023 (48 = 16·3) | cycle |
| 72 / 80 / 100 | 71 | 0.5026 / 1.0738 / 2.1050 | atoms 71, 73, 79, … |
| 39 / 41.5 / 45 / 50 | 227/6 = 37.83 | 0.2017 (38 = 2·19) / 0.5800 (atom 41) / 0.9331 / 1.4339 | cycles, atoms |
| 37.95 | 227/6 | 0.1169 (29·17/13 = 37.923) | rational cycle |
Uniform: every integer 2 ≤ N′ ≤ 100 has a cycle worth ≥ β₀ = 0.0795 (weakest 96 = 3·32; 64: 0.0866; 90 = 2·5·9: 0.0884), so
Λ*(c, c′) ≥ 0.0795 whenever c ≤ 101 and (c′, c) contains an integer. No cycle with ≤ 4 steps found for (10, 9.9), (47, 46.99).
Growth (hdx_growth.py, logs/growth.log): Λ*(c, √c) ≥ 1.6247, 5.1540, 7.9171, 11.7505, 19.8152, 28.5391 at c = 10, 47, 100, 200, 500, 1000
(arb), and Λ*(c, c^θ) ≥ S(c) = Σ_{c^θ≤n<c} 2w_n(1 − ln n/ln c) (constant envelope) ~ 8√c/ln c by the PNT, any fixed θ < 1. Floors must grow like
√c/ln c at least, while a link's eigenvalues grow only logarithmically in their index (archimedean symbol log(t/2π)), so the needed link
rank grows roughly like exp(Λ*) (heuristic Weyl count, not certified): super-polynomial in c.
**3.4 Certified link spectra (arb Galerkin in CCM's Fourier basis with ../WeilSimple's audited WeilMat (rigorous tail), inertia of A_N − θ by arb LDLᵀ; hdx_count.py; N = 100: 201 modes, prec 400; logs/count_N100_ws.log, identical to the wmat run count_N100.log).**
Lower bounds on N_{c′}(θ) = #{k : e_k(c′) < θ}, both parities (still growing with N: 41 → 63 → 80 → 96 at c′ = 47, θ = 0.5616,
N = 30/60/100/180; N = 180 runs: logs/count_N180_ws.log, and logs/count_N180_p1400.log at prec 1400 where prec 400 left a pivot undecided):
| c′ | θ = 1e−6 | θ = 1e−3 | θ = Λ* of the extension |
|---|---|---|---|
| √10 | | ≥ 3 | ≥ 13 at 1.6245 |
| √47 | | ≥ 10 | ≥ 197 at 5.1538 (197 of 201 modes) |
| 227/6 | | ≥ 71 | ≥ 73 at 0.2017, ≥ 82 at 0.5755 |
| 47 | ≥ 86 | ≥ 89 | ≥ 91 at 0.1023, ≥ 96 at 0.5616 |
| 71 | | | ≥ 96 at 0.5026 (N = 100) |
So e₈₆(47) < 10⁻⁶, e₈₉(47) < 10⁻³, e₇₁(227/6) < 10⁻³ and e₃(c′) < 10⁻³ for every c′ ≥ √10: Oppenheim's rank-1 floor e₂ is ≪ Λ* at all such links.

**3.5 Corollaries (certified: B + 3.3 + 3.4 + hdx_shortcert.py, logs/shortcert_k1.log, shortcert_k2.log).** "Equal-length" = all links
of length ln c′ (translates or tilings of one window), the regular setting of Garland and Oppenheim; with all ranks ≤ r, μ̄ ≤ e_{r+1}(c′).
- Short links. Grid g = 1.05, 1.1, 1.2, 1.35, 1.5, 1.75, 2, 2.5, 3.162277 with e_k(g_j) < L_{j+1} ≤ Λ*(10, g_{j+1}) certified by inertia
  (L = 2.5098 for g ≤ 2, 2.0321 at 2.5, 1.6245 at √10). By the monotonicities, e_k(c′) < Λ*(c, c′) on every [g_j, g_{j+1}] and every c ≥ 10:
  equal-length Garland (k = 1) fails for every c′ ∈ [1.05, √10], Oppenheim (k = 2) for every c′ ∈ [1.1, √10].
- Long links. For c′ ≥ √10, e₃(c′) ≤ e₃(√10) < 10⁻³ < β₀ = 0.0795: equal-length certificates of rank ≤ 2 fail for every window c ≤ 101 and
  link c′ ∈ [√10, c) with an integer in (c′, c); the slivers without one need rational cycles (3.3: positive where tested).
- N1 record. Equal-length local-to-global extension of Q ≥ 0 beyond c = 47 (all f) from links of length ln 47 needs floors ≥ 0.5616, i.e.
  ≥ 96 resolved directions on every link with floor near μ̄ (N₄₇(0.5616) ≥ 96, 0.5616 < w₄₇; heuristic count of the sub-10⁻³ cluster ≈ 94:
  per parity, Fourier modes below 2πc minus N(2πc), 181 − 134).
- N13. Equal-length trickling-down of CCM's simple even ground state (e₂ > e₁) from λ² = 227/6 to any λ′² ≥ 493/13 = 37.923 needs floors
  ≥ 0.1169, i.e. ≥ 71 directions of QW_λ (e₇₁(227/6) < 10⁻³); beyond 38, floors ≥ 0.2017 (≥ 73 directions).
- Mixed lengths. B(a) still forces a high-floor set (weighted average over every cycle's packets ≥ the cycle bound), not a floor on every
  link; the rank counts above are then not forced link by link. Floors above Λ* come cheaply only from links shorter than 1.05, through
  the archimedean ln(1/a′) growth (float: e₁(1.01) = 3.14, e₁(1.005) = 3.82 vs Λ*(10, ·) = 2.51, Λ*(47, ·) = 7.30; logs/small.log); those hold
  no atom (D), so the arithmetic stays in R_M, and by C the first 89 min-max values of G stay < 10⁻³ at c = 47.
Freedoms covered: any links of diameter ≤ ln c′, any weights ν ≥ 0, any admissible cutoffs, any floors μ_b plus compact corrections, the
remainder kept exactly. Not covered: certificates that use more than G ≤ Q_c (keeping parts of Q − G: non-compact link bounds such as Arch
on links, cross-link couplings, preconditioners), and multi-window links (unions of windows at atom distances; their floors are
restrictions of Q, hence ≤ those of their component windows, so 3.4's counts apply to them).

## 4. Adversarial: Davenport–Heilbronn (hdx_dh.py, logs/dh.log)
Same localization of DH's Weil-type form (Γ-factor (5/π)^{s/2}Γ((s+1)/2), no pole, signed atoms b(n)/√n at all n). Long-range ceiling
sup_θ λ_max Σ_{c′≤n<c} b(n) cos(k(n)·θ)/√n (T_n + T_n*) (600 random Kronecker phases θ, arb Rayleigh quotient at rational θ):
| c, c′ | 10, √10 | 47, √47 | 47, 31.98 | 49, 47 | 60, 47 |
|---|---|---|---|---|---|
| ζ (θ = 0) | 1.6245 | 5.1538 | 1.3715 | 0.5613 | 1.2744 |
| DH | 1.7234 | 3.6018 | 1.9036 | 0.1661 | 1.4025 |
The barrier is method-intrinsic: it holds for DH with comparable constants, so it does not tell ζ from DH (and B holds for any form
whose atoms reach beyond the links). ζ-specific only in form: Λ(n) ≥ 0 makes θ = 0 extremal (all prime powers align at once, the
configuration the pole balances at T = 0) and lets the cycle bounds (ii) drop all other terms; DH's signed b(n) need phase-consistent signs.

## 5. Novelty (logs/novelty_arxiv.log, logs/novelty_listings.log)
arXiv API (10 phrase queries: Garland/zeta, trickling down/quadratic form, Weil positivity/local, localization/explicit formula, Kronecker/
Weil, semilocal Weil, domain decomposition/zeta, …), monthly listings 2026-10 and 2026-11 (math.NT, CO, SP, FA, CA, PR, cs.DM; 8 title hits,
none relevant: 2611.13129 and 2610.32704 seen, 2610.37721 is graph colouring), Zenodo (63 hits; nearest 24377825 single-window lowest
eigenvalue law, 24051794 conditional sampling bound, 23635192 screw-function viability): no localization of Weil's form over links, no
local-to-global certificate for it, no Kronecker or integer-cycle ceiling. Classics: Garland 1973, Oppenheim 2018, IMS localization
(Ismagilov, Morgan, Simon, Sigal), Lieb–Yau's nonlocal localization (1988), Bohr–Landau alignment for Ω-results: none concerns Weil's form.
Verdict: the ceiling seems unrecorded, but it is a barrier (elementary ingredients), not a discovery.

## 6. Verdict
BARRIER, unconditional (no RH input; numbers arb-certified except §3.1–3.2, float illustrations). Garland's method and Oppenheim's
trickling-down transfer to Weil's form exactly (identity A), but cannot move N1 or N13: by B every certificate built from link spectral
data (floor + compact, remainder exact) needs a floor function dominating the Kronecker-aligned operator P_M of the prime atoms the links
do not hold in full, so μ̄ ≥ Λ*(c, c′) > 0 as soon as (c′, c) holds a prime power or a rational with a short prime-power factorization
(integer and rational cycles): ≥ 0.0795 across any integer (c ≤ 101), ≥ 0.5616 for every extension of c = 47, ≥ 0.1169 past λ² = 37.923.
Links c′ ≥ √10 have ≥ 3 (c′ = √10) to ≥ 89 (c′ = 47) eigenvalues below 10⁻³, so equal-length scalar and rank-one link data fail
for every window c ≥ 10 and every link c′ ∈ [1.1, √10] (certified grid) or c′ ∈ [√10, c) with an integer in between; equal-length
finite-rank data need ≥ 96 link modes (N1, past c = 47) and ≥ 71 (N13, past λ² = 37.923). Links shorter than 1.05 pass only through the
archimedean log growth and hold no atom: every atom stays in R_M and the first 89 min-max values of G stay below 10⁻³ (c = 47).
TREE N1 record (c ≤ 47 all f; c ≤ 71 even) and N13 record (λ² ≤ 227/6) unchanged.
Scope (reviewer, logs/review1.txt): B binds certificates that use only G ≤ Q_c with link spectral data on links of bounded diameter;
it does not exclude methods that keep parts of Q − G (Arch on links, cross-link couplings, preconditioning) or multi-window links.
Not repeated: R474 (Zuk via a group action), explore39 (Markov/Cheeger tools), explore28 (structured elimination; a deflated version of
this localization is a two-level domain decomposition, an engineering variant of a global solve).

## Suggested NOGO line
- explore64 (Orca, exploration: high-dimensional expanders — Garland's localization and Oppenheim's trickling-down → N1, N13) BARRIER: Weil's form is exactly an average of window forms plus K(1 − M), but on Kronecker-aligned high-frequency tests every certificate from link spectral data (links of diameter ≤ ln c′, any cutoffs, floor + compact, using only G ≤ Q) has all e_k ≤ floor − Λ*(c, c′), Λ* = aligned prime atoms no link holds, plus integer/rational cycles (38 = 2·19, 29·17/13, …): Λ* ≥ 0.5616 past c = 47, ≥ 0.1169 past λ² = 37.923, ≥ 0.0795 across any integer (c ≤ 101); links have e₈₉(47), e₇₁(227/6), e₃(√10) < 10⁻³ (arb), so equal-length Garland/Oppenheim fail for all links c′ ∈ [1.1, c) with an integer in (c′, c) (c ≥ 10) and equal-length data need ≥ 96 (N1) / ≥ 71 (N13) link modes; shorter links are atom-free; DH gives comparable ceilings (method-intrinsic) [ExploreHDX](Riemann/Deep/ExploreHDX/REPORT.md)

RH relevance: N1 (window extension of Weil positivity) and N13 (CCM's simple-ground-state step on nested windows): a certified method
barrier; no change to either record.

## Files
hdx_garland.py, hdx_decomp.py (float Garland test), hdx_kron.py (Kronecker witnesses), hdx_ceiling.py, hdx_path.py, hdx_rat.py,
hdx_beta0.py, hdx_growth.py (arb ceilings), hdx_count.py (arb link spectra), hdx_shortcert.py (short-link grid), hdx_small.py,
hdx_nested.py (float), hdx_dh.py (DH). Reviews: logs/review1.txt (gaps, all addressed), logs/review2.txt (second pass, remaining items fixed).
Reused: ../ExploreFunctional/code/xform.py (float x-space form), wmat.py (arb Fourier-basis form, closed forms of ../WeilSimple).
Every run under .scratch/memcap.sh 290; sampled peaks 0–36 MB; one process at a time, single-threaded.
