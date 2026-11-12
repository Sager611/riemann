BARRIER: ζ's windowed Weil form is a killed Lévy jump process with prime jumps (its Markov part, needs Λ ≥ 0) plus anti-Markov terms (the pole on even functions, reflected prime atoms on odd ones), so Beurling–Deny tools (Barta/Allegretto–Piepenbrink, Cheeger, Muckenhoupt-type Hardy bounds, Holley–Stroock) certify positivity where a parity sector is Markov (even c ≤ ρ² = 1.7549, odd c ≤ 2) and, after flipping the anti-Markov terms, only to c ≈ 1.923 (even) / 2.06 (odd), below Yoshida's c = 2; beyond, an exact ground-state identity turns an even-sector certificate into Barta's pointwise condition plus a Poincaré inequality for the Doob-transformed process against the rank-one pole, whose relative margin even at the true ground state is 0.13, 1.2·10⁻⁴, 1.4·10⁻⁸ at c = 2, 3, 4 and ≤ 3·10⁻⁷³ at c = 20 (Galerkin; about the second even Weil eigenvalue), so a Poincaré estimate must be that sharp, and Davenport–Heilbronn is certified further (c ≈ 3.73) (node N1).

# explore39: functional inequalities → N1 (Orca exploration, 12 Nov 2026, 08:43–09:42 UTC)

Derivation: [DERIVATION.md](DERIVATION.md). Code: [code/](code/). Logs: [logs/](logs/). Conventions of ../WeilZeta and ../ExploreEntropy:
I = [−a, a], a = ½ ln c (CCM: λ² = c), Q = pole + Arch − primes, w_n = Λ(n)/√n, κ(z) = e^{z/2}/(e^z − e^{−z}), J(z) = κ(|z|) − 2cosh(z/2),
D = ψ(1/4) − ln π = −5.3722, ρ = 1.324718 (ρ³ = ρ + 1). All numerics float64 or mpmath Galerkin; nothing here is interval-certified.

## 1. Transfer (bridge measured before choosing)
| node | tool | reading | verdict |
|---|---|---|---|
| N1 | Cheeger/Buser, Barta, Poincaré, Muckenhoupt | Q = E_c − S(c) + pole, E_c the Dirichlet form of a killed Lévy process with Lévy measure κ(\|z\|)dz + Σ_{n<c} w_n δ_{±ln n} (explore17 T2: Markov because Λ ≥ 0; pole anti-Markov) | chosen |
| N5 | Muckenhoupt (Hardy operator, Mellin version) | the Hardy operator H*g = ∫_x^∞ g/u is the Mellin multiplier 1/s that produces the NB weight \|s\|⁻²; d_N is a distance to integer dilations of a non-L² measure, not an operator norm; the zero-kernel side is explore38's | not chosen |
| N13 | Perron–Frobenius, Poincaré gap | step one is RH-neutral (explore32) and certified to L ≤ 2 (Botha) | side remark only (§2 A, §3) |

Tools transferred: (b) Barta's inequality / Allegretto–Piepenbrink positivity principle (λ₀ ≥ inf Hφ/φ for Markov forms, exact at the ground
state), Cheeger's λ₀ ≥ h²/(2d) and its converse λ₀ ≤ h; (a) Muckenhoupt-type two-sided Hardy criteria, whose home is the nearest-neighbour
(birth–death) minorant; (c) Poincaré inequality with explicit constant (and Holley–Stroock's e^{±osc} transfer) via the ground-state transform.

## 2. Statements (proofs in DERIVATION.md)
A (sign structure, exact). Reduced to [0, a], the even sector has kernel J(x−y) + J(x+y) + atoms, the odd sector J(x−y) − J(x+y) + atoms.
 - Even: min of the continuous kernel is 2J(a), so it is Markov ⟺ c ≤ ρ² = 1.754878 (explore32's ρ is the full-space threshold); even
   prime atoms are Markov (w_n ≥ 0).
 - Odd: continuous kernel [κ(\|x−y\|) − κ(x+y)] + 4 sinh(x/2) sinh(y/2) ≥ 0 for all c (the pole is Markov-compatible here); the reflected
   prime atoms −w_n δ(x + y − ln n) are anti-Markov: Markov ⟺ c ≤ 2.
 - DH (no pole, b(3), b(4) < 0): even Markov ⟺ c ≤ 3, odd ⟺ c ≤ 2.
 - Hence Perron–Frobenius per sector: simple positive even ground state for c ≤ ρ², odd for c ≤ 2.
B (flip minorant). Q(f) ≥ Q_flip(\|f\|), Q_flip = Q − 2∬ j⁻ f f Markov; so λ₀(Q) ≥ λ₀(Q_flip) ≥ inf H_flip φ/φ. This is the canonical
  term-wise Beurling–Deny comparison of the explicit-formula terms.
C (even-sector certificate = Barta condition + Poincaré inequality; exact identity). For even ψ > 0, f = ψg:
  Q(f) − ∫(Qψ/ψ) f² = E_ψ(g) − 2⟨cosh(x/2), ψ⟩² Var_ν(g),   E_ψ(g) = ½∬ j_M ψ(x)ψ(y)(g(x) − g(y))²,   ν ∝ cosh(x/2) ψ dx,
  j_M = κ(\|x − y\|) + Σ_{n<c} w_n δ(x − y ∓ ln n). So Q^ev ≥ λ ⟸ (i) Qψ ≥ λψ pointwise and (ii) Poincaré P(E_ψ, ν) ≥ R := 2⟨C, ψ⟩².
  At the ground state ψ₀, (ii) is automatic with P − R = inf [Q(f) − λ₀‖f‖²]/Var_ν(f/ψ₀) ≤ (λ₁ − λ₀)/W(φ₁) and positivity is (i), λ₀ ≥ 0; for a
  test function ψ, (i) needs inf Qψ/ψ ∈ [0, λ₀] and (ii) must be proved with relative accuracy (P − R)/R. The arithmetic enters as prime-jump
  rates w_n ψ(x)ψ(x ± ln n).
C′ (odd-sector certificate beyond c = 2). For a test function φ > 0 on (0, a], the generalized Barta certificate needs Hφ ≥ 0 pointwise and
  E⁺_φ ≥ E⁻_φ, where E⁻_φ collects the reflected prime jumps w_n φ(x)φ(ln n − x) on x + y = ln n (negative rates) and E⁺_φ the Markov ones; at
  the ground state the domination is automatic and positivity is λ₀ ≥ 0 (DERIVATION §5′).
D (pinning, interlacing). Q^ev = M^ev + 2CC^T gives λ₁(E_c^ev) − S(c) ∈ [λ₀(Q^ev), λ₁(Q^ev)], S(c) = 2Σ_{n<c} Λ(n)/√n − ψ(1/4) + ln π,
  and λ₀(E_c^ev) = S(c) + μ₀ < S(c) wherever μ₀ < 0 (numerically c ≥ 1.24, μ₀ ≤ −1.28 at every window used; not for c ↓ 1, where exterior
  killing makes μ₀ > 0; DERIVATION §6): the killed Lévy process
  has a cluster of decay rates exponentially close to S(c) above an isolated lowest rate. Equivalently, the h-transform by E's Perron vector is
  a conservative Markov process whose spectral gap λ₁(E) − λ₀(E) = |μ₀| + μ₁ is pinned at |μ₀| to relative accuracy μ₁/|μ₀| ≤ λ₁(Q^ev)/|μ₀|
  (1.3e−12, 4e−36, 2.6e−73 at c = 5, 10, 20, Galerkin): a Cheeger/Buser window [h²/(2d), 2h] for this gap would have to be that narrow.

## 3. Tests
Smallest cases (Markov ranges; logs/t2.log; Barta-LP = best Barta bound over the Legendre span (n = 20) on a 400-point grid; Galerkin = upper bound):
| c | sector | Barta-LP | Galerkin |
|---|---|---|---|
| 1.2 | even / odd | 0.519901 / 1.423094 | 0.520066 / 1.423348 |
| ρ = 1.3247 | even | 0.251467 | 0.251597 |
| 1.5 | even / odd | 0.090195 / 0.620909 | 0.090282 / 0.621165 |
| ρ² = 1.7549 | even | 0.014178 | 0.014212 |
| 2 | odd | 0.072845 | 0.073109 |
Two-sided within 0.4 %: in the Markov ranges Barta is an explicit, checkable lower bound (one positive test function).
Same window, other tools (c = 1.2, full form, N cells; logs/t6.log): the killing is positive (min 0.2994), Cheeger ≤ h²/(2d) = 0.044, 0.033,
0.026 at N = 50, 200, 800 (d = ln N + 0.29: the log-singular Lévy measure has infinite intensity, so Cheeger is vacuous in the continuum);
the nearest-neighbour minorant, where Muckenhoupt's two-sided criterion applies, gives 0.323, 0.305, 0.301 → inf W = 0.2994 (the killing
bound), and inf W = 0 at c = 1.356980. True λ₀ = 0.52007.

Flip minorant at scale (logs/t3_zeta.log, t4.log; λ₀(Q) Galerkin, n = 20–28; for c ≥ 5 float64 cannot resolve λ₀(Q), true values 9.2e−18
(c = 5, ChekAudit), 1.25e−43 (10, LostLaw certified), 1.9e−242 (47, Chekmachev)):
| c | λ₀(Q^ev) | λ₀(Q^ev_flip) | λ₀(Q^odd) | λ₀(Q^odd_flip) |
|---|---|---|---|---|
| 1.9 | 3.36e−3 | +1.50e−3 (Barta 1.45e−3) | 0.153 | (Markov) |
| 2.0 | 1.33e−3 | −4.39e−3 | 0.0731 | (Markov) |
| 2.05 | | | 0.0502 | +0.0143 (Barta 7.3e−3) |
| 2.1 | 5.73e−4 | −0.0151 | 0.0354 | −0.0602 |
| 3 | 5.6e−8 | −0.356 | 1.5e−5 | −0.900 |
| 5 | 9.2e−18 | −1.552 | | −1.979 |
| 10 | 1.25e−43 | −3.917 | | −3.986 |
| 20 | ≈1e−96 (fit) | −7.180 | | −6.695 |
| 47 | 1.9e−242 | −12.855 | | −11.517 |
Flip-certified ranges (Barta crossing / Galerkin crossing): even c ∈ [1.92281, 1.92375], odd c ∈ [2.05483, 2.06027]. The anti-Markov mass that
must be flipped grows (‖F‖ = 1.7, 4.0, 7.2, 12.9 at c = 5, 10, 20, 47; even sector) while the margin falls like e^{−4πc}.

Poincaré margin of statement C (logs/t8.log float64 Legendre; logs/t5b.log mpmath CCM basis, upper bound via φ₁, edge zone below 1e−10·max
dropped from W, which keeps it an upper bound):
| c | 1.2 | ρ² | 2 | 2.5 | 3 | 3.5 | 4 | 5 | 10 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|
| R = 2⟨C,ψ₀⟩² | 0.348 | 0.986 | 1.108 | 1.237 | 1.304 | 1.347 | 1.379 | 1.419 | 1.497 | 1.544 |
| (P − R)/R | 3.27 | 0.429 | 0.134 | 6.7e−3 | 1.2e−4 | 1.6e−6 | 1.4e−8 | ≤ 6.6e−13 | ≤ 3.4e−36 | ≤ 3.2e−73 |
P − R tracks (λ₁ − λ₀)/W(φ₁) (1.5392e−4 vs 1.5587e−4 at c = 3): the margin is the second even Weil eigenvalue, not λ₀.
Odd analogue C′ (logs/t11.log): relative margin of the domination E⁺_ψ ≥ E⁻_ψ at the odd ground state, η = 2.65, 1.82, 1.06, 0.29, 1.5e−2,
5.0e−4, 8.0e−6 at c = 2.05, 2.1, 2.2, 2.5, 3, 3.5, 4 (≈ odd gap/5); looser than the even sector at equal c, as explore32's 'even fails first' predicts.

Pinning (logs/t5.log, wmat CCM basis at 4πc/ln 2 + 128 bits; interlacing checked in every Galerkin space):
| c | S(c) | λ₀(E^ev) = S + μ₀ | λ₁(E^ev) − S (Galerkin) | interval [λ₀, λ₁](Q^ev) (Galerkin) |
|---|---|---|---|---|
| 5 | 8.314157 | 4.980302 | 4.45e−12 | [1.36e−17, 6.96e−12] |
| 10 | 12.447189 | 7.328454 | 2.10e−35 | [5.73e−42, 3.27e−35] |
| 20 | 18.387837 | 11.172048 | 1.89e−72 | [1.75e−79, 2.94e−72] |
The Galerkin λ₀ are not converged (truth 9.2e−18, 1.3e−43, ≈1e−96; at c = 47 N = 129 gives 7.1e−176 against the certified 1.9e−242), but
λ₁(Galerkin) is a valid upper end and the certified ε_c a valid lower end.

Parity (logs/t7.log): the reflection identity Q(Eg) − Q(Rg) = −4∬ j(x+y) g g with g the positive odd ground state gives λ_ev ≤ λ_odd − 2I(c),
I = 0.285, 0.210, 0.113 at c = 1.2, ρ, 1.5, but I < 0 from c = ρ² on (−0.011; −0.118 at c = 2). So CCM's step one follows from Markov tools
for λ² ≤ ρ² only with a numerical comparison (Barta(odd) > Rayleigh(even)); analytically explore32's ρ stands.

## 4. Adversarial
- Davenport–Heilbronn (logs/t3_dh.log, t4.log): no pole, Markov for c ≤ 3 (even) and c ≤ 2 (odd); flip-certified even c ≤ 3.734–3.765, odd
  c ≤ 5.620–6.047, so the whole DH form is certified to c ≈ 3.73 against ζ's 1.92. At DH's failing window c = 31 (λ_min = −6.1e−30) the flip
  minorants are −5.26 (even), −4.91 (odd), as uninformative as ζ's −7.2 at c = 20. Markov certifiability tracks conductor and gamma factor
  (DH's margins are larger at equal c), not RH.
- Moved ζ pairs (explore32's 14 models, first failure at c ∈ [2.03, 7.21]): their forms are Q_ζ + Δ with Δ a smooth finite-rank kernel; the
  Markov ranges end at c = 1.92 (even), and at c ∈ [2.03, 7.21] ζ's own Poincaré margin is ≤ 0.13 → ≤ 10⁻¹², so no constant-slack bound
  separates ζ from a model that fails there.
- The Markov structure itself is Euler-product information (Λ ≥ 0; explore17), absent for DH and complex characters, but it does not feed
  positivity: beyond c ≈ 2 positivity is carried entirely by the anti-Markov terms (pole, reflected prime atoms). No RH information.

## 5. Verdict and novelty
- Barrier, in two strengths. (i) Exact for the canonical comparisons: the term-wise flip certifies nothing beyond c ≈ 1.923 (even), Cheeger's
  bound is 0 (infinite jump intensity), the nearest-neighbour minorant tends to the killing bound (c ≤ 1.357). (ii) Worst-case for every other
  estimate (reviewer: a warning, not an impossibility theorem): certifying P ≥ R from a lower bound L ≤ P needs L/P ≥ 1/(1 + ε), ε = (P − R)/R ≈
  λ₁(Q^ev)/(R·W), so a method whose guaranteed slack is a constant (Muckenhoupt 4, Cheeger's square, Holley–Stroock e^{osc}) gives no
  guarantee once ε is below it, i.e. from c ≈ 2 on. This is not impossibility (reviewer, round 2): a sharp analytic bound could exist, but it
  would have to resolve the second even Weil eigenvalue, and in addition (i) needs a test function with inf Qψ/ψ ∈ [0, λ₀]. Sampled windows
  only; the Galerkin margins are not validated bounds.
  Holley–Stroock in particular: transferring a Poincaré constant from the flat weight to ψ₀ costs e^{osc ln ψ₀²} ≥ (max ψ₀/min ψ₀)² ≈ 10¹⁶ at c = 5
  (ψ₀ min/max = 9.9e−9, logs/t5b.log) against a needed accuracy 6.6e−13.
- New as far as found: the per-sector Markov thresholds (even ρ², odd 2; DH even 3), Barta certificates in those ranges, the flip ranges
  1.923/2.055 (DH 3.73/5.62), the Barta + Poincaré certificate form of even positivity with its margin table, the pinning of the killed Lévy process's decay
  rates at S(c). All routine consequences of explore17's Markov decomposition; none moves a node quantity (analytic windows stay ≈ log 2).
- Prior work: Suzuki, arXiv 2606.09096 v3 (TeX read: Thm 4 proves a positive simple even ground state for sufficiently small a via the
  Markov property and irreducibility of a limiting Dirichlet form, citing Lenz–Stollmann–Veselić's Allegretto–Piepenbrink theorem; no
  explicit range, no flip, no Poincaré margin); explore17 T2 (Markov part, pole anti-Markov), explore32 S2 (full-space PF exactly for c ≤ ρ).
  arXiv search (5 queries, one rejected) found nothing else; rg of our NOVELTY files and reports: no Barta/Allegretto/Doob usage.
- Reviewer: see §6.

## Suggested NOGO line
- explore39 (Orca, exploration: functional inequalities — Barta/Allegretto–Piepenbrink, Cheeger/Buser, Muckenhoupt-type Hardy criteria, Poincaré/Holley–Stroock → N1, N13; N5 not chosen: no operator-norm structure, explore38 covers the kernel side) BARRIER: Q = E_c − S(c) + pole with E_c a killed Lévy process with prime jumps (Λ ≥ 0); Beurling–Deny tools certify positivity where a parity sector is Markov (even c ≤ ρ² = 1.7549, odd c ≤ 2; Barta brackets λ₀ within 0.4 %) and, flipping the anti-Markov terms (pole on even, reflected prime atoms on odd functions), to c ≈ 1.923 (even) / 2.06 (odd), below Yoshida's 2; Cheeger is vacuous (infinite jump intensity), the nearest-neighbour (Muckenhoupt) minorant stops at c = 1.357; beyond, an exact ground-state identity makes an even-sector certificate Barta's pointwise condition (inf Qψ/ψ ∈ [0, λ₀]) plus a Poincaré inequality for the Doob-transformed process against the rank-one pole, whose relative margin even at the true ground state, 0.13, 1.2e−4, 1.4e−8 at c = 2, 3, 4 and ≤ 3e−73 at c = 20 (Galerkin, ≈ second even Weil eigenvalue), is the sharpness a Poincaré estimate needs (worst-case, not impossibility; reviewer); the process's second decay rate is pinned within λ₁(Q^ev) of S(c) = 2Σ_{n<c}Λ(n)/√n − ψ(1/4) + ln π; DH is certified further (c ≈ 3.73): no RH information [ExploreFunctional](Riemann/Deep/ExploreFunctional/REPORT.md)

RH relevance: N1 (analytic window positivity: why functional inequalities stop at Yoshida's log 2), N13 (per-sector Perron–Frobenius ranges for CCM's step one).

## 6. Reviewer
Verdict-only review (claims inline; INCOMPLETE): (1) Markov decomposition CORRECT; (2) sector thresholds and atom signs CORRECT, with
ground-state positivity a.e. and an analytic convexity certificate for J (now in DERIVATION §0); (3) flip domination CORRECT on the folded
sector space, GAP: Barta needs admissible φ (stated in DERIVATION §3: polynomials positive on I, in the form domain); (4) Poincaré form CORRECT
(unit-norm φ₁, W > 0 added); (5) pinning bracket CORRECT, but μ₀ < 0 is INCORRECT as stated (fails as c ↓ 1): fixed, μ₀ < 0 proved for
c ≥ 1.24 (DERIVATION §6); (6) barrier logic a fair worst-case warning, not an impossibility theorem, and Galerkin margins are not validated:
wording changed (§5). Round 2 (revised claims): Barta admissibility CORRECT (sector-adapted φ, φ(0) = 0 odd); μ₀ < 0 for c ≥ 1.24 GAP
(needs a validated negative Rayleigh quotient: now stated as numerical); 'an accurate estimate is a numerical eigenvalue computation' INCORRECT
(sharp analytic bounds can exist): removed; C′ identity CORRECT but 'odd positivity ⟺ E⁺ ≥ E⁻' INCORRECT (automatic at the ground state;
positivity is λ₀ ≥ 0): C and C′ restated as certificates (pointwise Barta condition + Poincaré/domination), the margins as margins of the
second condition at the true ground state. Robustness (logs/t10.log): flip crossings 1.92376/1.92375/1.92375 (even) and 2.06032/2.06027/2.06024 (odd) at
n = 14/20/28; Barta-LP unchanged from 400 to 800 grid points (1.44774e−3 → 1.44763e−3 at c = 1.9 even; 7.331e−3 → 7.326e−3 at c = 2.05 odd).

## Files
code/ (run with uv run --no-project --with numpy --with scipy [--with mpmath --with python-flint] python3 <script>, OMP etc. = 1, each under
.scratch/memcap.sh 290):
- xform.py: x-space parity Legendre Galerkin (autocorrelation form of Arch), pointwise operator H, flip matrices, Hmat, barta_lp.
- t0_valid.py (validation: c = 2 even 1.3293e−3, odd 0.07311; explore32 1.330e−3, 0.075), t1_sign_barta.py (kernel signs, first Barta),
  t2_barta_lp.py (Barta-LP brackets), t3_flip.py zeta|dh, t4_bisect.py (flip ranges), t5_pinning.py, t5b_vectors.py (mpmath, CCM basis),
  t6_cheeger.py (cell chain, Cheeger, nearest-neighbour), t7_parity.py, t8_poincare.py, t9_consts.py, t10_robust.py (n and grid
  convergence), t11_odd_margin.py (C′).
- wmat.py copied unchanged from ../ExploreEntropy/code (explore19's CCM-basis form, ζ and DH).
logs/: one log per script; peak RSS 41–116 MB (/usr/bin/time -l), every process single-threaded, one at a time.
