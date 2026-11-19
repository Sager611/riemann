
# Thm 5.2 audit (follow-up, 19 Nov 2026) — working notes

## A. Statement and source
- BTY = Bellotti, Trudgian, Yang, "Zero-free regions inspired by work of Heath-Brown", arXiv 2603.21490 (v1, 23 Mar 2026;
  e-print Bellotti_Trudgian_Yang_v1.tex), main.tex l.1182.
- Padhi §4.1 (l.376-417): w (BTY Def 1, θ = 1.1338, supp [0, D], D = 2θ cot θ), W(z) = ∫_0^D e^{-zu} w(u) du, w0 = 5.672787598;
  κ_m (BTY (18)), κ = 433/859; P = Σ_{k<=16} a_k cos(kx) >= 0 (BTY §4), a0 = 1, a1 = 865534/497079, a = 2919857/828465; K = 16, T0 = 1e10.
  f(u) = η w(ηu) Σ_m κ_m e^{-(2σ-1)mu}, F(z) = Σ_m κ_m W((z + (2σ-1)m)/η) (= Laplace transform of f: checked).
- eq:bty (l.393): for a zero ρ0 = β0 + it, t >= H0, η = 1-β0, λ0 = η log t, zero-free hypothesis β < 1 - A/log γ on
  [H0, Kt+T0], σ = 1 - A/log(Kt+T0), η <= η0 <= 0.0145, σ >= σ0 >= 0.993:
  c10 η w0 <= η w0((aκ/2) log t - G) + a0 F(σ-1) - a1(F(σ-1+η) + F(σ-η)) - Z + 1e-10 η + 2e-7 + E2 λ0 η^2 + E3 η^3, Z >= 0.
- M1 (l.400) = a1(F(σ-1+η)+F(σ-η)) - [η w0((aκ/2) log t - G - c10) + a0 F(σ-1) + 1e-10 η + 2e-7 + E2 λ0 η^2 + E3 η^3]; M1 > 0 => no zero.
- Derivation (l.405-417): BTY Lemma 6 identity Σ_n Λ(n) n^{-σ} f(log n) P(t log n) = Σ_k a_k Re(F(s_k-1) - Σ_ρ F(s_k-ρ) + Γ_k + E_k),
  s_k = σ + ikt; LHS >= c10 η w0 (Padhi L10-11); pole k>=1 <= 1e-10 η (BTY L13); k = 1 keeps ρ0, 1-ρ̄0; other pairs ρ, 1-ρ̄
  contribute <= 0 to -ΣF by L5 (Z = minus their total); Γ_k by L12, E_k by L8-9, tails <= 2e-7.
- Lemma twozero (l.558): ρ1 = β1 + i(t+d), β1 > 1/2, σ - β1 > 0, η1 = 1-β1: Z >= a1 Re(F(σ-1+η1-id) + F(σ-η1-id)).
  Check: s1-ρ1 = σ-β1-id = (σ-1+η1)-id; s1-(1-ρ̄1) = σ-1+β1-id = σ-η1-id. OK. Pair ≠ ρ0's pair (β1 > 1/2 > 1-β0). OK.
  Other k: ρ1's pair still >= 0 by L5 (same real parts). OK.
- M2 = M1 + a1 Re(F(σ-1+η1-id) + F(σ-η1-id)) > 0 for one ordering => the pair cannot exist (l.577).
- Thm 5.2 (l.580): 29 <= l_k <= 37 integer, λ*(l_k) < k·1e-4·l_k <= 0.0125(l_k - 1); two zeros with H0 < γ_j <= e^{l_k},
  1-β_j <= k·1e-4 => |γ0-γ1| >= d*(l_k, k). Proof (l.585): e_min = λ*(l_k)/l_k, e_max = k·1e-4; boxes
  (log t, η0, η1, |d|) ∈ [l-1, l+1e-12] × [e_min, e_max]^2 × [0, d*] for each l <= l_k; coarse k grid, next larger k for the rest;
  final d*(l_k, k) = min_{l <= l_k} d*(l, k).

## B. Inputs from BTY vs their hypotheses
- BTY Lemma 6 (F_sum_lem, TeX l.507): ζ(β0+it) = 0, β0 > 1/2, t >= H, A ∈ (0, A0], η ∈ (0, A0/log t], no zeros with
  β >= 1 - A/log γ, γ >= H; s_k = σ + ikt, σ = 1 - A/log(Kt+T0) => Re Σ_ρ F(s_1-ρ) >= F(σ-β0) + F(σ-1+β0) - 1e-8,
  Re Σ_ρ F(s_k-ρ) >= -1e-8 (k ≠ 1). A0 = 1/4.8596 (TeX l.404), η0 = A0/log H = 0.0071093, σ0 = 1 - A0/log(KH+T0).
- Proof (l.531-590): |γ| <= T = Kt+T0: zeros grouped {ρ, 1-ρ̄}; each pair Re F(z) + Re F(2σ-1-z̄) >= 0 (Lemma 5, needs
  0 <= Re z <= 2σ-1, from the zero-free hypothesis); k = 1 isolates ρ0's pair. |γ| > T: Lemma 5 (ii) termwise,
  > -44η(1e-10/(2π) + 6e-11)(log t + 1) > -1e-8 using η <= A0/log t.
- Padhi's range exceeds BTY's hypotheses (η up to 0.0145 > η0 = 0.0071; λ0 = η log t up to 0.35 > A0 = 0.206; σ down to 0.993
  < σ0 = 0.99347). He re-proves what depends on it: Lemma L5 (l.442: pair non-negativity for σ >= 0.993, η0 <= 0.0145 via
  B(y) > 0 with ε0 = C_w η0^2/(w0(2σ-1)^2) <= 1/500; check: 51.9·0.0145^2/(5.67279·0.986^2) = 0.0019786 < 0.002 OK);
  the tails (l.500): η(log t+1) <= 0.015·46 = 0.69 -> 44·0.69·7.6e-11 = 2.3e-9 per k <= 3e-9, 17 k's + k=0 term
  (88·0.015·log T0/T0 = 3.0e-9) <= 6e-8 < 2e-7 OK. Lemma 5 (ii) constant 44 is scale-free in η (W(ζ) = w0/ζ + O(|ζ|^-3)); not re-derived.
- E2, E3 (L89, l.459), c10 (L1011, l.473), G (L12, l.492): re-derived by Padhi for (η0, σ0); computed in his code by mpmath
  (bty2.py), certified only via the outward factor 1 ± 1e-6 (rep_cert.consts): directions correct (E2, E3 up; c10, G down).
  Paper l.416 says "evaluated in ball arithmetic": for these four constants that is not what the code does (minor).
- RamRepair (23 Oct) re-derived the BTY-type iteration independently: λ(29) = 0.21494 vs Padhi 0.21384 (his is conservative).

## C. Lemma twozero and Thm 5.2, line by line
- l.558-571 (Lemma twozero): valid. Keeping a second pair in BTY's k = 1 sum only replaces a dropped non-negative term by its
  value; the pairs {ρ0, 1-ρ̄0}, {ρ1, 1-ρ̄1} are distinct (β1 > 1/2 > 1-β0); all other pairs >= 0 by L5 because
  0 <= σ-β <= 2σ-1 for zeros below Kt+T0 (zero-free hypothesis + symmetry). Multiplicity: a double zero counts twice, still valid.
- l.573-577 (M2): M2 = M1 + a1 Re(F(σ-1+η1-id) + F(σ-η1-id)), F with the primary's η: matches rep_cert.margin (r1 + r2 - lhs).
- l.585 (proof): e_min = λ*(l_k)/l_k (code: lam_at(lk)/lk, value at the next grid point >= lk: <= λ*(lk), safe);
  e_min > A/log(Kt+T0) needed for σ-β1 > 0: e.g. 0.00737 > 0.0067 (l_k = 29).
- Boxes [l-1, l+1e-12] × [e_min, e_max]^2 × [0, d*] (code: Lt = ball(l0, l1+1e-12), η's in [emin, de], d in [0, dcert]).
  Only |d| matters (F(z̄) = conj F(z)). Label t = the lower zero: the upper one is at t+d, log(t+d) <= l1 + 1e-13.
- σ = 1 - A/log(Kt+T0), A = lam_at(log(K e^{l1} + T0)) (λ* non-increasing in l: checked on all 61 grid values, so A at the
  box top is the minimum over [H0, Kt+T0]: valid). σ-ball spans t ∈ [e^{l0}, e^{l1}]; in the swapped ordering the primary is
  the upper zero, whose σ can exceed the ball by ~1e-15 when the pair straddles an interior box edge (formal nit; margin effect
  ~1e-13, far below the certified margins). Constants at η0 = ceil(1e4 e_max)/1e4, σ0 = floor(1e5 σ_lo)/1e5: safe directions.
- λ*/l non-increasing (claimed l.538): checked on the 61-point grid (lamstar2c.json): true; λ* itself non-increasing: true.

## D. Recomputed certifications (his rep_cert.certify, ball arithmetic, maxdepth 40, his d* values; code/repcert_cells.py)
- (29, 76): d* = 0.061013 (directly certified cell): ok = True, 2979 boxes, 29 s (logs/repcert_29_76.log).

## E. Independent checks (code/m2corner.py; mpmath quadrature of BTY's w, constants parsed from the BTY TeX)
- BTY TeX: c_k (l.757-780) = his cc; Σ c_k^2 = 14912370 -> a0 = 1, a1 = 865534/497079, a = 2919857/828465 (diff 7e-21);
  κ_m (18) (l.242) = his km, Σ = 433/859; w0 = sec^2θ(θ tanθ + 3θ cotθ - 3) = 5.672787599 = his W0; w (Def 1, l.203) = his w.
- d* table (dstar2c_all.json, 507 entries) <= min over l <= l_k of the certified value inherited from the next certified k
  above (monotone in e_max), and each l's first certified k = floor(1e4 λ*(l)/l) + 1: no entry exceeds what is certified.
- Constants at cell (29, 76) (η0 = 0.0076, σ0 = floor(σ_lo)): E2 = 904.0, E3 = 7994.6, c10 = 0.05980, G = 1.56862
  (paper: G(0.993) = 1.56839, G(0.9935164) = 1.56888; E2 λ0 = 904·0.2058 = 186 ≈ his 185.4 at BTY's point; c10 = 0.0598 as l.843).

## F. Corner margins and falsification (logs/m2corner.log; independent float M2 vs his ball margin rep_cert.margin)
- W by 400-point Gauss-Legendre on [0, D] (checked vs mpmath quad: 4e-14). At every corner below, independent M2 = his ball
  margin to all printed digits, so his code evaluates the paper's M2 (l.577) with a correct W.
- Corners log t = l_k, (η0, η1) ∈ {e_min, e_max}^2, d = d*: M2 > 0 everywhere; d_fail = first d with M2 = 0 at that corner:
  | cell (l_k, k) | d* | worst corner (η0, η1) | M2(d*) | d_fail/d* |
  | (29, 76) | 0.061013 | (e_max, e_min) | 4.08e-3 | 1.010 |
  | (31, 69) | 0.129986 | (e_max, e_min) | 1.25e-3 | 1.031 |
  | (30, 72) | 0.092150 | (e_max, e_min) | 4.78e-3 | 1.034 |
  | (32, 68) | 0.055368 | (e_max, e_min) | 5.52e-3 | 1.019 |
  | (29, 120) | 0.058174 | (e_max, e_min) | 3.41e-1 | 1.116 |
  Corners with η0 = e_min have M1 > 0 (primary alone excluded) and never fail. (37, 121) has d* = 0 in the table (set to 0).
- Mechanism (numbers): at η0 = e_max, M1 < 0 (e.g. -0.0346 at (29,76), -0.707 at (29,120)); the second pair's term
  a1 Re(F(σ-1+η1-id) + F(σ-η1-id)) (0.047, 1.13) is what makes M2 > 0; it falls with d, hence a finite d*.
- e_min > 1 - σ (Lemma twozero's σ - β1 > 0) at the box bottom: l_k = 29: 0.007374 > 0.00676; l_k = 37: 0.005707 > 0.005425.
- Reruns: (29,76) ok 2979 boxes 29 s; (32,68) ok 5927 boxes 65 s; (31,71) ok 7833 boxes 110 s; (30,72) killed at 420 s;
  (31,69) killed after ~10 min (both large d*, not reproduced; certified in his run).
## G. Interior grid (logs/m2grid.log): 5×7×7×11 points per box, both orderings: min of max(M2, M2 swapped) at the top corner
  (log t = l_k, η0 = η1 = e_max, d = d*): (29,74) 1.069e-3, (29,76) 1.247e-2, (30,72) 5.605e-3, (31,69) 1.370e-3, (31,71) 1.496e-2,
  (32,68) 1.161e-2. Symmetric-corner d_fail/d*: 1.031, 1.036, 1.041, 1.035, 1.037, 1.048; (29,120) 1.166; (37,120) 1.262.
  The mixed corner (e_max, e_min) fails first in one ordering (1.010 d* at (29,76)) but the swapped ordering has M1 > 0.
## Verdict: THM 5.2 HOLDS (REPORT.md, section 'Thm 5.2 audit (follow-up)').
