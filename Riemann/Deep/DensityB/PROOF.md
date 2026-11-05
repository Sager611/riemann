# PROOF: zero-density exponents A(σ) below the Lean-checked half-plane

Node N7 (zero-free regions and zero density). Certificate: verify.py (this directory; exact rational arithmetic, standard
library only). Reproduction of the numerics: repro/ (our scripts; they download nothing and need a local copy of the
ANTEDB Python code, Apache-2.0, github.com/teorth/expdb, which is not included here). Paper proofs, not refereed, no Lean.

## 0. Setting and sources

B := 0.87338196661320752 = 87338196661320752/10^17.

(H) ζ(s) ≠ 0 for Re s > B. Lean-checked (QuasiRiemannTracker, J. Bebel, 31 Oct 2026; since 1 Nov the Lean-checked
half-plane is Re s > B_2 + 2.2·10^−254 with B_2 = (777 − √1081)/852 < B, so (H) holds a fortiori; verify.py checks
B > B_2).

Notation and lemma names follow Tao–Trudgian–Yang, arXiv 2501.16779v2 ("TTY"), which is the same text as the zero-density,
large-value and moment chapters of the ANTEDB blueprint (github.com/teorth/expdb, main branch read 2026-11-05):
- μ(σ): least exponent with ζ(σ + it) ≪ |t|^{μ(σ)+o(1)}; convex, non-increasing, μ(1 − σ) = μ(σ) + σ − 1/2.
- LV(σ, τ): large-value exponent of Dirichlet polynomials Σ_{n∈[N,2N]} a_n n^{−it} with |a_n| ≤ 1, values ≥ N^σ on a
  1-separated set W ⊂ [T, 2T], T = N^τ (TTY's normalisation); |W| ≪ N^{LV(σ,τ)+o(1)}.
- LV_ζ(σ, τ): the same for a_n = 1_I(n), I ⊂ [N, 2N] an interval ("zeta large values"), again with W ⊂ [T, 2T]
  (the restriction matters here: near t = 0 the sum is ≈ |I|; orchestrator fix after review, 5 Nov 2026).
- M(σ, A): ∫_T^{2T} |ζ(σ + it)|^A dt ≪ T^{M(σ,A)+o(1)}.
- A(σ): infimum of A with N(σ − δ, T) ≪ T^{A(1−σ)+o(1)} (δ > 0 infinitesimal); in particular N(σ, T) ≪ T^{A(σ)(1−σ)+ε}.
- Montgomery conjecture (MC) at (σ, τ): LV(σ, τ) ≤ 2 − 2σ.

Published record on the range of this note (TTY Table 3 = ANTEDB Table "Current best upper bound on A(σ)"):
A(σ) ≤ 3/(2σ) for 1867/2347 ≤ σ < 7/8 (Bourgain 2000; Ivić 1984, Theorem 11.4 of Ivić's book). The record comes from
the Jutila (k = 2) large-value theorem, which gives MC for τ ≤ 3σ − 1, fed into Corollary zero-large-cor3 below.

TTY results used (all classical in substance; references inside TTY):
- [HM] Theorem montgomery-lv (Halász–Montgomery): MC holds at (σ, τ) if every sum Σ_{n∈[N,2N]} n^{iu} with
  N^{1−o(1)} ≤ |u| ≤ T is ≪ N^{2σ−1−δ} for a fixed δ > 0. (Its proof uses only these zeta sums.)
- [cor3] Corollary zero-large-cor3: if MC holds for 0 ≤ τ ≤ τ0 + σ − 1, τ0 ≥ 3 − 3σ, and
  LV_ζ(σ, τ) ≤ (3 − 3σ)τ/τ0 for 2 ≤ τ < 4τ0/3, then A(σ) ≤ 3/τ0. (Zero detection + Huxley subdivision + raising to
  powers. Sanity: MC ranges τ_M = 1, 4σ − 2, 3σ − 1 give, via τ0 = τ_M + 1 − σ, Ingham's 3/(2 − σ), Huxley's
  3/(3σ − 1) and Ivić's 3/(2σ).)
- [lvz-mu] Corollary lvz-mu: LV_ζ(σ, τ) = −∞ whenever σ > σ0 + τ μ(σ0).
- [mad] Lemma mad: M(σ0, A) = sup_{τ≥2, σ≥1/2} (A(σ − σ0) + LV_ζ(σ, τ))/τ (1/2 ≤ σ0 ≤ 1, A ≥ 1).
- [zfl] Lemma zero-from-large: A(σ)(1 − σ) ≤ max(sup_{τ≥2} LV_ζ(σ,τ)/τ, limsup_{τ→∞} LV(σ,τ)/τ).
External theorems: μ(1/2) ≤ 13/84 (Bourgain, JAMS 30 (2017)); M(1/2, 12) ≤ 2 (Heath-Brown 1978).

## 1. The three consequences of (H)

(C1) μ(σ) = 0 for σ ≥ B.
Proof (Titchmarsh §14.2 with 1/2 replaced by B). Fix σ > B, δ := (σ − B)/3, t large. The disc D around 2 + it of
radius R = 2 − B − δ lies in Re s ≥ B + δ, contains no zero (H) and, for t ≥ 3, not the pole; so log ζ has a
holomorphic branch on D with |log ζ(2 + it)| ≤ log ζ(2). On D, Re log ζ = log|ζ| ≤ C log t (polynomial growth of ζ).
Borel–Carathéodory on the concentric disc of radius r = 2 − B − 2δ gives |log ζ(s)| ≤ C_δ log t there. Hadamard's three
circles, centred at σ1 + it with σ1 = 2, radii r1 = 1 − η (where |log ζ| ≤ log ζ(1 + η)), r2 = 2 − σ and
r3 = 2 − B − 2δ, gives |log ζ(σ + it)| ≤ C (log t)^a with a = log(r2/r1)/log(r3/r1) < 1 (choose η small so that
r1 < r2 < r3). Hence log|ζ(σ + it)| = o(log t), i.e. ζ(σ + it) ≪ t^ε. So μ(σ) = 0 for σ > B, and μ(B) = 0 because μ is
convex, hence continuous, and μ ≥ 0. ∎

(C1') μ(1 − B) = B − 1/2 (functional equation: |ζ(1 − B + it)| ≍ |t|^{B−1/2}|ζ(B − it)|), and μ(u) ≤ B − 1/2 for
u ≥ 1 − B (μ non-increasing).

Lemma 1 (zeta sums of medium length). Assume (H). For every ε > 0, N ≥ 2, interval I ⊂ [N, 2N], and real t with
πN ≤ |t|:  Σ_{n∈I} n^{−it} ≪_ε N^{1−B} |t|^{B−1/2+ε}.
Proof. Truncated Perron (Titchmarsh Lemma 3.12) with c = 1 + 1/log N and height U = |t|/2, for x ∈ [N, 2N + 1] a
half-integer: Σ_{n≤x} n^{−it} = (1/2πi)∫_{c−iU}^{c+iU} ζ(w + it) x^w dw/w + O(N log N/U). Move the segment to
Re w = 1 − B. The pole of ζ(w + it) is at w = 1 − it, outside the rectangle since |Im w| ≤ |t|/2 < |t|. On the new
vertical side |ζ(1 − B + i(t + v))| ≪ |t|^{B−1/2+ε} for |v| ≤ |t|/2 (C1'), contributing ≪ N^{1−B}|t|^{B−1/2+ε} log|t|.
On the horizontal sides Im(w + it) = t ± t/2 stays ≍ |t| and |ζ| ≪ |t|^{B−1/2+ε} (C1'), contributing
≪ N |t|^{B−3/2+ε} ≤ N^{1−B}|t|^{B−1/2+ε} (as N^B ≤ |t|). The Perron error N log N/U ≪ N^{1+ε}/|t| is of the same
or smaller size. Subtract two such x (and adjust O(1) endpoint terms). ∎
In TTY's exponent language, with |t| = N^{τ'} (1 ≤ τ' ≤ 2): the sum is ≪ N^{g_B(τ')+o(1)}, g_B(τ') = 1 − B + (B − 1/2)τ'.
This is the B-process image (B − 1/2, 1/2) of the "zeta exponent pair" (0, B) given by (C1); it is valid for zeta sums
only (not a general exponent pair). For |t| ≤ πN the first-derivative test gives ≪ N/|t| as usual.

(C2) μ(a) ≤ (13/84)(B − a)/(B − 1/2) for 1/2 ≤ a ≤ B: convexity of μ between μ(1/2) ≤ 13/84 and μ(B) = 0.

(C3) If M(s1, m1) = 1 (1/2 ≤ s1 < B, m1 ≥ 1), then M(s, m) = 1 for s1 < s < B and m = m1(B − s1)/(B − s).
Proof via [mad]. M(s1, m1) = 1 gives LV_ζ(σ, τ) ≤ τ − m1(σ − s1) for all σ ≥ 1/2, τ ≥ 2; (C1) and [lvz-mu] give
LV_ζ(σ, τ) = −∞ for σ > B. For σ ∈ [1/2, B]: m(σ − s) + LV_ζ(σ, τ) ≤ τ + (m − m1)σ − ms + m1 s1, increasing in σ
(m > m1), so at most its value at σ = B, which is τ + m(B − s) − m1(B − s1) = τ. Hence the sup in [mad] is ≤ 1. ∎
(Exact form of the classical tool, which the proof above does not need: for F holomorphic in α ≤ Re s ≤ β and
decaying as |Im s| → ∞, ‖F(γ + i·)‖_r ≤ ‖F(α + i·)‖_p^{1−θ}‖F(β + i·)‖_q^θ with θ = (γ − α)/(β − α) and
1/r = (1 − θ)/p + θ/q (Gabriel's convexity theorem; proof by duality and the three-lines lemma). Applied to
F(s) = ζ(s)·exp((s − 3iT/2)²/T²) with q = ∞ at β = B it gives (C3) again. Ivić's moment table (his Theorem 8.4) is of
this kind: its entries 4/(3 − 4σ), 10/(5 − 6σ), 19/(6 − 6σ) have 1/m affine in σ. The variant that interpolates m
itself linearly is false for p ≠ q (replace F by F·e^{λs} and let λ → ±∞). We could not access Ivić's book to quote
his Lemma 8.3 verbatim.) Example (verify.py): from Ivić's m(3/4) = 14.27, (C3) gives m(0.78) ≥ 18.85 versus Ivić's 17.17.

## 2. What does not help (and why the obvious routes fail)

Proposition 2. Let 1/2 < σ < B.
(a) Through [lvz-mu], the new μ data (C2) on 1/2 ≤ a ≤ B and its hull with old bounds give no new τ with
LV_ζ(σ, τ) = −∞: on a segment μ*(a) = c(B − a) (c > 0) of the new hull, (σ − a)/μ*(a) increases as a decreases
(derivative sign σ − B < 0), so max_a (σ − a)/μ*(a) is attained at an old point. The same holds for Corollary mu-lv
(MC for τ ≤ (2σ − 1 − a)/μ(a)) with a ∈ [1/2, B]: its target x = 2σ − 1 is < B. It fails on the reflected side
a ∈ [1 − B, 1/2], where the new point (1 − B, B − 1/2) lies far below the old value μ_old(B) + B − 1/2; Theorem 1
below is exactly Corollary mu-lv at a = 1 − B (whose proof needs [lvz-mu] at σ0 = 1 − B < 1/2, i.e. Lemma 1).
(b) Through [mad], every LV_ζ bound obtainable from (C3), LV_ζ(σ, τ) ≤ τ − m(s)(σ − s), is weaker than the one at s1,
because m(s)(σ − s) = m1(B − s1)(σ − s)/(B − s) ≤ m1(σ − s1) ⟺ (s − s1)(σ − B) ≤ 0. More generally, data
(s, 1/m, M/m) interpolated linearly towards (B, 0, 0) give m(Mτ/m − σ + s) = m1[(ν1τ − σ + s1) + θ(B − σ)/(1 − θ)],
increasing in θ; with finite q at B the bound is linear-fractional in θ, so extremal at θ ∈ {0, 1}, and θ = 1 gives
τ + q(B − σ) ≥ τ.
(c) Hence the LV_ζ(σ, ·) part of [zfl] at σ < B is unchanged by (C1)–(C3): the only new zeta large-value information
for τ ≥ 2 is LV_ζ(σ', τ) = −∞ for σ' > B. Inside the Halász–Montgomery inequality (Lemma hl-improv) the levels σ' > B
are removed for τ' ≥ 2, but every known bound has σ' + LV_ζ(σ', τ') decreasing in σ' (each term of the minimum has
slope ≤ −4, e.g. the fourth moment), so the supremum sits at σ' = 2σ − 1 < B and nothing changes there either.
So the moment route proposed in round R348 gives nothing. The gain comes from τ' ∈ (1, 2), i.e. zeta sums longer
than T^{1/2}: there (C1) reflects to Lemma 1, which beats all known bounds (Section 4).

## 3. The gain

Theorem 1 (longer Montgomery range). Assume (H). For 3/4 < σ < (1 + B)/2 and 0 ≤ τ < τ_B(σ),
  τ_B(σ) := (2σ − 2 + B)/(B − 1/2),
one has LV(σ, τ) ≤ 2 − 2σ.
Proof. For τ ≤ 1 this is the L² mean value theorem. Let 1 < τ < τ_B(σ) ≤ 2 and run [HM]: we must bound
Σ_{n∈[N,2N]} n^{iu} for N^{1−o(1)} ≤ |u| ≤ T = N^{τ+o(1)}. If |u| ≤ πN, the first-derivative test gives ≪ N^{o(1)}.
If πN ≤ |u| ≤ T, Lemma 1 gives ≪ N^{g_B(τ)+O(ε)} with g_B(τ) = 1 − B + (B − 1/2)τ, since g_B is increasing. Now
g_B(τ) < g_B(τ_B(σ)) = 2σ − 1 (verify.py (b)), and 2σ − 1 > 1/2 = g_B(1). Choosing ε small, all far terms are
≪ N^{2σ−1−δ} with fixed δ > 0, and [HM] gives LV(σ, τ) ≤ 2 − 2σ. ∎

Theorem 2 (zero density). Assume (H). Put σ* := (5 − 4B)/(7 − 6B) = 0.856092011887… For σ* ≤ σ < B,
  A(σ) ≤ 3/(τ_B(σ) + 1 − σ) = 3(B − 1/2) / ((5/2 − B)σ − (5/2 − 2B)),
that is N(σ, T) ≪_ε T^{3(B−1/2)(1−σ)/((5/2−B)σ−(5/2−2B)) + ε}. The bound equals 3/(2σ) at σ = σ* and is smaller for
σ > σ*. (For σ > B, N(σ, T) = 0 by (H).)
Proof. Fix σ ∈ [σ*, B) and small η > 0; set τ0 := τ_B(σ) + 1 − σ − η. Then (verify.py, on [σ*, B]):
(i) τ0 ≥ 3 − 3σ (check (c)); (ii) MC holds for 0 ≤ τ ≤ τ0 + σ − 1 = τ_B(σ) − η by Theorem 1 (σ* > 3/4 and
B < (1 + B)/2); (iii) for 2 ≤ τ < 4τ0/3 we have τ < (σ − 1/2)·84/13 (check (e), with margin ≥ 0.018), so
LV_ζ(σ, τ) = −∞ by [lvz-mu] with σ0 = 1/2 and Bourgain's 13/84. [cor3] gives A(σ) ≤ 3/τ0; let η → 0 (A(σ) is an
infimum). The closed form is check (g); the comparison with 3/(2σ) is τ_B(σ) ≥ 3σ − 1 with equality exactly at σ*
(checks (h) and "tauB(sigma*) = 3 sigma* − 1"). Alternative to (iii): Heath-Brown's LV_ζ ≤ 2τ − 12(σ − 1/2) suffices
because τ0 ≤ 3σ − 3/4 (check (f)). ∎

Values (verify.py): A(0.86) ≤ 1.73489 (record 1.74418); A(0.865) ≤ 1.71331 (1.73410); A(0.87) ≤ 1.69226 (1.72413);
A(0.8725) ≤ 1.68193 (1.71919); A(B−) ≤ 3(B − 1/2)/(9B/2 − B² − 5/2) = 1.678315 (3/(2B) = 1.717461);
A(B−)(1 − B) = 0.212504.

General half-plane. The same proof under "ζ ≠ 0 for Re s > θ" (3/4 < θ < 1) gives A(σ) ≤ 3(θ − 1/2)/((5/2 − θ)σ −
(5/2 − 2θ)) on [(5 − 4θ)/(7 − 6θ), θ); use Heath-Brown's twelfth moment for (iii): τ0 ≤ 3σ − 3/4 holds at σ = σ*(θ)
because σ* > 3/4, and at σ = θ because 4θ² − 27θ/4 + 23/8 > 0 (negative discriminant); the 13/84 route can fail at the
left end for smaller θ (e.g. θ = 0.8). For θ = 7/8 (OpenAI 2026) the range is [6/7, 7/8). With Bebel's newer
T_* = B_2 + 2.2·10^−254 in place of B the bound improves by O(10^−254) (immaterial).

## 4. Is this really better than everything known? (comparison)

(1) Literature large-value theorems recorded in the ANTEDB (L², Huxley, Heath-Brown, Jutila k = 1…19, Guth–Maynard
and its k-variants, Bourgain's optimized polytopes, Heath-Brown's second theorem, Ivić's): their best MC range on
[0.80, 0.875] is Jutila's 3σ − 1 (repro/explore.py reproduces 3/(2σ) on [0.80, 0.869] with all of them plus
raising to powers and subdivision).
(2) [HM] with the best known exponential-sum exponents β (ANTEDB compute_best_beta_bounds, including TTY's 2025
exponent pairs): on τ ∈ [1.485, 1.68] ⊃ [τ_B(σ*), τ_B(B)] = [1.568, 1.661] the old zeta-sum exponent
g_old(τ) = 1 − τ/2 + τβ(1 − 1/τ) is a piecewise linear function with four pieces, each strictly above g_B (verify.py,
exact rationals). So τ_B(σ) exceeds the old HM range. Side remark: the old HM range itself beats Jutila for
σ > 325/374 = 0.86898 (pair (178/641, 365/641) = B-image of TTY's (89/1282, 997/1282)), giving the unrecorded routine
bound A(σ) ≤ 267/(552σ − 325) on (0.86898, 0.88166); Theorem 2 beats it on (0.86898, B) (check (i)), and (H) makes
it moot above B.
(3) ANTEDB's own pipeline (zero_density_estimate.compute_large_value_region + compute_sup_rho_on_tau, i.e. Corollary
11.7 with raising to powers; repro/antedb_check.py) returns 3(1 − σ)/(2σ) without our hypothesis and exactly
7000911873997641(1 − σ)/(10166362708667453σ − 4707725417334906) = A_new(σ)(1 − σ) with it, on windows covering
[0.8561, 0.8625] ∪ [0.8650, 0.8733] (seven runs, τ0 between 2.29 and 2.39; near the right end of a window a fixed τ0
is too small and the code returns a constant (4 − 4σ)/τ0 piece, an artefact of fixing τ0, also present without the
new hypothesis).
(4) Full numerical optimization of [cor-0] over τ0 with all estimates (repro/explore.py, τ-grid 0.004) gives the same
values to within grid error (e.g. 1.7341 at 0.86); a scout's earlier "A(0.86) ≤ 1.72134" is not reproduced.

## 5. Self-review (two passes)

Pass 1 (hidden hypotheses, ranges, directions).
- Lemma 1 needs |t| ≥ πN; smaller |u| in [HM] are handled by the first-derivative test. The horizontal Perron sides
  avoid the pole because U = |t|/2. Sharp cutoffs only cost N^{o(1)}.
- [HM] needs the far-off-diagonal bound uniformly for |u| up to T (not only at height T): g_B is increasing, so the
  worst case is |u| = T. Large-value patterns allow J anywhere; only differences u = t − t' enter.
- Theorem 1 is used only for σ ∈ [σ*, B) ⊂ (3/4, (1 + B)/2) and τ ≤ τ_B(σ) ≤ 2, where Lemma 1 applies.
- [cor3] needs MC on a closed range; we use τ0 − η and let η → 0. τ0 ≥ 3 − 3σ holds.
- The LV_ζ hypothesis of [cor3] is covered by a proven theorem (Bourgain 13/84), with Heath-Brown's twelfth moment as a
  second, independent proof.
- Direction: smaller A is better; A_new < 3/(2σ) ⟺ τ_B(σ) > 3σ − 1 ⟺ σ > σ*.
Pass 2 (inputs that are not theorems; newer records).
- (H) is a Lean-kernel-checked statement (unrefereed but machine-verified); everything else used is refereed or
  classical (TTY's lemmas are standard zero-detection, subdivision and power arguments; TTY itself is a preprint).
- Newer records: arXiv math.NT listings 2026-10 (866 entries) and 2026-11 (141) and arXiv API queries for 2025–2026
  zero-density and large-value papers show no improvement of A(σ) for ζ on [0.856, 0.874] (2610.01035 de Faveri–Pandey
  improves large values of ζ on the critical line, which only affects the LV_ζ side, already −∞ here; 2507.08296
  Chen–Gupta–Li is the q-aspect). ANTEDB main (2026-11-05) still lists 3/(2σ).
- Our own earlier rounds R299/R304 state that the half-plane cannot improve Halász–Montgomery bounds or density below
  it; that is wrong for zeta sums of length > T^{1/2} (Lemma 1). R348's proposed moment route is correctly dead
  (Proposition 2).
