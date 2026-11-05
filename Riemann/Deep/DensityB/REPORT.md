RECORD: A(sigma) <= 3(B − 1/2)/((5/2 − B)·sigma − (5/2 − 2B)) = 1.12014589983962256/(1.62661803338679248·sigma − 0.75323606677358496) for (5 − 4B)/(7 − 6B) = 0.856092… < sigma < B, B = 0.87338196661320752 (with the Lean-checked half-plane Re s > 0.87338196661320752)

# DensityB report (2026-11-05)

In words: N(σ, T) ≪ T^{A(1−σ)+ε} with this A beats the published record A = 3/(2σ) (Ivić 1984, Bourgain 2000; still the
entry in Tao–Trudgian–Yang's table and in the ANTEDB on 2026-11-05) on 0.8561 < σ < B, by up to 2.3% at σ → B.
Proofs: [PROOF.md](PROOF.md). Certificate: [verify.py](verify.py) (exact rationals, ALL PASS, 0.12 s). Numerics and
the ANTEDB cross-check: [repro/](repro/). Paper proofs, not refereed, no Lean.

## How it works (one paragraph)
The half-plane gives Lindelöf on Re s = B (C1). Its mirror image under the functional equation, |ζ(1 − B + it)| ≪
|t|^{B−1/2+ε}, bounds zeta sums of length N ∈ [T^{1/2}, T] by N^{1−B}T^{B−1/2+ε} (a "zeta exponent pair"
(B − 1/2, 1/2)), better than every known bound there. In the Halász–Montgomery inequality this extends the range of
the Montgomery large-value conjecture from T ≤ N^{3σ−1} (Jutila) to T < N^{τ_B(σ)}, τ_B(σ) = (2σ − 2 + B)/(B − 1/2);
the standard zero-detection corollary (TTY Cor. zero-large-cor3) turns this into A(σ) ≤ 3/(τ_B(σ) + 1 − σ). The
zeta-sum side of that corollary is empty here by Bourgain's μ(1/2) ≤ 13/84. The three routes named in the brief,
(C1) at Re s = B, the convexity line (C2) and moment interpolation (C3), are all proved but give nothing below B
(Proposition 2): they only add information at levels above B, and the gain needs the reflected line Re s = 1 − B.

## Table
| σ range | old record (source) | new bound | derivation |
|---|---|---|---|
| [0.7954, 0.856092] | 3/(2σ) (Bourgain 2000; Ivić 1984 = Ivić's book Thm 11.4) | unchanged | τ_B(σ) ≤ 3σ − 1 there |
| (0.856092, B) | 3/(2σ) (same; TTY Table 3, ANTEDB) | 3(B − 1/2)/((5/2 − B)σ − (5/2 − 2B)) | Lemma 1 (reflected half-plane) → Thm 1 (Montgomery conjecture for τ < τ_B) → Thm 2 via TTY Cor. zero-large-cor3, τ0 = τ_B + 1 − σ; LV_ζ part −∞ by μ(1/2) ≤ 13/84 (or HB's 12th moment) |
| [B, 1) | 3/(10σ − 7), 24/(30σ − 11), … | −∞ (N(σ, T) = 0) | the half-plane itself |

Values: A(0.86) ≤ 1.73489 (was 1.74418); A(0.865) ≤ 1.71331 (1.73410); A(0.87) ≤ 1.69226 (1.72413);
A(0.8725) ≤ 1.68193 (1.71919); just below B: A ≤ 3(B − 1/2)/(9B/2 − B² − 5/2) = 1.678315 (3/(2B) = 1.717461),
i.e. N(σ, T) ≪ T^{0.21250+ε} as σ → B−.

## Answers to the brief's questions
- Record improved? Yes, on (0.856092, B), as above. The best unpublished "old" bound there (Halász–Montgomery with the
  B-image (178/641, 365/641) of TTY's 2025 exponent pair, which beats Jutila for σ > 325/374 = 0.86898 and gives
  A ≤ 267/(552σ − 325)) is also beaten (verify.py check (i)).
- Density-hypothesis range: unchanged (A ≤ 2 already holds for σ ≥ 25/32; the gain starts at 0.856).
- Just below B: A(σ) ≤ 1.6783 and A(σ)(1 − σ) → 0.2125 (old 0.2175); above B it is −∞. No continuity at B.
- (C1)–(C3): proved (PROOF.md §1; (C3) via TTY's moment/large-value duality Lemma mad; Gabriel's interpolation is the
  version with 1/m linear, the brief's formula). They give no gain (PROOF.md Prop. 2): the interpolated moment
  m(s) = m1(B − s1)/(B − s) yields LV_ζ ≤ τ − m(s)(σ − s), always weaker than at s1 when σ < B. Example of (C3) itself:
  m(0.78) ≥ 18.85 vs Ivić's 17.17 (newer: de Faveri–Pandey 2610.01035 improve moments at σ = 1/2 for 4 < k ≤ 12;
  irrelevant here).
- The scout's A(0.86) ≤ 1.72134: not reproduced; the full numerical optimization gives 1.7349 (= Theorem 2 up to grid
  error).

## Cross-checks
- ANTEDB's own code (Corollary 11.7 large-value computation, all literature large-value estimates, raising to powers):
  without our hypothesis it returns 3(1 − σ)/(2σ); with it, exactly
  7000911873997641(1 − σ)/(10166362708667453σ − 4707725417334906) = A_new(σ)(1 − σ) on windows covering
  [0.8561, 0.8625] ∪ [0.8650, 0.8733] (repro/antedb_check.py; one fixed τ0 per window).
- Own pipeline (repro/explore.py): reproduces 3/(2σ) on [0.80, 0.869] without the hypothesis (sanity) and the new
  values with it.
- The old best zeta-sum exponents (ANTEDB best β, 4 linear pieces on τ ∈ [1.485, 1.68]) lie strictly above the new
  line g_B(τ) = 1 − B + (B − 1/2)τ (exact check in verify.py).
- The formula A = 3/(τ_M + 1 − σ) reproduces Ingham (τ_M = 1), Huxley (4σ − 2) and Ivić (3σ − 1).

## verify.py output (uv run --no-project python3 verify.py; runtime 0.124 s)
```
PASS B above Liu's contact B2 = (777-sqrt1081)/852 (0 < 777-852B and (777-852B)^2 < 1081)
PASS B < 7/8
PASS gB(1) = 1/2
PASS gB(2) = B
B      = 0.87338196661320752
sigma* = 0.856092011887297  (5-4B)/(7-6B)
PASS 3/4 < sigma* < B
PASS B < (1+B)/2 (HM range for Theorem 1 contains [sigma*, B])
PASS (a) 1 <= tauB(s) <= 2 on [sigma*, B]
PASS (b) gB(tauB(s)) = 2s-1 exactly (HM threshold)
PASS (c) tau0 > 3-3s
PASS (d) (lvo): LV <= (3-3s) tau/tau0 on [2tau0/3, tau0]
PASS (e) (lvoz) via mu(1/2)<=13/84: 4tau0/3 <= (s-1/2)*84/13
PASS (f) (lvoz) alternatively via HB twelfth moment: tau0 <= 3s-3/4
PASS (g) 3/tau0 equals closed form 3(B-1/2)/((5/2-B)s-(5/2-2B))
PASS (h) Anew <= 3/(2s), equality only at sigma*
PASS (i) tauB(s) > (641s-414)/89 (beats routine HM with (178/641,365/641))
PASS (j) 4tau0/3 > 2
PASS tauB(sigma*) = 3 sigma* - 1 (crossover with Jutila k=2)
PASS tauB(s) > 3s-1 for s in (sigma*, B]
PASS beta pieces cover tau in [tauB(sigma*), tauB(B)]
PASS g_old(tau) > gB(tau) on [1.485, 1.68] (new zeta-sum bound strictly better there)
C2 mu(0.8) <= 0.030415  mu(0.85) <= 0.009691
C3 m(3/4) = 14.2702  m(0.78) >= 18.8547  Ivic m(0.78) = 17.1665
PASS C3 example m(0.78) >= 18.8 > Ivic 17.1
PASS C3 no-gain: (B-s1)/(B-s)*(sig-s) <= sig-s1 for s1<s<sig<B
PASS C2 no-gain: (sig-a)/mu_line(a) maximal at the left end a=1/2

sigma        old 3/(2s)     routine-HM 267/(552s-325)   new 3/tau0       tauB
  0.856092      1.75214          (n/a)               1.75214   1.56827
  0.860000      1.74418          (n/a)               1.73489   1.58920
  0.865000      1.73410          (n/a)               1.71331   1.61599
  0.868983      1.72615          (n/a)               1.69649   1.63733
  0.870000      1.72413        1.71991               1.69226   1.64277
  0.872500      1.71919        1.70476               1.68193   1.65616
  0.873381      1.71746        1.69948               1.67831   1.66088
A(B-) <= 1.678315  = 3(B-1/2)/(9B/2 - B^2 - 5/2);  3/(2B) = 1.717461 ; A(B-)(1-B) = 0.212504
PASS A(B-) closed form

runtime 0.124 s
ALL PASS
```

## Novelty pass
- arXiv math.NT listings 2026-11 (141 entries) and 2026-10 (866), titles searched for zero density, zero-free,
  half-plane, 7/8, large value, Dirichlet polynomial, Lindelöf, moment. Read: 2610.01035 (sixth moment; large values on
  the critical line only), 2610.15989 (short strips), 2610.18923 (trilinear sums), 2610.12234, 2610.15291, 2610.26835,
  2611.00924 (half-plane constants, no density consequence), 2611.03356 (Landau's local method), 2611.05000,
  2610.37371. None states a zero-density consequence of a zero-free half-plane below it.
- arXiv API, 2025–2026: "zero density"/"zero-density" in math.NT (120 newest), "large value(s)" + "Dirichlet
  polynomial": 2508.02041 (near σ = 1), 2507.08296 (q-aspect), 2609.22624, 2607.04632 (survey), 2609.17875,
  2509.09771, 2503.07410: no improvement of A(σ) for ζ on our range.
- Zenodo ("zero density"; "zero-density AND zeta", 25 newest each): nothing relevant.
- GXL evidence API (1 query, "zero-free half-plane improves zero density below it via Halász–Montgomery and the
  functional equation"): 0 supporting; the 2 "refuting" passages are about Beurling zeta functions (no functional
  equation), so they do not apply.
- ANTEDB (main, 2026-11-05): 3/(2σ) on [0.7954, 0.875); no lemma deriving density from a zero-free half-plane.
- Our files: DISCOVERIES.md, COROLLARIES.md, discoveries/, REFUTATIONS.md: nothing similar. Our rounds R299 ("the 7/8
  half-plane … Halász–Montgomery bounds do not improve: HOLDS") and R304 ("no mechanism improves density just below
  θ = 7/8") say the opposite; both are wrong for zeta sums longer than T^{1/2}. R348's moment route is dead
  (Proposition 2) and should be logged in NOGO.md.

## Significance (for the reviewer)
A record on node N7 in the literal sense: the 1984 bound 3/(2σ) falls on (0.8561, 0.8734). The method is a short,
classical chain (functional equation + Halász–Montgomery + Huxley-type zero detection); the new ingredient is the
2026 half-plane. Range and gain are small (≤ 2.3% of A, only within 0.0173 below B), and every zero above B is already
excluded, so it does not change the density-hypothesis range or the sup of A that drives primes in short intervals.
For any half-plane Re s > θ (3/4 < θ < 1) the same proof gives A(σ) ≤ 3(θ − 1/2)/((5/2 − θ)σ − (5/2 − 2θ)) on
[(5 − 4θ)/(7 − 6θ), θ); for OpenAI's θ = 7/8 the range is [6/7, 7/8). Recommendation: DISCOVERIES.md only if the
reviewer accepts a record by routine method; otherwise COROLLARIES.md. Correct R299/R304 in the tree either way.

## Leads
- Use Lemma 1 inside averaged large-value arguments (Jutila's k-th powers, Heath-Brown's double zeta sums, Bourgain's
  lemma) instead of only pointwise; a longer Montgomery range there would lower the threshold 0.856.
- The Lean-checked half-plane is uniform for Dirichlet L-functions: L(1 − B + it, χ) ≪ (q(|t| + 1))^{B−1/2+ε}, so the
  same Halász–Montgomery step should lengthen the q-aspect Montgomery range for Σ_χ N(σ, T, χ) just below B.

RH relevance: N7 (zero-density exponents: how many zeros can lie off the critical line).

Runtime of this agent: about 1 h (19:08–20:05 UTC). Downloads deleted.
