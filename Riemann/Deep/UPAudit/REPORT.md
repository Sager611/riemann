AUDIT PASSED — Theorem U′ (D(λ) > 0 for every λ > 7/11, unconditionally) survives a line-by-line referee: complete proofs of Lemma Z, Lemma M, Step 6e′ and Step 7′ (FULLPROOF.md) close every step; four cosmetic corrections, no error.

# upaudit — line-by-line referee of Theorem U′ (dpcbeyond, Riemann/Deep/DPCBeyond/PROOF.md)

## Short answer

The proof is correct at θ ∈ (1/φ − 1, min{κ/2, 1}) with Conrey's κ = 8/7 − δ, so every φ > 7/11 works.
- **Lemma Z** (rough Möbius polynomial bounded pointwise through the zeros) is correct as stated. Its new idea checks out. The cutoff f(Σu_i/L) stays inside the r-fold integral, so an r-fold product of prime sums costs X^{max d} once, not X^{r·max d}.
- **Lemma M, 6e′, 7′.** Every exponent and constant checks. With A = 12 fixed, no constant chosen later feeds back. The "sketch-level" parts (M4, 6e′(c)) now have complete proofs.
- **Unchanged steps of Theorem U.** None uses θ < ½. DPC's post-shift main term needs only T^θ ≤ T/L^{20}, i.e. θ < 1.
- **Inputs.** Conrey's slope-8/7 density holds uniformly for σ ≥ ½ with a single log T (Baluyot audit). Fujii's moments and Lamzouri's DoubleSumZeros are quoted correctly.
- **Numerics.** Nothing falsifies the proof. The identity behind Lemma Z holds to 10⁻¹² with real zeros, planted deep clusters cost only a bounded factor, and the toy moments diverge exactly past q = κ/θ.

## Brief items

| Item | Verdict | FULLPROOF.md |
|---|---|---|
| (1)(i) Lemma Z | correct, complete proof | §1 |
| — explicit formula (no RH) | correct; contour to Re w = −∞, residues 1, −m_ρ, −1 | Z1 |
| — inversion of squarefree n | correct. Newton power sums replace the set-partition Möbius inversion and give the same coefficients (−1)^ℓ/z_λ, all of modulus ≤ 1. Collision packets carry ρ_η^h q^{−hs}. | Z2 |
| — prime-power packets | correct: ρ_η(h log q/L)/h, support T^{η/h} < q < T^{b/h}, as reviews 1–2 had repaired | Z1, Z3 |
| — X^{max d}, not X^{r·max d} | correct: Σu_i ≤ log(X/𝔟) on the support, then IBP gives |I| ≤ C(X/𝔟)^{M_z}∏κ_A(L|z_i|) | Z4 |
| — R* bound | correct: 𝔟^{−½+y}(X/𝔟)^{y+max d⁺} ≤ X^y𝔟^{−½}X^{max d⁺} | Z5 |
| — gamma factor, pole and trivial zeros | negligible: one pole entry costs ≤ X·2T(LT/2)^{−A} ≤ 1 for A ≥ 3 | Z5(b), Z6 |
| — uniformity | uniform: C = C(f, ψ, ω, η, θ, b, A); uses only t ≥ T/2 and 0 ≤ y ≤ ½ | Z6 |
| (1)(ii) Lemma M, M1–M4 | correct | §2 |
| — Fujii only at comparable heights | holds: the integer trick, restricted to heights [T/4 − 1, 3T], is uniform in m | M1, M2 |
| — M4 | correct: covering, then Hölder twice, giving exponent κ(1 − 1/w)². Large windows need the unconditional count (my first write-up missed this; see R1) | M4 |
| (1)(iii) 6e′ | correct: pointwise bound from (Z±), y-integral with the 1/L gain, all-zero moments, then Hölder over S_a with exponents (v, u, w) and Φ(ε) → 0 | §3 |
| (1)(iv) 7′ and order of choices | correct: c₂ is independent of ε. The order is φ → θ → P → f, η → A = 12 → q, s → u, w, v, r → C → ε → T | §4 |
| (1)(v) θ-check of unchanged steps | no θ < ½ use anywhere. Only 5a and its three uses (6e A₁ and A₂, the 7 Hölder) needed it, and U′ replaces all of them. DPC Prop 1 (l.274, 368–377) and Lemma 6 (l.673, 702) need θ < 1 | §5 |
| (2) Numerics | no counterexample; N-A to N-E below | §6 |
| (3) Inputs | correct, including uniformity, σ − ½ ≪ 1/log T, N(σ, 3T) and multiplicity | §7 |

## Numerical falsification (all under memcap 290 MB; peak 33 MB)

- **N-A, identity (real zeros).** Tested at T = 10⁴ and 1.8·10⁴, three (θ, η, b) settings, up to 5 prime factors, prime-power packets 4, 8, 9, 16, 25, 27, …, and s = ½ ± y + it.
  - Newton power sums vs direct M(s): error 3·10⁻¹³.
  - With the prime-power binomial expansion: 5·10⁻¹².
  - Zero side with the cutoff inside: 10⁻⁷, set by the truncation of the zero window.
  - Control: with no zeros in the table, the error is O(|M|).
  - The agent's n3 had checked only Λ-weighted sums with j ≤ 2. N-A is the first test of the full inversion.
- **N-B, bound (real zeros).** sup |M|/bound stays at 0.4–2.7, flat in T from 10³ to 10⁴. This is a weak test: these polynomials have ≤ 76 terms.
- **N-C, toy with planted off-line zeros** (the zero side of Z1 without primes, scale-invariant in units of 1/L; L = 20, 30, 44; depth a ≤ 0.48L = 21).
  - Three coincident pairs cost e^{0.9…1.6} and five clustered pairs e^{1.5…2.0}. If the contributions multiplied, the cost would be X^{2d} or X^{4d}, up to e^{47} at a = 21.
  - The ratio to Lemma Z's bound is ≤ 0.41 (A = 2) and ≤ 8.4 (A = 4), and does not grow with a.
  - With A = 12 the ratio climbs to 4.4·10³. This does not contradict the lemma: the constant C(A) of Gevrey-class bumps is huge.
- **N-D, M3 threshold in the toy.** The excess q-moment of one zero of depth a grows at rate ≈ 0.69·qθ at a ≈ 17, approaching qθ.
  - Under the profile e^{−κa}, κ = 8/7: the bin series converges for q ≤ 2.5, is borderline at q = 3 and diverges at q = 4.
  - So q < κ/θ is the asymptotic threshold, not just a proof artefact. This is heuristic sharpness, consistent with M3.
- **N-E, dense random configuration in the toy.** 354 random pairs drawn from the profile e^{−κa}, plus three forced deep pairs (a = 8.5, 11.5, 14.5).
  - The ratio to Lemma Z's bound is ≤ 1.7 (A = 2).
  - The q = 2 moment in the window of the deepest pair matches N-D's isolated-pair prediction to within 4 %, so the per-zero contributions add up and clustering does not amplify them.

## Corrections (none changes the theorem)

C1: supp ρ_η ⊂ [η, b] (DPC l.673). U's "> T^{η/2}, ≤ ⌊2θ/η⌋" is weaker but true. C2: M4's power of |m| is immaterial. C3: the letter q is used for two objects. C4: the Ñ window, the same as UAudit d2, is harmless here because Ψ_ρ counts the near zeros. R1, R2 are gaps in my own first write-up that the reviewer found, and both are repaired: M4a for unbounded windows; the far Poisson shells need q > 1, which holds.

## Review

I ran one `reviewer` on the most delicate step, M4 plus 6e′(c)–(d) pasted inline. Its verdict was GAP in my write-up (R1). I repaired it, and the re-review returned M4a CORRECT, with the far part correct for q > 1 (R2, adopted). The reviewer also confirmed the final Hölder bookkeeping and the threshold θ < κ/2.

## Residual risks

- **Inputs used as published.** These are Conrey's density via Baluyot's thesis (separately audited, Riemann/Deep/Baluyot), Fujii's moments as quoted by Lz (unconditional), and DPC Prop 1 and Lemma 6 (audited in DPCAudit and UAudit, and re-read here at θ ∈ (½, 4/7)).
- **Asymptotic constants.** Lemma Z is asymptotic, with C(A) astronomically large for smooth Gevrey bumps. No computation at reachable heights can see the theorem's regime. N-A tests the identity exactly, and N-C/N-D test the mechanism only in the toy.
- **Weaker fallback.** If Conrey's slope were unavailable, Jutila's κ = 1 − δ gives λ > 2/3 (Theorem U), and any κ ∈ (1, 2] gives λ > 2/(2+κ).
- **Novelty** was not re-assessed; it was not in the brief.

## Files

FULLPROOF.md (§0–§8). Scripts: na_identity.py, nb_ratio.py, nc_toy.py, nd_moments.py, ne_config.py. Logs: logs/na_identity.{log,json}, nb_ratio.{log,json}, nc_toy.{log,json}, nd_moments.{log,json}, ne_config.{log,json}.

RH relevance: N8
