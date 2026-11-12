NO IDEA: Bérard's Ehrenfest-time cutoff, run with the L¹-optimal (Beurling–Selberg) Tauberian kernel, is exactly Carneiro–Chandee–Milinovich's 2015 one-endpoint proof of |S(t)| ≤ (¼ + o(1)) log t/log log t, so microlocal tools reproduce ¼ and nothing better; going below needs uniform cancellation of the orbit sum Σ_{n<e^Δ}Λ(n)n^{−1/2−it} past Δ_E ≈ 2 log log t, which neither field has, and on Odlyzko's windows (2.5·10³–1.4·10²¹) that sum is statistically a random-phase sum (node N9).

# ExploreMicrolocal: microlocal and semiclassical analysis → RH tree (explore47, Orca, 12 Nov 2026)

## Choice (measured before choosing; DERIVATION §1)
- rg (synonyms) in NOGO/Tree/TREE: no microlocal, Egorov, Gårding, Fefferman–Phong, Bérard, Ehrenfest hits; Berry–Keating only in R51,
  prolatefals (N13) and aclgh2's semiclassical count (N13).
- N1, tools (a)+(d), not chosen: Q_c is a Toeplitz form on [−a, a], a = ½ ln c, whose multiplier is band-limited to exactly [−2a, 2a]
  (n < c ⇔ log n < 2a). Its symbol oscillates on the window's own Heisenberg scale (h_eff ≍ 1 for every c), its Weyl symbol jumps at
  |x| = a − ½ log n, and near ξ = 0 it is negative already for c < 2, where Yoshida proves positivity (the pole is rank 2, not a symbol).
  Sharp Gårding (−Ch), Fefferman–Phong (−Ch²) and anti-Wick (≥ inf σ) errors are O(1) against margins ε_c ≈ e^{15.9}c⁵e^{−4πc}. Measured
  (n1_symbol.py, |ξ| ≤ 400): inf σ_c = −2.54 at c = 2 and −5.05 at c = 47, and σ_c < 0 on 1.1 and 117 phase-space cells, while Q_c has no
  negative eigenvalue. The only exact symbol certificate, plane-wave anti-Wick, is Krein's positive extension ⇔ Q_c ≥ 0 (already in
  explore19's dictionary).
- N13, two-term Weyl law, not chosen: one degree of freedom ⇒ Bohr–Sommerfeld ⇒ N_H − N̄ = O(1), while N_ζ − N̄ = S(T) is unbounded
  (Tsang Ω±; Selberg's mean square (1/2π²) log log T): prolatefals' barrier in general form, known since Berry–Keating 1999.
- Chosen: tool (b), Bérard's Ehrenfest-time cutoff → N9 quantity S(t), with consequences for N1 (Mori's ϖ) and N8 (gap constant).

## Transfer (DERIVATION §2–3)
Compact hyperbolic surface ↔ ζ: Selberg trace formula ↔ Guinand–Weil formula; closed geodesics ℓ, #{ℓ ≤ L} ~ e^L/L, amplitudes ~e^{−ℓ/2}
↔ log n with Λ(n)n^{−1/2}; self-adjointness ↔ RH; Bérard's O(λ/log λ) ↔ S(t) = O(log t/log log t); Ehrenfest time ↔ Δ_E, where the
orbit sum e^{Δ/2} reaches the density log t.

Statement (RH; every t > 0 and Δ > 0). R^± = Selberg's majorant/minorant of χ_[−t,t] of type Δ (primes n < e^Δ; CCM 1503.00955 Lemma 2), η = log n/Δ:
  −log(t/2π)/(2Δ) + O⁻(t,Δ) − p⁻ ≤ S(t) ≤ log(t/2π)/(2Δ) + O⁺(t,Δ) + p⁺,
  O^±(t,Δ) = −(1/Δ)Σ_{n<e^Δ} Λ(n)n^{−1/2}[(1/π + (1−η)cot πη) sin(t log n) ± (1−η) cos(t log n)],
with the τ-term log(t/2π)/(2Δ)·(1 + o(1)) (within 0.1–0.4 % at t = 100 and 1000) and p^± = O(e^{Δ/2}/(Δt)²) (≤ 10⁻⁴ in every test). Bérard's step bounds the orbit sum by
TB(Δ) = Σ|weights| ~ e^{Δ/2}/Δ, so Δ_E = 2 log log t − 2 log log log t and |S(t)| ≤ (¼ + o(1)) log t/log log t. This is CCM's
"new and simple proof" (1503.00955) verbatim (Goldston–Gonek's ½ came from S(t+h) − S(t), two endpoints). Bérard's own kernel
(Hörmander's Tauberian theorem) has more L¹ excess than Beurling's H^±, which are L¹-extremal, so the microlocal version gives ≥ ¼.

## Tests (logs/; every process < 290 MB, peaks 18–228 MB)
- Smallest case (bt_validate.py, first 10⁵ zeros, t = 100 and 1000, Δ = 3–7, both signs): Σ_γ R^±(γ) over the zeros matches the
  explicit-formula side to ≤ 4·10⁻³; the pole term equals 2 to 10⁻⁴; τ_exact is within 0.1–0.4 % of log(t/2π)/(2Δ). Example brackets:
  t = 100, Δ = 5: −0.0178 ≤ S = −0.0024 ≤ 0.0112.
- Scale (bt_window.py): six Odlyzko windows, 2000–2500 zeros each, Δ = 2…14 (93,371 prime powers, phases mod 2π in mpmath), 364,000 bound
  checks at γ_k^±. All hold; minimum margin 0.0037.

| window | t | Δ_tri | B_tri | ¼log t/loglog t | max abs S | at Δ_E: TB / sup O_ζ / sup O_rand | at Δ = 14: TB / sup O_ζ / sup O_rand |
|---|---|---|---|---|---|---|---|
| zeros #1001–3000 | 2.5·10³ | 3.00 | 1.88 | 0.95 | 1.31 | 0.88 / 0.69 / 0.72 | 24.7 / 1.21 / 1.34 |
| #20001–22000 | 1.9·10⁴ | 3.25 | 2.20 | 1.08 | 1.46 | 0.88 / 0.69 / 0.72 | 24.7 / 1.34 / 1.34 |
| #90001–10⁵ | 7.2·10⁴ | 3.50 | 2.39 | 1.16 | 1.63 | 0.88 / 0.70 / 0.72 | 24.7 / 1.32 / 1.34 |
| #10¹² + 1… | 2.7·10¹¹ | 5.00 | 4.13 | 2.01 | 1.94 | 1.69 / 0.98 / 0.97 | 24.7 / 1.38 / 1.34 |
| #10²¹ + 1… | 1.4·10²⁰ | 6.00 | 5.94 | 3.02 | 2.00 | 2.23 / 1.09 / 1.08 | 24.7 / 1.27 / 1.34 |
| #10²² + 1… | 1.4·10²¹ | 6.00 | 6.13 | 3.13 | 2.03 | 2.23 / 1.08 / 1.08 | 24.7 / 1.32 / 1.34 |

B_tri = min_Δ[log(t/2π)/(2Δ) + TB(Δ)] is the best bound the transferred argument gives. Δ_tri is its Ehrenfest time.

## What it gives on the nodes (TREE N9, N1, N8)
- N9: no new bound; CCM's ¼ (2013/2015) stands. At every computable height the method sits at about twice its own leading term:
  B_tri/(log t/log log t) = 0.48–0.52 on [10³, 10¹²] and 0.42 at 10¹⁰⁰ (bt_range.py), and B_tri is 1.4–3.0 × the observed max|S|.
- Where it breaks: the orbit sum, not the Tauberian step. At Δ_E the triangle inequality loses only ×1.3–2.1, and at Δ = 14 it loses
  ×18–20, so the gain has to come from Δ ≫ Δ_E. That needs uniform-in-t cancellation for Dirichlet polynomials of length (log t)^{2/θ},
  which neither field has: we found no improvement of Bérard's rate on any compact negatively curved manifold (Keeler 1905.05136,
  Canzani–Galkowski 2010.03969 and Sussman 2109.09926 widen its hypotheses, not its rate). Near-alignments for primes ≲ log t log log t exist (the resonance method).
- Hypothesis |O^±(t,Δ)| ≪ e^{θΔ/2} for all t (θ < 1): with Δ = (2/θ)log log t one gets S-constant θ/4, multiplicity θ/2, GG gap θπ/log log t,
  and Mori's RH-alone rate λ_min ≥ exp(−exp(μ^ϖ)) for ϖ > θ/2. Today θ = 1, so ϖ > ½, which matches MoriRate's "natural limit of RH +
  S-bounds"; ϖ = θ/2 is our reading of his cluster step and was not rechecked in his text.
- Explicit unconditional statement in a range (asked): routine. RH is verified to 3·10¹², so the inequality holds unconditionally for
  t ≤ 10¹² up to ≈ 10⁻¹⁰. There B_tri(10¹²) = 4.26, against Trudgian 2014's 6.53 and HSW 2022's 13.09. Not certified (τ, p^± and the
  off-line term need explicit constants); at most a range record (RECORDS.md).

## Adversarial
- The statement assumes RH (γ real) and uses Λ(n)n^{−1/2} weights. Davenport–Heilbronn violates both premises, so the statement is a
  consequence of RH and carries no RH information.
- Random-phase model (iid θ_p per prime, p^k ↦ kθ_p, same sample counts): at every Δ = 2–14, ζ's sup orbit sums are within −11 % … +5 %
  of the model's and its rms within −3 % … +14 % (the +14 % is at t = 2.5·10³, Δ = 14, where ζ is sampled at its own zeros). Examples:
  sup 0.69 vs 0.72 at 2.5·10³ and 1.08 vs 1.08 at 1.4·10²¹. The slack above max|S| is generic Bohr–Jessen cancellation, not arithmetic.
- GUE and Poisson spectra with ζ's density have no prime orbit sum, so the transferred inequality does not apply to them. For the
  orbit sum, the matching adversary is the random-phase model above.

## Novelty
Nothing is claimed new. The statement is in CCM 1309.1526 and 1503.00955 (explicit RH versions in 2010.13307; GG 2007). The dictionary
is in Hejhal (LNM 548, 1976), Randol 1978 and Bérard 1977; for the Ω side see Jakobson–Polterovich–Toth math/0612250, and for
modular-surface Weyl-remainder data Then 1212.3149. arXiv API searches on 12 Nov (Ehrenfest / Bérard / Weyl remainder × zeta,
S(t) × Selberg, Gårding × RH) found no transfer beyond these.

## Files
DERIVATION.md (choice, dictionary, statement, tests, range, N1 scan); btlib.py; bt_validate.py; bt_window.py; bt_summary.py; bt_range.py;
n1_symbol.py; logs/bt_validate.txt, logs/bt_zeros*.txt|json, logs/bt_summary.txt, logs/bt_range.txt, logs/n1_symbol.txt. The Odlyzko tables (zeros1, zeros3–5 from
www-users.cse.umn.edu/~odlyzko/zeta_tables/) were downloaded to .scratch/explore47 and then deleted. To rerun, fetch them into
./zeros/ (curl -s -L -o zeros/zeros1 <url>/zeros1, likewise zeros3–5), then run bt_validate.py, bt_window.py zeros1 1000 3000 | 20000 22000 |
90000 100000, bt_window.py zeros3 | zeros4 | zeros5, bt_summary.py, bt_range.py and n1_symbol.py.

## Suggested NOGO line
- explore47 (Orca, exploration: microlocal and semiclassical analysis — Bérard's Ehrenfest-time Weyl remainder, sharp Gårding/Fefferman–Phong, anti-Wick, two-term Weyl laws → N9; N1, N13 measured) NO IDEA: Bérard's cutoff with the L¹-optimal kernel is CCM's 2015 one-endpoint proof of |S(t)| ≤ (¼+o(1))log t/log log t (Δ_E ≈ 2 log log t); below ¼ needs uniform cancellation of Σ_{n<e^Δ}Λ(n)n^{−1/2−it} past Δ_E (the same limit sets Mori's ϖ > ½ and GG's gap π); on Odlyzko windows 2.5·10³–1.4·10²¹ (364,000 checks) the transferred bound holds, is 1.4–3.0 × max|S| and its orbit sums are random-phase-generic; Gårding/FP/anti-Wick have no small parameter on Weil windows (h_eff ≍ 1) [ExploreMicrolocal](Riemann/Deep/ExploreMicrolocal/REPORT.md)

RH relevance: N9 (S(t) constant), with consequences for N1 (Mori's ϖ) and N8 (GG gap constant).
