THEOREM: (RH) Σ_{n≤N} R(n)² ≥ (0.0558 − o(1)) N^{3/2}(log N)^{1/2} and Σ_{N<n≤2N} R(n)² ≪ N^{3/2}(log N)^{1/2} log log N (Montgomery's pair-correlation theorem), where R(n) = λ_n − (n/2)(log n − 1 − log 2π + γ); hence R(n) = Ω((n log n)^{1/4}) unconditionally, the smoothed mean of R(n) is exactly 7/4, and RH ⟺ Σ_m m^δ y_m² < ∞ for some δ < 1/2, the weight m^{1/2} never being summable.

# LiRate report (node N3, rate of Li's criterion), lirate agent, 2026-11-11

## Short answer
Under RH the Li remainder is, in mean square, of size N^{1/4}(log N)^{1/4} (two-sided, up to log log N); it is not
o((n log n)^{1/4}) (this part holds without RH); its mean is 7/4 (not 0, not 3/4); pointwise only O(√n log n/√(log log n)) is
provable, the gap to n^{1/4} being an exponential sum over zeros below √n that RH does not control. The prime-side oscillation
−π^{−1/2}Λ(m)m^{−1/2}(log m)^{−3/4} n^{1/4} cos(2√(n log m) − 3π/4) is confirmed in the data for n ≤ 10⁵ to 3–4 digits, phase included.
A published conjecture (Merlini–Sala–Sala, |λ_tiny(n)| ≤ a log n, a ≈ 2) is false: first failure n = 5083; false for every a.

## Results (proofs: proofs.md; heuristics: derivation.md)
| # | Statement | Hypothesis | Status |
|---|---|---|---|
| Lemma A | λ_n = (n/2)(log n − κ) + 7/4 + [−2 Re Σ_{γ≤T} e^{inθ_γ} + 2S(T)cos nθ(T) − ∫_T^∞ S g_n' + O(T² log T/n)]; the smooth part is t(n) − 1/4 + O(log n/n), −1/4 = 2·(−π/8)/π | RH | rigorous; reviewer CORRECT (after repair) |
| Thm 1 | Σ_{n≤N} R² ≥ (c_1 − ε) N^{3/2}√(log N), c_1 = (3/5)^{1/2}(4/5)²/(2π√2) = 0.0558 | RH only (no pair correlation, no simplicity) | rigorous; reviewer CORRECT |
| Cor 1 | limsup |R(n)|/(n log n)^{1/4} ≥ 0.289 (RH); R(n) = Ω((n log n)^{1/4}) | none (¬RH: |R(n)| ≥ e^{cn} i.o.) | rigorous |
| Thm 2 | Σ_{N<n≤2N} R² ≪ N^{3/2}(log N)^{3/2}; ≪ N^{3/2}(log N)^{1/2} log log N with Montgomery | RH | rigorous sketch; reviewer CORRECT |
| Cor 3 | Σ_m m^δ y_m² < ∞ ⟺ δ < 1/2; RH ⟺ Σ m^δ y_m² < ∞ for some/every δ ∈ [0, 1/2) (δ = 0: Arias de Reyna 2011) | RH / equivalence | rigorous |
| Thm 3 | |R(n)| ≤ (0.202 + o(1)) √n log n/√(log log n) (0.168 with the reviewer's centring) | RH | rigorous sketch (constant fixed after review) |
| Thm 4 | Σ φ(n/N) R(n) = (7/4) Σ φ(n/N) + O(log N), φ ∈ C_c^∞(0,∞); Maślanka's "tiny oscillations" λ̃_n have mean 1 | RH | rigorous; reviewer CORRECT |
| Cor 4 | Merlini–Sala–Sala (arXiv 1903.06233 §2) "|λ_tiny(n)| = a log n, a ≈ 2": false for every a | none | rigorous + numerical counterexample; reviewer CORRECT |
| Barrier | pointwise exponent c* = inf{c: R = O(n^c)} ∈ [1/4, 1/2] under RH; lowering 1/2 needs cancellation in Σ_{γ≤√n} e^{inθ_γ} (phases ≈ n/γ ≥ √n) | — | BARRIER (von Koch analogue) |
| H1 | (1/N)Σ_{n≤N}(R − 7/4)² ~ (2√2/(3π)) √(N log N) = 0.300 √(N log N) | RH + GUE pair correlation | PREDICTION (heuristic) |
| H2 | lines of (R(u²) − 7/4)/√u at 2√(log m), amplitude π^{−1/2}Λ(m)m^{−1/2}(log m)^{−3/4}/2, phase π/4 | RH | PREDICTION, confirmed numerically |

Proof idea of Thm 1: two-sided weight φ = ψ*ψ (Φ(x) = Σ φ(n/N)e^{inx} ≥ 0), λ_{−n} = λ_n, Cauchy–Schwarz against
W(n) = Σ_{γ≤T} cos nθ_γ; all off-diagonal zero pairs then count with the favourable sign, the trend's transform is
log(2πθ)/θ² (FT of |x| log|x|), and T² log(T/2π) = (3/5) N ∫φ is optimal. Thm 2: Montgomery for the zeros below
√(N/log N), Schur + Selberg's ∫S² for the rest.

## Numerical tests (logs/)
* λ_n, n ≤ 10⁵: ball-arithmetic DFT of log ξ(1/(1−z)) (python-flint, fft_lambda.py) on |z| = 0.999 (M = 2¹⁷, 360 bits, 28 s, 117 MB)
  and 0.998 (M = 2¹⁸, 300 bits, n ≤ 5·10⁴, 77 s, 198 MB): all printed digits agree (≥ 24 decimals); r = 0.99/0.995 agree for n ≤ 10⁴.
  λ_100 = 118.603775376791, λ_500 = 991.900092992245 (= R436), λ_1000 = 2326.05316168646 (Maślanka 2326.05316169),
  λ_100000 = 462580.782406902231409416038 (= Johansson 1309.2877 §5, all 27 digits). Zero-sum check with Odlyzko's 10⁵ zeros
  (check_zeros.py): n = 100, 500, 1000, 2000 agree up to the boundary term g_n(H)S(H) (1.8e−6 … 7.3e−4), as Lemma A predicts.
* R(n): 1.378 (n = 100), 3.413 (500), 2.506 (10³), −11.815 (10⁴), −2.778 (5·10⁴), −32.421 (10⁵); max |R|/(n log n)^{1/4} = 1.558 at n = 79302 (R = 47.92).
* Thm 4 (tiny.txt): smoothed mean of R (bump on (N/2, N)) = 1.8874, 1.7704, 1.7443, 1.75115, 1.749898 at N = 10³, 3·10³, 10⁴, 3·10⁴, 10⁵
  (7/4; deviation 1.0e−4 at 10⁵ vs log N/N = 1.2e−4); λ̄_n − t(n) − 3/4 = −4.167e−4, −4.167e−5, −4.167e−7 at n = 10², 10³, 10⁵ (= −1/(24n); DFT of the Γ-part); smoothed mean of λ̃_n → 1.
  Sharp means 1.689 (R) and 0.939 (λ̃): the old reading "λ_n − trend is the mean-zero oscillating part" (R436) misses the constant.
* Mean square (analyze.txt): (1/N)Σ_{n≤N}(R − mean)²/√(N log N) = 0.041, 0.085, 0.087, 0.083, 0.100, 0.106, 0.124, 0.135 at N = 500, 10³,
  2·10³, 5·10³, 10⁴, 2·10⁴, 5·10⁴, 10⁵; dyadic blocks 0.14 → 0.184 (local exponents 0.6–0.78). H1 finite-N: data/pred = 0.76, 0.64, 0.74
  at 10³, 10⁴, 10⁵ (growth confirmed, constant ~30% high, asymptotic constant reached only at speed 1/log N).
* Thm 1: Σ_{n≤N} R²/(c_1 N^{3/2}√log N) = 1.82, 2.12, 2.03, 1.71, 1.96, 2.01, 2.30, 2.46 (N = 500 … 10⁵); the weighted two-sided version
  exceeds its bound 5.1× (N = 2·10⁴) and 4.5× (10⁵) (check_tau.txt, which also confirms τ(θ) ≈ log(2πθ)/θ² and Φ ≥ 0).
* Thm 3: max_{100≤n≤10⁵} |R(n)|/(0.202 √n log n/√(log log n)) = 0.37.
* H2 (spectral lines, analyze.txt): see derivation.md H2; m = 2: |C| = 0.1820 (pred 0.1820), phase 0.787 (π/4 = 0.785).
* Cor 4 (tiny.txt): λ̃_5083 = 17.158913 > 2 log 5083 = 17.067314 (first n ≥ 30; two radii agree to 25 digits); 7584 failures n ≤ 10⁵;
  max |λ̃_n|/log n = 4.18 (n = 79302); a = 3 fails from 46992, a = 4 from 79286.

## Literature and novelty
Checked (scouts, TeX/PDF): Arias de Reyna 2011 (Funct. Approx. 45; Thm 1/Cor 2/Thm 6: ℓ² criterion only), 2204.01036; Lagarias
math/0404394 Thm 1.1 (O(√n log n); remarks that S_f may be smaller, "remains to be explained"); Mazhouda–Voros 1506.01755 Thm 4.1;
Voros math/0506326, 1403.4558 (o(n), O(√n log n)); Coffey math-ph/0505052, 0706.0343 (per-term Laguerre cosine asymptotic (120),
|S_γ(n) + n| = O(n^{1/4}) Cor. 16, conjecture |S_Λ| = O(n^{1/2+ε})); Maślanka math/0402168 (tiny oscillations, no growth law);
Johansson 1309.2877; 0912.2390 (conjectures O(n^{1/2+ε})); Merlini–Sala–Sala 1903.06233, 1904.07005, 2012.07980; Sekatskii 1403.4484;
Odžak–Smajlović 1410.4384; 1511.04603. None has an Ω-result, a mean-square bound, the mean 7/4, or a refutation of the a log n
conjecture. arXiv listings Oct–Nov 2026 (7 cats): 0 title hits; Zenodo "Li coefficients" (69 records; 18902196, 21087821, 17957840
have no files, 22393505 scanned): nothing; GXL evidence API (3 claims): 0 support for Ω / mean square / Laguerre oscillation.
Our files: no prior entry (TREE N3 R430/R436 list only O(√n log n)). Novelty check by own search only; not adversarially reviewed.
Before → now: pointwise O(√n log n), no lower bound, no mean-square law, "oscillating part" taken as mean zero → two-sided mean
square N^{1/2+o(1)} (sharp exponent), Ω((n log n)^{1/4}), mean 7/4, explicit lines verified, one conjecture refuted.

## Review
reviewer 1 (proofs.md v1): Thm 1 CORRECT; Thm 2 CORRECT modulo Lemma A; Lemma A GAP (c_0 only numerical; Riemann–Lebesgue gives
o(1)) → repaired by the exact computation c_0 = −1/4 (Fejér kernel + θ_RS constant −π/8); Thm 3 GAP (CCM's S_1 bound is
asymmetric: −π/24 … π/48) → constant changed 0.170 → 0.202. Reviewer: "Fourier positivity and Schur arguments are standard;
novelty not established by the proofs". reviewer 2 (Lemma A repair, Thm 4, Cor 4): all three CORRECT ("∫r = −1/8; the exact Fejér sum has no constant term; c_0 = −1/4"; Thm 4 normalized error O(log N/N); Cor 4 follows from Cor 1; numerical first crossings not re-verified by the reviewer, they are by two DFT radii).

## RH relevance: N3
Rate/threshold sub-node of Li's criterion: the remainder's true mean-square order under RH is now fixed (exponent 1/2, two-sided),
giving the sharp weighted form of the Arias de Reyna equivalence (weights m^δ, δ < 1/2) and an unconditional Ω-bound; the
pointwise exponent is pinned to [1/4, 1/2] with an explicit obstruction. Positivity of λ_n is not advanced (λ_n ~ (n/2) log n dwarfs R).

## Files
proofs.md (Lemma 0, A, Thms 1–4, Cors 1–4, barrier); derivation.md (H1–H3); fft_lambda.py (λ_n and Maślanka trend by DFT;
`uv run --no-project --with python-flint python3 fft_lambda.py 0.999 131072 360 100000 logs/lam_r999.txt [gamma]`);
analyze.py, tiny.py, check_tau.py, check_zeros.py (needs Odlyzko zeros1, deleted), pred_ms.py; logs/*.txt outputs;
logs/lam_sample.txt (λ_n, 30 digits, n ≤ 1000 and every 100th n ≤ 10⁵; the full tables are regenerated in 30 s).
Peak memory (/usr/bin/time -l): fft_lambda.py 117 MB (r = 0.999) / 198 MB (r = 0.998, M = 2¹⁸), analyze.py 120 MB (328 s),
tiny.py 65 MB, check_tau.py 56 MB, check_zeros.py 54 MB, pred_ms.py 55 MB; all single-threaded. Scout evidence: .scratch/lirate/answers/.
