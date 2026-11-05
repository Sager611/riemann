NO IDEA: every transfer either dies at once or is already in use. Cluster expansions stop at Re s = 1. Lee–Yang lifts, Basu–Midha flows and the prime-factor-gas zero are RH-equivalents that give no new handle. The approximate-functional-equation main term already violates RH (1.2 % of zeros in [10³, 3·10⁴] are off-line, the first at Lehmer's pair 7005.084). The zero-attraction ODE is Hakimoglu-Brown's Λ ≤ 0.13963.

# Exploration explore1: zero-free regions in statistical mechanics → RH tree (2026-11-05)

Field: Lee–Yang/Asano, Barvinok interpolation, polymer/cluster expansions, Borcea–Brändén stability preservers, Peters–Regts complex dynamics, heat-flow/Hamilton–Jacobi zero motion.
Sources: arXiv search/listing and e-print pages (no API), the TeX of 2611.03624, our NOGO/TREE/HBDBN. Five transfers were made; three were tested numerically and two argued (D, E). None survives.

## Survey: the 10 most relevant results of 2024–2026
1. **Basu–Midha, 2611.03624 (2 Nov 2026).** ⟨0|e^{z(a_q+a_q*)}|0⟩, the Fourier transform of the q-Gaussian, has only imaginary zeros for every q ∈ [−1, 1). Its Jensen polynomials are "q-matching polynomials", proved real-rooted by a (q, z) PDE plus bounded root repulsion. Equivalently they are e^{−tD_n}zⁿ, where D_n = (4+z²)∂² − (2n+1)z∂ + n(n+2) is a Borcea–Brändén preserver. Consequences: annealed DSSYK has the Lee–Yang property, and quenched DSSYK has a zero-free strip.
2. **Gonçalves–Kerr–Klurman–Vedana, 2610.26705.** Lee–Yang representations (Kurasov–Sarnak; Alon–Cohen–Vinzant) make the zero measure of a real-rooted trigonometric polynomial crystalline. At most 2/3 of its zeros are zeta zeros, and |ζ| is large on its zeros.
3. **Wei, 2610.10652.** Fourier quasicrystals with positive integer masses are exactly the real sections of effective Lee–Yang cycles. **Wei, 2609.26802:** generators of stability-preserving semigroups have differential order ≤ 2 (a full classification), with sharp spectral gaps.
4. **Shao–Shi, 2509.06623 and 2610.12456.** A Weitz-type FPTAS on the whole Lee–Yang zero-free region, built on edge-deletion ratios. Coupling independence gives a graph-uniform zero-free neighbourhood of [0, 1] (Potts, Holant).
5. **Peters–Regts–Reppekus, 2603.17858.** Very strong spatial mixing implies no zeros near λ > 0 for the hard-core model. This is the converse direction to Regts' "zero-free ⇒ SSM".
6. **Barvinok, 2608.03687.** E exp(Σφ_i) ≠ 0 for Lipschitz multi-spin interactions once Σ L_i ≤ 1/(10√(rp)). A linear increase in interaction energy needs only a logarithmic increase in the field.
7. **Hlushchanka–Peters, 2609.35102.** Independence-polynomial zeros of recursively built graph sequences avoid a neighbourhood of [0, ∞). Proof by complex dynamics of the recursion.
8. **Bencs–Huang–Lee–Liu–Regts, 2610.40598.** W.h.p. no Fisher zeros in a strip around [0, β_KS] for the antiferromagnetic Ising model on random regular graphs.
9. **Heat flow of zeros.**
   - Höfert–Jalowy–Kabluchko 2512.17808: polynomial powers under heat flow, via a self-consistent Stieltjes-transform equation.
   - Jalowy 2609.24867: transport maps; proves three conjectures.
   - Assiotis 2608.05185: proves cases of the heat-flow conjecture.
   - Campbell–Jalowy 2605.31356: links Pólya–Schur theory to free infinitely divisible laws.
10. **Bahri–Jauslin–Lebowitz, 2607.26975.** If Lee–Yang zeros stay away from z₀, the thermodynamic limit commutes with all derivatives at z₀. **Kozitsky, 2603.18675:** the Lee–Yang property holds for isotropic D-vector models on ℤ for every even D.

Background used: Kurasov–Sarnak 2020; Alon–Cohen–Vinzant 2023 (every real-rooted exponential polynomial is a Lee–Yang restriction); Biskup–Borgs–Chayes–Kleinwaks–Kotecký 2004 (zeros at first-order transitions lie on phase-coexistence curves); Kabluchko 2022 (Curie–Weiss Lee–Yang zeros via backward heat flow); Rodgers–Tao 2018; Spira's approximate-functional-equation (AFE) conjecture, as reported by Jerby 2405.12557.

## Transfers and tests

### A (N6, N2): Basu–Midha flow for Riemann's Jensen polynomials (least dead)
- **Object.** Let X ~ Φ(|u|)du/∫Φ (the Riemann law) and R_n(z) := E(z+iX)ⁿ = Σ_ℓ (−1)^ℓ C(n,2ℓ) m_ℓ z^{n−2ℓ}.
  - These are the reversed z-variable Jensen polynomials of Ξ, in Appell form: R_n = f(∂)zⁿ with f = E e^{ixX}.
  - RH ⟺ every R_n is real-rooted (Jensen–Pólya).
- **Transfer.** Basu–Midha's mechanism needs a *closed lowering* moment recurrence ∂_p b_k = −(4a k² + 2βk) b_k + 2γ k(2k−1) b_{k−1}.
  - On the Fourier side this reads ∂_p f = −[a(x∂)² + β x∂ + γ x²] f.
  - On R_n it becomes ∂_p R_n = −[a(n−θ)² + β(n−θ) + γ∂²] R_n, with θ = z∂.
  - This is a second-order generator with leading coefficient −(a z² + γ) ≤ 0, i.e. the preserver class of Borcea–Brändén and Wei 2609.26802.
- **What would have to be proved.** A family X_p with X_0 = 0 (or Gaussian) and X_{p₀} = the Riemann law, satisfying this recurrence with a(p), γ(p) ≥ 0. All R_n(Ξ) would then be real-rooted, i.e. RH.
- **Tests** (`a_basu_midha.py`):
  1. *Constant coefficients are impossible.* For a > 0, b_k^{1/2k} stays bounded (about 2√(γ/a); a compact law). For a = 0 the law is Gaussian. The Riemann law is neither: m_ℓ^{1/2ℓ} = 0.107, 0.200, 0.264, 0.342, 0.394 at ℓ = 1, 5, 10, 20, 30 (unbounded), against 0.508 for the Gaussian of the same variance at ℓ = 30 (sub-Gaussian).
  2. *The arithmetic family, the de Bruijn–Newman (dBN) flow, raises instead of lowering.*
     - Its recurrence is ∂_t b_k = b_{k+1}, so ∂_t R_n = −(R_{n+2} − 2zR_{n+1} + z²R_n), which is not closed in degree n.
     - It is provably not a preserver. As t → −∞, R_n(t)/σⁿ → He_n(z/σ) with σ² = 1/(2|t|), which is real-rooted.
     - Yet for every t < Λ some R_n(t) has non-real roots: Jensen's theorem plus Newman's Λ > −∞ (Rodgers–Tao: Λ ≥ 0).
  3. *The Ornstein–Uhlenbeck flow lowers and closes but points the wrong way.* ∂_p R_n = −(∂² − z∂ + n)R_n has Hermite as its attractor, so the preserving direction runs Riemann → Gaussian.
  4. *Time-dependent coefficients are circular.* Reaching the Riemann law (unbounded support, Fourier transform of order 1, maximal type) needs a(p) → 0 at p₀ in a tuned way. That is the Pólya–de Bruijn kernel program (NOGO R7.6, R34).
  - Sanity check: R_n is real-rooted for n ≤ 30 (60 digits).
- **Verdict:** dies.

### B (N8, N7): Lee–Yang lift of the Riemann–Siegel main term
- **Transfer.** By Alon–Cohen–Vinzant, an exponential polynomial is real-rooted exactly when it is the restriction of a Lee–Yang polynomial; GKKV used this for ζ.
  - Locally Z = Z_N + R, with Z_N = 2Σ_{n≤N} n^{−1/2} cos(θ − t log n) and N = ⌊√(t/2π)⌋. After linearising θ, Z_N is an exponential polynomial.
  - The two "phases" ζ_N(s) and χ(s)ζ_N(1−s) coexist exactly on Re s = 1/2 (the BBCKK picture).
- **What would have to be proved.**
  1. F_N(s) = ζ_N(s) + χ(s)ζ_N(1−s) has only critical-line zeros for 2πN² ≤ t < 2π(N+1)². This is equivalent to the local model having a Lee–Yang lift.
  2. A Rouché bound for the remainder R.
- **Test** (`b_lehmer.py`, `b_scan.py`). Item 1 is false:
  - At Lehmer's pair, Z(7005.08) = +0.0039 but Z_33 = −0.0660 (the remainder is ≈ +0.070), so Z_33 has no real zeros there.
  - F_33 has zeros at 1/2 ± 0.07648 + 7005.08414i (|F| < 2·10⁻²⁷). The off-line distance matches √(2|R|/|Z''|) = 0.078.
  - Scan of [10³, 3·10⁴]: Z_N has 34620 sign changes against 35024 zeros of ζ (equal to N(3·10⁴) − N(10³)). There are 280 lost close pairs; 24 of the first 40 have a verified off-line zero of F_N, with σ − 1/2 = 0.060–0.202.
- **Known in kind.** Spira: the classical AFE approximation fails RH (Jerby 2405.12557).
- **Consequence.** Every fixed-order Riemann–Siegel model fails at pairs whose local max |Z| is below its remainder, and such pairs are expected at every order (GUE small gaps). A two-phase / Lee–Yang argument must use ζ itself.
- **Verdict:** dies at once.

### C (N10): Lee–Yang zero of the squarefree prime-factor gas
- **Object.** P_x(z) = Σ_{n≤x} μ²(n) z^{ω(n)}: a hard-core "knapsack" gas of primes with Π p ≤ x and fugacity z. P_x(1) = Q(x) (the count of squarefree n ≤ x) and P_x(−1) = M(x).
- **Structure.** Selberg–Delange gives P_x = x(log x)^{z−1}[G(z)/Γ(z) + O(1/log x)], with G(z) = Π_p (1 + z/p)(1 − 1/p)^z.
  - Every term vanishes at z = −k, so zeros of P_x converge to −1, −2, …, and P_x'(−1) ≍ x/log²x.
  - Hence RH ⟺ the real zero r(x) near −1 satisfies |r(x) + 1| ≪ x^{−1/2+ε}: a "Lee–Yang edge" with finite-size exponent 1/2.
- **Tests** (`c_pgf.py`):
  - r(x) + 1 tracks −M(x)/P_x'(−1): at x = 10⁷, 0.01906 against 0.02003, with M = 1037.
  - Real-rootedness would give a Poisson-binomial law for ω(n) and M(x)/Q(x) = Π(1 − 2p_i), but it fails:
    - first at x = 114 (9z³ + 32z² + 30z + 1, discriminant −28139);
    - for 42.5 % of squarefree x ≤ 10⁶;
    - with a complex pair near −1.59 at x = 10⁷;
    - and Σ_{n≤x} z^{ω(n)} fails from x = 42732.
- **Mismatch.** Statistical-mechanics tools give zero-*free* regions; here RH is a zero-*attraction* rate.
- **Verdict:** dies; a reformulation only.

### D (N7): polymer / cluster expansion of log ζ over primes (argued)
- Primes are non-interacting polymers (the Euler product), so the Kotecký–Preiss condition reduces to Σ_p |p^{−s}| < ∞ ⟺ σ > 1.
- In the strip, ζ is an analytic continuation through the pole at 1, so no small-activity polymer representation exists there.
- The mollified version ζM − 1 is the zero-density detector N7 already uses.
- **Verdict:** dies at once (R48).

### E (N2): zero dynamics / Hamilton–Jacobi attraction (argued; known)
- **Inequality.** For H_t, the zero with the largest imaginary part Y obeys dY/dt ≤ −1/Y − Σ_j 2Y/((x−x_j)² + 4Y²). Every other zero attracts it, real or not, because for a conjugate pair min_A [f(A) + f(2Y−A)] ≥ f(2Y) with f(u) = u/(c+u²).
- **Lower bound without S(T).** The sum is bounded below through −Im(H'/H) evaluated above the zeros: ¼ log(x/4π) + ½ Re ζ'/ζ((1+η)/2 + ix/2).
- **Heuristic.** With real neighbours, starting from Y = 3/4 (the 7/8 theorem) at x ≈ 6·10¹², Y reaches 0 by t ≈ 0.06.
- **Already done.** This is Hakimoglu-Brown's mechanism: log-derivative at height 61/20, a conjugate-pair inequality and a 24-step certificate giving Λ ≤ 0.13963 (audited in [HBDBN](../HBDBN/REPORT.md)).
- **Verdict:** known; any further gain is N2-internal ([DBN6](../DBN6/REPORT.md)).

## Ranking (RH relevance × chance × gain)
| rank | idea | node | relevance | chance | gain | result |
|---|---|---|---|---|---|---|
| 1 | A: Basu–Midha lowering flow to the Riemann law | N6 (N2) | high | ≈ 0: constant coefficients excluded; general case = RH | RH | dies (test 1) / circular (test 4) |
| 2 | E: zero-attraction ODE for Λ | N2 | medium | done by Hakimoglu-Brown | record exists | known |
| 3 | B: Lee–Yang lift of the AFE main term | N8, N7 | medium | 0 | — | false at Lehmer's pair; 1.2 % of zeros off-line |
| 4 | C: Lee–Yang zero of the prime-factor gas | N10 | medium | ≈ 0 | — | not real-rooted from x = 114; RH-equivalent only |
| 5 | D: cluster expansion over primes | N7 | high | 0 | — | absolute convergence only for σ > 1 |

Live lead for the N6 owners: Basu–Midha's q-matching method is the nearest working real-rootedness technique. The exact obstruction is the need for an arithmetic deformation whose moments obey a *lowering* recurrence starting from z^n.

## Suggested NOGO line (orchestrator)
- explore1, stat-mech zero-freeness: cluster expansions stop at Re s = 1. Lee–Yang lifts, Basu–Midha flows (the dBN Jensen flow raises; Ornstein–Uhlenbeck points to Hermite) and the prime-factor-gas zero (P_x not real-rooted from x = 114) are RH-equivalents. The AFE main term fails RH (F₃₃ zero at 1/2 + 0.0765 + 7005.084i; 1.2 % of zeros in [10³, 3·10⁴]). The zero-attraction ODE is Hakimoglu-Brown's. [report](Riemann/Deep/ExploreStatMech/REPORT.md)

## Files
`a_basu_midha.py` (53 s), `b_lehmer.py` (17 s), `b_scan.py` (1000 s, dominated by mpmath findroot), `c_pgf.py` (29 s, about 0.3 GB); raw numbers in `outputs.txt`. All runs were single-threaded at nice 19 under a hard alarm.
