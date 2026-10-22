/-!
# Round 136 — independent checks of A > 1 for CM fields; the Smale mean-value claim; the Gauss–Fourier gap

* CM fields, every A > 1 (Riemann/Deep/CM3): reviewer CORRECT for the degree-one theorem; the sharpness
  sentence needed "degree-one" (in Q(√−163), (2) is a principal prime of norm 4 < 163/4), fixed. Scouts (5/5 OK):
  - checks/lemmaP_rel.gp rerun: 141/141, 372 blocks, 12073 vanishing pairs; it tests the identity as stated;
  - poisson_e2: brute/main 0.9415–1.0161 (finite-size caveat at the smallest case, as stated);
  - Theorem 1: exponents check, no N³/X term, stated for detector sequences only (false for arbitrary b);
  - Theorems 2–3: Markov margin ≥ 1.125e−5 on a grid A ∈ (1, 1.2], σ ∈ (1/2, 9/10]; exceptional set O(X^(1−c));
  - sharpness: F = Q(√5), K = F(√−1), F(√−3), F(√−7): least principal degree-one prime norms 5, 19, 11 vs
    bounds 1, 9/16, 49/16.
* Smale's mean value constant (arXiv 2610.21472): 𝓜 := sup over P(0) = 0, P′(0) = 1, deg P ≥ 2 of
  min over critical points c of |P(c)/c| — exactly Smale's normalized constant. The paper proves 𝓜 = 𝓑 = 𝓐
  (𝓑 = sup 1/B(0) over admissible univalent Blaschke pairs) and claims "a new disproof of Smale's mean value
  conjecture" with 𝓜 > 1.1018 from the pair (9/16, 3/8). No explicit polynomial: existence via Runge–Hilbert
  approximation; the auxiliary normalized quartic z − (3/2)z² + z³ − z⁴/4 has S = 1/4. A knowledge scout
  claimed the bound is compatible with the conjecture; that contradicts the TeX and is discarded.
  → Orca: audit 𝓜 ≥ 𝓑 and build an explicit certified counterexample.
* Gauss–Fourier conjecture (arXiv 2610.29790): real gap. Γ is defined circularly (l. 1238–1241), "directly
  intersect" is informal (l. 1248–1264), Theorems 3.22/3.24 (l. 1622–1670) are conditional, and l. 1691–1694
  applies them to arbitrary loci. No numerical counterexample (degree 5–7 search timed out). → reviewer.
-/

namespace Fayed.R136
end Fayed.R136
