# Fasi–Higham positivity: the diagonal case, m = 3, and where the conjecture first fails

**Notation.** q_{k,m}(z) = Σ_{j=0}^{m} (−1)^j C(m,j) (k+m−j)!/(k+m)! z^j is the denominator of the
[k/m] Padé approximant of e^z, p_{k,m} its numerator, 1/q_{k,m}(z) = Σ c_n(k,m) z^n, M = k+m,
σ = k/M, τ = n/M. Fasi–Higham (SIMAX 40 (2019), Conj. 3.2): c_n(k,m) > 0 for all n when m is odd.
Scardino (arXiv 2610.39051, read in TeX) refutes it for k ≫ m (c₃₈(66,5), c₃₇(44,7), c₄₂(45,9) < 0;
infinitely many for m = 5), proves k = 0 and m = 1, eventual positivity under D1 (the positive zero R
is the unique zero of least modulus) and D1 for fixed m, large k. Open there: the diagonal k = m
(the case used in expm's error bound) and every fixed m ≥ 3 with all k; m = 3 is delicate because its
Hermite profile Ψ₃(λ) = 1 − cos(√3 λ) touches 0.

## 1. Results

1. **Theorem (m = 3).** c_n(k,3) > 0 for every k ≥ 0 and n ≥ 0; in fact
   c_n(k,3) ≥ (R₀ − 2|R₊|) z₀^{−n} with 1 − 2|R₊|/R₀ = 1/(3(k+1)) + O(k^{−2}). Proof in §2.
2. **Diagonal, all n (rigorous, computer-assisted).** c_n(m,m) > 0 for all n ≥ 0 for every odd
   m ≤ 201 (Arb balls + certified roots; `arb_cert.py diag 201`). `verify.py` re-proves odd m ≤ 101
   with mpmath only (exact integer signs for n < n₁, interval-arithmetic pole bound for n ≥ n₁).
   The previous check (m ≤ 61, n ≤ 12m) did not reach the tail threshold n₁ ≈ 0.73 m² (n₁ = 2685 at
   m = 61, 7478 at m = 101, 29970 at m = 201).
3. **Whole columns below the first failure, all n.** For every odd m ≤ 101, every k < k_min(m)
   gives c_n(k,m) > 0 for all n (`arb_cert.py col`), in particular every k ≤ 2m + 11 and so the
   whole near-diagonal (also run separately: k = m ± 1, ±2, ±4 for odd 59 ≤ m ≤ 101).
4. **First failures** (table §5). k_min(m) := least k with some c_n(k,m) < 0, exact for every odd
   m ≤ 101 (negative certified exactly, every smaller k positive for all n by 3).
   k_min(m) = 2m + c with c = 18 (m = 13, 15), 16 (17, 19), 15 (21–33), 13 (35–71), 14 (73–93),
   12 (95–101); the first negative is at n = M − d with
   d = 6, 5, 4, 3, 2 on the same kind of ranges (d = 3 for 35 ≤ m ≤ 93), just below the Padé order M+1.
   So k_min(m) − m = m + c > 0 and grows: it never changes sign; failures recede from the diagonal.
5. **Asymptotic picture (heuristic, matches all data).** Failures occur for large M iff
   σ > σ*(τ), with σ*(τ) = 0.900, 0.863, 0.816, 0.758, 0.720, 0.692, 0.676 at
   τ = 0.6, 0.7, 0.8, 0.9, 0.95, 0.98, 0.994 and σ*(τ) → 2/3 as τ → 1⁻. Hence k_min(m) = 2m + O(1)
   with the first failure at n = M − O(1). The diagonal (σ = 1/2) is outside with exponential room
   for n ≤ 2m and an O(1) margin beyond (§3). No counterexample is expected on the diagonal.

## 2. Proof of the m = 3 theorem

Put N = k+1, y = N − z. Then (k+1)(k+2)(k+3) q_{k,3}(z) = Q := y³ + 3Ny + 2N. The cubic has one real
root y₀ = −d and two complex roots y_± = d/2 ± i(√3/2)s, where d is the root of
F(d) = d³ + 3Nd − 2N in (0, 2/3) and s² = d² + 4N (Vieta: Re y_± = d/2, |y_±|² = 2N/d = d² + 3N).
F is increasing, F(1/2) = 1/8 − N/2 < 0 and F(2/3) = 8/27 > 0, so 1/2 < d < 2/3 for N ≥ 1.
The zeros of q are z₀ = N + d > 0 and z_± = N − d/2 ∓ i(√3/2)s.

*D1.* |z_±|² − z₀² = (N − d/2)² + (3/4)(d² + 4N) − (N + d)² = 3N(1 − d) > 0.

*Residues.* c_n = Σ_j R_j z_j^{−n} with R_j = −1/(z_j q'(z_j)) = Q(0)/(3 z_j (y_j² + N)), and the cubic
gives y_j² + N = −2N(y_j + 1)/y_j, so R_j = −Q(0) y_j / (6N z_j (y_j + 1)). Thus
R₀ = Q(0) d /(6N(N+d)(1−d)) > 0 and, using |y₊|² = 2N/d, |y₊ + 1|² = d² + d + 1 + 3N,
|z₊|² = (N+d)² + 3N(1−d):

  (2|R₊|/R₀)² = 8N(N+d)²(1−d)² / ( d³ [(N+d)² + 3N(1−d)] [3N + d² + d + 1] ).

With N = d³/(2 − 3d) (equivalent to F(d) = 0) one has the polynomial identity (checked in verify.py)

  (2−3d)³ ( d³[(N+d)² + 3N(1−d)][3N+d²+d+1] − 8N(N+d)²(1−d)² ) = d⁵ (1−d)² (3d−2) g(d),
  g(d) = −3d³ + 12d² − 30d + 12.

g'(d) = −9d² + 24d − 30 < 0 everywhere and g(1/2) = −3/8, so g < 0 on [1/2, ∞); with 3d − 2 < 0 the
right side is positive, i.e. 2|R₊| < R₀. Therefore
c_n = R₀ z₀^{−n} + 2 Re(R₊ z₊^{−n}) ≥ R₀ z₀^{−n} − 2|R₊| |z₊|^{−n} ≥ (R₀ − 2|R₊|) z₀^{−n} > 0. ∎
(The argument works for every real N > 1/4. As N → ∞, d = 2/3 − 8/(81N) + …, 2|R₊|/R₀ → 1:
the inequality is tight to first order, which is why Ψ₃ touches 0.)

## 3. Structure behind the diagonal and the failures

* **Padé identity.** p/q = e^z + O(z^{M+1}) gives c_n(k,m) = [z^n] e^z / p_{k,m}(z) for n ≤ M.
  For odd m the exact expansion is 1/q = Σ_{ℓ≥0} e^z z^{(M+1)ℓ} E(z)^ℓ / p(z)^{ℓ+1},
  E(z) = (1/M!) ∫₀¹ (1−v)^k v^m e^{zv} dv (positive coefficients), so
  c_n = Σ_{ℓ ≤ n/(M+1)} T_ℓ(n), T_ℓ(n) = [z^{n−(M+1)ℓ}] e^z E^ℓ / p^{ℓ+1}. On the diagonal p(z) = q(−z).
  (For even m the signs alternate in ℓ: this is where oddness of m enters.)
* **Partial fractions mislead.** On the diagonal the residues are huge and cancel: R₁ ≈ 2.3·10¹⁸ at
  m = 61 while c₀ = 1. The ratio c_n R^n / R₁ is a theta-function-like profile with deep positive
  minima at n ≈ 2.16m, 4.46m, 7.08m (≈ 10⁻¹⁴, 4·10⁻⁸, 2·10⁻⁵ at m = 61). The "first complex pair"
  (2|w₁|/w_r → 2) therefore says nothing about signs; the relevant objects are saddle points.
* **Saddles.** With w = z/M, q(Mw) ∝ ∫₀^∞ e^{M f(u,w)} du, f = σ log u + (1−σ) log(u−w) − u,
  u_± the roots of u² − (1+w)u + σw = 0. The main term is the positive saddle
  w* = τ(1+τ)/(1+τ−σ) (u₊ = 1+τ): positive, size e^{−M[f₊(w*) + τ log w*]}. The oscillating competitor
  is the first Poisson-dual term of the pole sum of e^z/p along p's zero arc; its saddles solve
  2w² + (5 + τ − 9σ) w + (2+τ)(1−τ) = 0, exponent gap
  Δ(σ,τ) = Re[f₋(w₁) − 2f₊(w₁) − τ log w₁] + f₊(w*) + τ log w*. Prediction: failure iff Δ > 0,
  giving σ*(τ) above; at (σ, τ) = (2/3, 1) the two dual saddles merge at w = 0. Checks: first failure
  (95,41) at (σ,τ) = (0.699, 0.978) vs σ*(0.978) ≈ 0.692; (214,101) at (0.679, 0.994) vs 0.677;
  negatives of (150,50) start at τ = 0.92 vs τ*(0.75) = 0.911; σ = 0.6 (k = 61, m = 41): none.
* **Diagonal margin.** Δ(1/2, τ) < −0.17 for all 0 < τ ≤ 1 and ≤ −0.27 for τ ≥ 0.3 (−0.22, −0.65,
  −0.71, −0.52, −0.31 at τ = 0.01, 0.3, 0.5, 0.9, 0.99): for n ≤ 2m the oscillating part is below
  e^{−0.34 m} of the main term.
  Beyond n = 2m the closest approach to zero is the dip at n ≈ 2.17m where T₁ takes over. With
  A_n = (c_n − √(c_{n−1}c_{n+1}))/(c_n + √(c_{n−1}c_{n+1})) (A_n ≤ −1 iff c_n ≤ 0 between positive
  neighbours), D(m) = 1 + min_{n ≤ 8m} A_n = 0.83 (m = 5), 0.76 (21), 0.75 (41), 0.73 (89–101),
  0.75 (149), 0.76 (197): bounded away from 0, no trend. From n₁ ≈ 0.73 m² the D1 pole bound holds.

## 4. What is still missing for all odd m (diagonal)

A proof needs uniform-in-m versions of §3: (i) n ≤ 2m: rigorous Stokes geometry of p_{m,m}'s zero
arc with explicit error terms (the exponential room Δ < −0.17 is there); (ii) the crossover windows
n ≈ 2ℓm + 0.17m (only an O(1) margin, D ≈ 0.73); (iii) 4m ≲ n ≲ 0.73m², a theta-function regime;
(iv) an explicit D1 tail (R_j/R₁ and the gap |z₁|/R − 1 ≈ 1.04/m² as functions of m). Uniform
Airy/Bessel expansions of reverse Bessel polynomials make (i)–(iv) plausible but long. No simple
positivity mechanism (total positivity, the θ_m continued fraction, h(t) ≥ 0 for the Bessel-filter
impulse response: false, it rings) was found.

## 5. First failures k_min(m) (first negative index n in brackets)

m=5: 66 (38) · 7: 44 (37) · 9: 42 (41) · 11: 43 (46) · 13: 44 (51) · 15: 48 (57) · 17: 50 (62) ·
19: 54 (68) · 21: 57 (74) · 23: 61 (80) · 25: 65 (86) · 27: 69 (92) · 29: 73 (98) · 31: 77 (104) ·
33: 81 (110) · 35: 83 (115) · 37: 87 (121) · 39: 91 (127) · 41: 95 (133) · 43: 99 (139) ·
45: 103 (145) · 47: 107 (151) · 49: 111 (157) · 51: 115 (163) · 53: 119 (169) · 55: 123 (175) ·
57: 127 (181) · 59: 131 (187) · 61: 135 (193) · 63: 139 (199) · 65: 143 (205) · 67: 147 (211) ·
69: 151 (217) · 71: 155 (223) · 73: 160 (230) · 75: 164 (236) · 77: 168 (242) · 79: 172 (248) ·
81: 176 (254) · 83: 180 (260) · 85: 184 (266) · 87: 188 (272) · 89: 192 (278) · 91: 196 (284) ·
93: 200 (290) · 95: 202 (295) · 97: 206 (301) · 99: 210 (307) · 101: 214 (313).
m = 3: no failure (theorem). Note k_min(9) = 42 < 45 of the paper's example.

## 6. Methods and files

* **Exact signs** (`verify.py`): Y_n = round(Σ_j −b_j s^j Y_{n−j} / b₀), Y₀ = 2^P, b_j the integer
  coefficients of b₀ q. The rounding errors ε_i (|ε_i| ≤ 1/2) propagate through the recurrence's own
  impulse response y_n = c_n s^n, so |Y_n − 2^P y_n| ≤ ½ Σ_{i<n} |y_i|, bounded on the fly; the sign
  of c_n is certain when |Y_n| exceeds the bound (P doubles otherwise). No cancellation problem.
* **Tail** (both scripts): certified root enclosures (Arb `acb_poly.roots`; in verify.py Aberth roots
  plus Gerschgorin discs for A = diag(z_i) − W 1ᵀ, whose characteristic polynomial is Q/b_m),
  interval residues R_j = −b₀/(z_j Q'(z_j)), D1 check, and Σ_{j≥2} |R_j/R₁| (R/|z_j|)^{n₁} < 1,
  which gives c_n ≥ (R₁ − Σ_j |R_j| (R/|z_j|)^n) R^{−n} > 0 for all n ≥ n₁.
* **Files:** `verify.py` (stdlib + mpmath; `python3 verify.py` re-checks §2, every counterexample of
  §5 exactly, c_n(k_min−1, m) > 0 for n ≤ 4M+200, and the diagonal for all n, m ≤ 101; ≈ 3 min);
  `arb_cert.py` (python-flint; reruns 2–4: diagonal to 201 ≈ 10 min, columns to m = 101 ≈ 2 h).
