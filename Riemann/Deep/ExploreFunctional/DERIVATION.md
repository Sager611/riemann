# explore39 derivation: Beurling–Deny structure of the windowed Weil form and what functional inequalities can certify

Conventions (../WeilZeta/DERIVATION.md, ../ExploreEntropy/REPORT.md): I = [−a, a], a = ½ ln c; Q(f) = pole + Arch − primes;
pole(f) = 2⟨f, e^{−x/2}⟩⟨f, e^{x/2}⟩ = 2⟨f, C⟩² − 2⟨f, s⟩², C = cosh(x/2), s = sinh(x/2); primes(f) = 2 Σ_{n<c} w_n C_f(ln n),
C_f(h) = ∫ f(x) f(x−h) dx, w_n = Λ(n)/√n; Arch(f) = D‖f‖² + ½∬_{ℝ²} κ(|x−y|)(f(x) − f(y))² dx dy (f = 0 off I; explore32 S0),
κ(z) = e^{z/2}/(e^z − e^{−z}) = Σ_m e^{−(2m+½)z}, D = ψ(1/4) − ln π = −5.3721834192,
K(r) := ∫_r^∞ κ = arctan(e^{−r/2}) + ½ ln coth(r/4) (differentiate: −1/(4cosh(r/2)) − 1/(4 sinh(r/2)) = −κ(r)).
J(z) := κ(|z|) − 2cosh(z/2) = Σ_{k≥1} e^{−(2k+½)z} − e^{z/2} (z > 0): decreasing on (0, ∞), and convex on (0, (ln 25)/3) ⊃ (0, 2 ln ρ]
since J'' ≥ (25e^{−5z/2} − e^{z/2})/4 > 0 there (reviewer's certificate; numerically min J'' = 5.05 on (0, 2 ln ρ], t9); J(ln ρ) = 0, ρ³ = ρ + 1.
"Markov" below means positivity preserving in the Beurling–Deny sense: a Dirichlet form (nonnegative jump kernel) plus a potential bounded
below, not necessarily sub-Markovian (D < 0). Kernel integrals ∬ j f f with the log-singular κ are formal; the form is ∫U f² + ½∬ j (f(x) − f(y))²,
and only bounded pieces (j⁻, atoms) are ever integrated separately.
Davenport–Heilbronn (DH): κ_DH(z) = e^{−z/2}/(e^z − e^{−z}), D_DH = ln(5/π) + ψ(3/4), no pole, atoms b(n)/√n with b(2) = 0.197,
b(3) = −0.312, b(4) = −1.442 (wmat.dh_coeffs; explore17 T2).

## 1. Markov part (needs Λ ≥ 0)
−2w C_f(h) = w‖f − τ_h f‖² − 2w‖f‖² (on ℝ) gives Q = E_c − S(c) + pole, with
  E_c(f) = ½∬_{ℝ²} j_M(x − y)(f(x) − f(y))², j_M(z) = κ(|z|) + Σ_{n<c} w_n (δ(z − ln n) + δ(z + ln n)) ≥ 0,   S(c) = 2Σ_{n<c} w_n − D.
E_c is the Dirichlet form of the symmetric Lévy process with Lévy measure κ(|z|)dz + Σ_{n<c} w_n(δ_{ln n} + δ_{−ln n}), killed on leaving I
(characteristic exponent φ(t) − D + 2Σ w_n(1 − cos(t ln n)), φ = Weil's archimedean symbol). S(2) = 5.3722, S(3) = 6.3524, S(5) = 8.3142,
S(10) = 12.4472, S(20) = 18.3878, S(47) = 26.9856. Only the pole is anti-Markov on the full space (explore17 T2).

## 2. Parity reduction and sign structure (Theorem A)
Write Q(f) = ∫ V f² − ∬ j(x − y) f(x) f(y), j = J + atoms w_n δ(· ∓ ln n). For f even/odd with g = f|[0,a] (‖f‖² = 2‖g‖²):
  ∬_{I²} j(x−y) f f = 2∬_{[0,a]²} [j(x − y) ± j(x + y)] g g,   so   j_ev = j(x−y) + j(x+y),  j_odd = j(x−y) − j(x+y).
(a) Even, continuous part: j_ev = J(|x−y|) + J(x+y) ≥ min over s = x + y of J(min(s, 2a − s)) + J(s) ≥ 2J(a) (J decreasing, then convex
    with midpoint a), attained at (a, 0). So j_ev ≥ 0 ⟺ J(a) ≥ 0 ⟺ a ≤ ln ρ ⟺ c ≤ ρ² = 1.754877666 (explore32's ρ = 1.3247 is the
    full-space threshold, where J(2a) changes sign; the even sector stays Markov to ρ²). Even atoms: +w_n at x − y = ±ln n and at
    x + y = ln n, Markov iff w_n ≥ 0 (true for ζ). Numerics t1: min j_ev = 3.60 (c = ρ), 2.42 (1.4), −4.1·10⁻⁴ (1.755).
(b) Odd, continuous part: J(|x−y|) − J(x+y) = [κ(|x−y|) − κ(x+y)] + 4 sinh(x/2) sinh(y/2) ≥ 0 for every c (the pole's −2⟨f,s⟩² is
    Markov-compatible on odd functions). Odd atoms: +w_n at x − y = ±ln n (n < √c), −w_n at x + y = ln n (n < c): the reflected prime
    atoms are anti-Markov. So the odd sector is Markov ⟺ no n < c ⟺ c ≤ 2.
(c) DH: j_ev = κ_DH(|x−y|) + κ_DH(x+y) > 0, j_odd = κ_DH(|x−y|) − κ_DH(x+y) > 0; atoms b(2) > 0, b(3), b(4) < 0: even Markov ⟺ c ≤ 3,
    odd ⟺ c ≤ 2.
Consequence (Perron–Frobenius per sector, Reed–Simon XIII.44): on the folded space L²([0,a]) the sector semigroup is positivity improving
(j_sec > 0 off a null set; irreducible), so the sector ground state is simple and positive a.e.: even for c ≤ ρ², odd for c ≤ 2.

## 3. Barta / Allegretto–Piepenbrink certificate
For a Markov sector form (j_sec ≥ 0) and φ > 0: with f = φ g,
  Q(f) = ∫ (Hφ/φ) f² + ½∬ j_sec φ(x)φ(y)(g(x) − g(y))²  ⟹  λ₀ ≥ inf_x (Hφ)(x)/φ(x),
(Hφ)(x) = (D + K(a − x) + K(a + x)) φ(x) + ∫_I κ(|x − y|)(φ(x) − φ(y)) dy + e^{−x/2}⟨φ, e^{x/2}⟩ + e^{x/2}⟨φ, e^{−x/2}⟩ − Σ_n w_n[φ(x − ln n) + φ(x + ln n)]
(φ = 0 off I). Admissible φ: polynomials on I, positive on [−a, a] (odd sector: on (0, a] with φ′(0) > 0), extended by 0; they lie in the form
domain (the jumps at ±a have finite log-energy), Hφ is finite on the open interval and ≈ ½ ln(1/(a − |x|)) φ at the ends, and the identity
holds on the folded sector space. Equality for the true ground state, so sup_φ inf_x (Hφ/φ) = λ₀ (two-sided, exact). Code: xform.Hmat, barta_lp (LP over the
Legendre span: maximise t subject to Hφ − tφ ≥ 0, φ ≥ 0 on a 400-point grid; bisection).

## 4. Flip minorant
Split j_sec = j⁺ − j⁻ (j⁻ ≥ 0 the anti-Markov part). Since ∬ j⁻ f f ≥ −∬ j⁻ |f||f| and ∬ j⁺ f f ≤ ∬ j⁺ |f||f|:
  Q(f) ≥ Q_flip(|f|),   Q_flip := Q − 2∬ j⁻ f f  (Markov),   hence λ₀(Q) ≥ λ₀(Q_flip) ≥ inf_x (H_flip φ/φ).
Q_flip changes only the anti-Markov part of the kernel, by the least pointwise amount that makes it Markov: it is the canonical
term-wise Beurling–Deny comparison of the explicit-formula terms (a finer minorant must use the sign pattern of f).
Anti-Markov parts flipped: ζ even: (J(x−y) + J(x+y))⁻ (the pole's excess), ζ odd: reflected atoms 2w_n ∫ g(x) g(ln n − x);
DH: negative-coefficient atoms (even: direct and reflected; odd: direct), positive reflected atoms (odd).

## 5. Even-sector certificate: Barta condition + Poincaré inequality (Theorem B)
On even functions Q = M + 2CC^T with M = E_c − S(c) Markov on the full space for every c (§1). For even ψ > 0 and f = ψ g the ground-state
identity for M (valid for any symmetric kernel; sign used only for ≥) and Mψ = Qψ − 2C⟨C, ψ⟩ give exactly
  Q(f) − ∫ (Qψ/ψ) f² = E_ψ(g) − R_ψ Var_ν(g),   E_ψ(g) = ½∬ j_M(x−y) ψ(x)ψ(y)(g(x) − g(y))²,   R_ψ = 2⟨C, ψ⟩²,   ν = Cψ dx/⟨C, ψ⟩.
Hence: Q^ev ≥ λ ⟸ ∃ ψ > 0 even with (i) Qψ ≥ λψ pointwise and (ii) the Poincaré inequality E_ψ(g) ≥ R_ψ Var_ν(g) ∀g.
Conversely, if the even ground state ψ₀ is positive (c ≤ ρ² by §2; numerically every c, explore32), (i) holds with λ = λ₀ and (ii) holds with
Poincaré constant P = R + μ_W, μ_W = inf_{f ≁ ψ₀} [Q(f) − λ₀‖f‖²]/W(f), W(f) = (1/⟨C,ψ₀⟩)[∫ (C/ψ₀) f² − ⟨C, f⟩²/⟨C, ψ₀⟩] = Var_ν(g).
(infimum over even f with W(f) > 0). Testing f = φ₁ (second even eigenvector, ‖φ₁‖ = 1): μ_W ≤ (λ₁ − λ₀)/W(φ₁). At ψ₀ condition (ii) is automatic (Q − λ₀ ≥ 0) and positivity is
condition (i), λ₀ ≥ 0 (reviewer, round 2). So the Poincaré route needs a test function with Qψ ≥ 0 pointwise, i.e. inf Qψ/ψ ∈ [0, λ₀] (whenever
(ii) holds, inf Qψ/ψ ≤ λ₀), and a proof of (ii) for the ψ-Doob transform of the prime-jump Lévy process against the rank-one pole; the relative
margin of (ii), even at ψ₀, is (P − R)/R, controlled by the second even Weil eigenvalue.

## 5′. Odd sector beyond c = 2 (comparison of two jump forms)
Folded odd form q(g) = ∫ U_odd g² + ½∬ k_odd (g(x) − g(y))², k_odd = k⁺ − k⁻, k⁻ = Σ_{n<c} w_n δ(x + y − ln n) (reflected atoms; the pole and the
continuous part are in k⁺). With ψ > 0 the odd ground state and g = ψh: q(g) − λ₀‖g‖² = E⁺_ψ(h) − E⁻_ψ(h),
E⁻_ψ(h) = ½ Σ_n w_n ∫ ψ(x)ψ(ln n − x)(h(x) − h(ln n − x))² dx, each integral over x, ln n − x ∈ [0, a]. At ψ the domination E⁺_ψ ≥ E⁻_ψ is
automatic and odd positivity is λ₀ ≥ 0 (reviewer, round 2: subtracting a constant leaves both energies unchanged). For a test function φ ≈ ψ
the generalized Barta certificate needs Hφ ≥ 0 pointwise and E⁺_φ ≥ E⁻_φ (the reflected prime jumps dominated by the Markov ones); the latter's
relative margin at φ = ψ is η = inf [q − λ₀‖·‖²]/E⁻_ψ (t11): 2.65, 1.82, 1.06, 0.29, 1.5·10⁻², 5.0·10⁻⁴, 8.0·10⁻⁶ at c = 2.05, 2.1, 2.2, 2.5,
3, 3.5, 4, about (λ₁ − λ₀)^{odd}/5: looser than the even sector at equal c (explore32: the even sector fails first).

## 6. Pinning (interlacing)
Even: Q^ev = M^ev + 2CC^T ⟹ μ_k ≤ λ_k ≤ μ_{k+1} (μ_k eigenvalues of M^ev, λ_k of Q^ev). With Q ≥ 0: λ₀ ≤ μ₁ ≤ λ₁, i.e.
  λ₁(E_c^ev) − S(c) ∈ [λ₀(Q^ev), λ₁(Q^ev)],   and generally λ_k(E^ev) − S ∈ [λ_{k−1}, λ_k] (Q^ev).
μ₀ < 0 is not automatic (reviewer: as c ↓ 1 the exterior killing makes M^ev = Arch^ev positive). It holds for c ≥ 1.24 given a negative
even Rayleigh quotient of Arch there (explore32's Galerkin crossing a = 0.10717, c = 1.2390, and ours; float64, not validated): the bottom
decreases in c by domain monotonicity, and for c > 2 the Perron vector ψ_A > 0 of Arch^ev has primes(ψ_A) = 2Σ w_n C_{ψ_A}(ln n) ≥ 0, so
μ₀(M^ev) ≤ Arch(ψ_A) < 0. At the windows used μ₀ = −1.285, −3.334, −5.119, −7.216 (c = 2, 5, 10, 20), far from 0; there
λ₀(E_c^ev) = S(c) + μ₀ < S(c).
Odd: Q^odd = M^odd − 2ss^T ⟹ λ_k ≤ μ_k ≤ λ_{k+1}. The killed Lévy process of §1 therefore has a cluster of decay rates exponentially close
to the arithmetic level S(c) (t5). Upper ends are Galerkin eigenvalues (valid upper bounds); lower ends need window positivity.

## 7. Cheeger and nearest-neighbour (Muckenhoupt) minorants
Cheeger (Lawler–Sokal, with killing): λ₀ ≥ h²/(2d), d = sup_x (killing + total jump rate). κ(z) ~ 1/(2z) has infinite mass, so d = ∞ and the
bound is 0; on N cells d_N = ln N + 0.29 (t6). Dropping all but nearest-neighbour jumps keeps the killing W(x) and a jump part O(h)∫|f′|²,
so the birth–death minorant (where Muckenhoupt/Chen two-sided Hardy criteria live) tends to inf W; W(0) = D + 2K(a) + 8 sinh(a/2) (the
explore32 S2 potential at x = 0 for c < 2) vanishes at c = 1.356980, so this route stops there.

## 8. Reflection identity (parity order)
For g ≥ 0 on [0, a], with Eg, Rg the even/odd extensions: Q(Eg) − Q(Rg) = −4 ∬_{[0,a]²} j(x + y) g(x) g(y). With g the positive odd ground
state (c ≤ 2): λ_ev ≤ λ_odd − 2I(c), I(c) = ∬ J(x + y) g g/‖g‖². I > 0 for c ≤ 1.5 (0.285, 0.210, 0.113 at c = 1.2, ρ, 1.5), < 0 from
c = ρ² (−0.011); there the parity order λ_ev < λ_odd follows instead from Barta(odd) > Rayleigh(even) (t2, t7).
