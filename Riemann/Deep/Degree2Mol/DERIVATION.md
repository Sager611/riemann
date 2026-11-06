# DERIVATION — degree2mol

Notation: ν = mollifier length for the degree-2 L-function (ψ has length T^ν), ξ = ν/2 the ζ-normalised length
(2508.11108 l.849; Bernard 1404.6432 Thm th:mm2 uses c(P,Q,2R,ν/2)); θ_KS = Kim–Sarnak exponent 7/64.
Line numbers: AT21 = arXiv 2004.03581 zeros-paper.tex; B = Bernard arXiv 1404.6432 main.tex; KRZ = arXiv 1609.03738
KRZ02_post_referee_new.tex; S = arXiv 2508.11108 1_Latest_9_Aug.tex.

## 1. Inputs (quoted, not re-proved)

1. Mean value (B l.770-783, Thm th:mm2; range l.774): for P ∈ R[x] with P(0) = 0, P(1) = 1 (B l.161) and Q a polynomial with Q(0) = 1,
   V = Q(−(1/(2 ln T)) d/ds) L(f,s), σ₀ = 1/2 − R/ln T:
   (1/T)∫₁ᵀ |Vψ(σ₀+it)|² dt = c(P,Q,2R,ν/2) + o(1) for ν < (1−2θ_KS)/(4+2θ_KS), where
   c(P,Q,r,ξ) = 1 + (1/ξ)∫₀¹∫₀¹ e^{2rs} [d/dx (e^{rξx} Q(s+ξx) P(x+u))|_{x=0}]² du ds.
   The bracket is rξQ(s)P(u) + ξQ'(s)P(u) + Q(s)P'(u) = e^{−rs}(w(s)P'(u) + ξw'(s)P(u)), w(s) = e^{rs}Q(s); this is AT21 l.205-211
   and S eq. (cpqroptim)/(cpqroptimwy) (S l.274-279). Bernard: ν < 5/27 at θ_KS = 7/64, ν < 1/4 at θ_KS = 0 (B l.796). AT21 restates it at l.194-218 (P ∈ R[x], l.194).
2. Range (AT21 Thm thm:nu, eq. (nu-interval) l.229): the same asymptotic holds for 0 < ν < 1/4 − θ_KS/2, for every holomorphic or Hecke–Maass
   newform on Γ₀(N), any N (AT21 Thm, l.138). θ_KS = 7/64 gives ν < 25/128; θ_KS = 0 (generalised Ramanujan for the Maass
   spectrum) gives ν < 1/4.
3. Zero count (B Prop. l.791-795; AT21 eq. (kappa_lower_bound) l.214-218): κ_f ≥ 1 − (1/R) log c(P,Q,R,ν/2) for Q in Bernard's class
   (B eq. (formepolynome), l.786: Q = 1 + Σ i^{n+1}λ_n[(1−2x)^n − 1], λ_n real). For real Q this is exactly
   Q(0) = 1 and Q(y) + Q(1−y) = const (odd polynomial in 1−2y), Conrey's condition, also KRZ l.172.
4. Our Q: Q(y) = 1 + Σ_{n=1}^{30} h_n (T_{2n−1}(1−2y) − 1), h_n rational: T_{2n−1} odd, T_{2n−1}(1) = 1, so Q(0) = 1 and
   Q(y) + Q(1−y) = 2 − 2Σh_n. It is in Bernard's class (finite real combination of (1−2y)^{odd} − 1). cert.py checks Q(0) = 1 and
   Q'(y) ≡ Q'(1−y) as exact polynomial identities. Our P: P = p/p(1), p = Σ_{j odd ≤ 25} λ^j x^j/j!; P ∈ Q[x], P(0) = 0, P(1) = 1 (checked).
   ν: we use ξ = 25/256 − 10⁻⁷, 1/8 − 10⁻⁷, 5/54 − 10⁻⁷, strictly inside the open ranges.

## 2. Expansion (all P, Q real)

c = 1 + (1/ξ)∫∫(wP' + ξw'P)² = 1 + (C/ξ)W₀ + ξB W₁ + 2(∫ww')(∫PP'),  B = ∫₀¹P², C = ∫₀¹P'², W₀ = ∫₀¹w², W₁ = ∫₀¹w'².
∫PP' = (P(1)² − P(0)²)/2 = 1/2 and ∫ww' = (w(1)² − w(0)²)/2 = (e^{2R}Q(1)² − 1)/2, so
  (E)  c = (1 + e^{2R}Q(1)²)/2 + (C/ξ)W₀ + ξB W₁.
For P = x (B = 1/3, C = 1) and Q(1) = 0 this is S eq. (cpr)-(J) (l.322-327): c = 1/2 + J(Q). cert.py does NOT use (E): it integrates
the original integrand (x first, exactly; then y, exactly), so (E) is only used by the untrusted optimiser.

## 3. Optimal P for fixed Q (Conrey's Section 4; AT21 l.996-1005, eq. (kappa-ABC) l.1003, for real Q)

Minimise (C/ξ)W₀ + ξBW₁ over P(0)=0, P(1)=1: Euler–Lagrange P'' = λ²P, λ = ξ√(W₁/W₀), P = sinh(λx)/sinh λ. Then
∫P'² = λ²(1/2 + sinh 2λ/(4λ))/sinh²λ, ∫P² = (sinh 2λ/(4λ) − 1/2)/sinh²λ, and α∫P'² + αλ²∫P² = αλ coth λ (α = W₀/ξ), so
  inf_P c = (1 + e^{2R}Q(1)²)/2 + √(W₀W₁) coth(ξ√(W₁/W₀))     (= AT21 (kappa-ABC) with A = W₀, C = W₁, ν/(2A)·√(AC) = ξ√(W₁/W₀)).
Joint optimum: alternate Q-step (quadratic in h for fixed λ: linear solve) and P-step (λ = ξ√(W₁/W₀)); each step lowers c.
The fixed point solves both Euler–Lagrange equations (Question 2 of S, l.954): P = sinh(λx)/sinh λ and S solving the ODE of §4
with c = −c₀/c₁ for B, C of that P. Complex Q (Bernard's class) with real P: c(Q_r + iQ_i) = c(Q_r) + (C/ξ)∫e^{2Ry}Q_i² +
ξB∫e^{2Ry}(RQ_i + Q_i')² + e^{2R}Q_i(1)²/2 ≥ c(Q_r) (cross terms vanish since P real), so complex Q never helps; with complex P as
well (not covered by B l.161), alternating minimisation from random complex starts returns to the real optimum (logs/cplx.log).

## 4. The S-functional of S §sect:gene (l.858-934), re-derived

y = 1/2 + t/(2R), Q(y) = S(t), Q'(y) = 2R S'(t), e^{2Ry} = e^{R+t}, dy = dt/(2R), β := Q(0) + Q(1) = 1 + Q(1), S(−t) = β − S(t).
W₀ = (e^R/(2R))∫_{−R}^{R} e^t S², W₁ = (e^R R/2)∫_{−R}^{R} e^t (S + 2S')².
∫e^t(S+2S')² = −∫e^tS² + 4∫e^tS'² + 2(e^R S(R)² − e^{−R}S(−R)²), S(R) = β−1, S(−R) = 1 (integration by parts of 4∫e^tSS').
So (C/ξ)W₀ + ξBW₁ = (e^R/(2R))[c₀∫e^tS² + c₁∫e^tS'²] + c₁(e^R/(2R))(e^R(1−β)² − e^{−R})/2,
c₀ = C/ξ − ξBR², c₁ = 4ξBR² (S l.873 ✓). Folding [−R,0] onto [0,R] gives K(S) of S (Kfuncgen) l.864 ✓, hence
(kappa: general) l.876-878 ✓ (checked numerically in cert.py to relative 10⁻⁶⁰ for all three certified (P,Q)).
Euler–Lagrange: c₁(2cosh t S')' = c₀(2cosh t S − βe^{−t}) ⇔ S'' + tanh t S' + cS = cβ/(1+e^{2t}), c = −c₀/c₁ = 1/4 − C/(4Bξ²R²)
(S l.894-899 ✓; e^{−t}/(2cosh t) = 1/(1+e^{2t})). At the optimum, integrating c₁∫2cosh t S'² by parts gives (generfinalK) l.880 ✓.
Boundary data S(0) = β/2, S(R) = β−1 ✓ (geninit l.868). With β free the optimum adds a natural boundary condition; our Ritz basis leaves β free.

c₀ = c₁ case (S eq. (simplipara) l.372: ξ = √(3/5)/R, P = x, β = 1): (e^R/(2R))·4/(5ξ) = 2e^R/(5√(3/5)) = (2/√15)e^R and ξR/3 = √(3/5)/3, so
  (Q5)  c(S) = 1/2 − √(3/5)/3 + (2/√15) e^R I_R(S),  I_R(S) = ∫₀ᴿ e^t(S² + S'²) + e^{−t}((1−S)² + S'²) dt,
for every admissible S (S(0) = 1/2, S(R) = 0, S ∈ H¹). κ > 2ξ/3 ⇔ log c < R − (2/3)√(3/5) ⇔ c < e^{−(2/3)√(3/5)} e^R.
Zero extension: S_a on [0,R_a], S_a(R_a) = 0, extended by 0 to [0,R] (R ≥ R_a) is in H¹ and I_R = I_{R_a}(S_a) + e^{−R_a} − e^{−R}. Then
c < e^{−(2/3)√(3/5)}e^R ⇐ (2/√15)(I_a + e^{−R_a}) + (1/2 − √(3/5)/3 − 2/√15)e^{−R} < e^{−(2/3)√(3/5)}; the e^{−R} coefficient is
−0.2746 < 0, so R = R_b is the worst point of [R_a, R_b] and the term may be dropped when R_b = ∞.
The paper's κ(ξ) uses the minimiser S_R of the strictly convex K*_R (unique solution of the BVP (mainODEbdry) l.433), so
κ(ξ) ≥ κ(S) for every admissible S. Admissible includes H¹ functions such as the zero extension: K*_R is continuous on the affine
H¹ space {S(0) = 1/2, S(R) = 0}, C¹ (indeed polynomial) functions are dense in it, and the infimum is attained at the smooth
BVP solution, so inf over C¹ = inf over H¹ = K*_R(S_R). Independent check (odecheck.py, mpmath ODE shooting of (mainODEbdry) and S eq. (simplcPRQ)):
κ(1/8) = 0.0854501, κ(1/4) = 0.1763812 (paper 0.0854…, 0.176… ✓) but κ(1/2) = 0.3584369 and κ(2/3) = 0.4463668, while S's table
(l.807-840; rows l.828, l.830) prints 0.364… (2/3) and 0.334… (1/2): those two entries are wrong (the direct functional I and (simplcPRQ) agree to 10 digits).
The closed form (eqn: SRT) with C₁(R) (const: C1R), evaluated with mpmath hyp2f1 (closedform.py), equals the shooting solution
(S'(0) = −0.4415648493, S(R/2) = 0.191015319 at θ = 1/2; Wronskian formula (spec: Wronsk) checked) and gives κ(1/2) = 0.3584369085.

## 4b. Families of arXiv 2607.00282 (main.tex l.545-625)

l.545-547: P, Q real polynomials, P(0) = 0, P(1) = 1, Q(0) = 1, Q'(y) = Q'(1−y) (same class as ours); l.558 (kappaprop) and
l.562 (Conreyformconst) are Conrey's c, identical to §1; l.565: κ(θ) := sup over (P, Q, R). Cor. fullresultpgl3per (l.587-603):
on average over the twist family, N₀/N > ϰ_r(ϖ) with ϰ₂(ϖ) = κ((1−ϖ)/(2(1+ϖ))), ϰ₃(ϖ) = κ((1−3ϖ)/(6(1−ϖ))). (ourproportionlev)
l.605-611 then asserts ϰ_r(ϖ) > 2θ/3 at those θ, i.e. Question 5 of 2508.11108 (our q5cert.py). Simple zeros (Cor. simplpergl2 l.619-623;
RatioGL3.tex l.185): Q linear. Monotonicity: for fixed Q, R, inf_P c = (1 + e^{2R}Q(1)²)/2 + √(W₀W₁)coth(ξ√(W₁/W₀)) decreases in ξ,
so κ(θ) is non-decreasing and κ(θ(ϖ)) ≥ κ(ξ_cert) once θ(ϖ) ≥ ξ_cert = 1/6 − 10⁻⁷ or 1/2 − 10⁻⁷ (with δ = 10⁻⁷: θ₃(ϖ) = (1−3ϖ)/(6(1−ϖ)) ≥ 1/6 − δ ⟺ ϖ ≤ 3δ/(1+3δ) ≈ 2.9999991·10⁻⁷, and θ₂(ϖ) = (1−ϖ)/(2(1+ϖ)) ≥ 1/2 − δ ⟺ ϖ ≤ δ/(1−δ) ≈ 1.0000001·10⁻⁷; corrected after review, the earlier text said 3·10⁻⁷ resp. 4·10⁻⁷; the claims are limits ϖ → 0, so only the cutoff wording changes).

## 5. Exact integration (cert.py, q5cert.py)

For a polynomial G and a ≠ 0: ∫₀¹ e^{ay}G(y)dy = Σ_{k≥0} (−1)^k (e^a G^{(k)}(1) − G^{(k)}(0))/a^{k+1} (finite; repeated integration by
parts), and the same with lower limit 1/2 (factor e^{a/2}). In cert.py G = a₁Q² + 2ξa₂QU + ξ²a₃U², U = RQ + Q', a₁ = ∫P'², a₂ = ∫PP',
a₃ = ∫P² (exact rationals), a = 2R. The two sums are exact rationals; the only rounding is exp(a), exp(a/2) and log in Arb
(4000 bits in both scripts: at small R the exact sums reach ~10³⁰⁰ before cancelling to O(1)), outward-rounded.

## 6. Constants and where they live

| constant | value | source | script |
|---|---|---|---|
| ν range Bernard | (1−2θ)/(4+2θ) = 5/27 at θ = 7/64 | B l.774, l.796 | cert.py 'Bernard_uncond' ξ = 5/54 − 10⁻⁷ (params.json) |
| ν range AT21 | 1/4 − θ/2 = 25/128 | AT21 l.229 | 'AT21_uncond' ξ = 25/256 − 10⁻⁷ |
| ν under Ramanujan | 1/4 | AT21 l.229 (θ=0), B l.796 | 'RSC' ξ = 1/8 − 10⁻⁷ |
| AT21 R, h_n | table l.1014-1063 | AT21 | at21_tables.json |
| family lengths | ξ → 1/6 (r = 3), 1/2 (r = 2) | 2607.00282 l.590-593 | cert.py 'family_PGL3', 'family_PGL2', 'family_PGL2_simple' (ξ = 1/6 − 10⁻⁷, 1/2 − 10⁻⁷) |
| κ formula | 1 − log c/R | AT21 l.214-218 | cert.py c_kappa |
| c(P,Q,r,ξ) | Conrey integral | B l.776-781, AT21 l.205-211 | cert.py c_kappa |
| c₀, c₁, c | C/ξ − ξBR², 4ξBR², −c₀/c₁ | S l.873, 894 | cert.py cross_checks |
| √(3/5)/3, 2/√15, e^{−(2/3)√(3/5)} | 0.258199, 0.516398, 0.596666 | §4 above | q5cert.py s35, s15, a0, b0, rhs |
