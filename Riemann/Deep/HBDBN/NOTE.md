# Working notes: audit of arXiv 2610.12803v1 (Hakimoglu-Brown), 2026-10-22

Sources read: the paper's TeX (main.tex, 1789 lines, every line), Polymath arXiv 1904.12438v2 (TeX and PDF numbering), and Polymath's dbn_upper_bound/arb/TloopSinglematv2.c at commit 5fde84e. All downloads were deleted after use.

## Normalisation and notation (as in the paper)
- H_t(z) = ∫₀^∞ e^{tu²}Φ(u)cos(zu)du and H₀(z) = ξ(½+iz/2)/8. A zero z = x+iy corresponds to ρ = (1−y)/2 + ix/2.
- S(x,p,t) = −Im(H'/H)(x+ip) = Σ over conjugate pairs of K_p(d,v), with real zeros at half weight.
- K_h(d,v) = (h−v)/(d²+(h−v)²) + (h+v)/(d²+(h+v)²).
- Zero motion (Polymath Prop 3.1(i)): ż = 2Σ'1/(z−z_k). Hence h' = −1/h − 2Σ'K_h and q' = −2 − 4hΣ'K_h, with q = h².
- Signs checked by hand: ∂_tH = −∂_z²H gives ż = H''/H'. The zero's own conjugate contributes 2/(2ih), i.e. −1/h vertically.

## Chain from "Re ρ ≤ 7/8" to Λ ≤ 0.13963
1. (eq:M0) M(0) ≤ 9/16. This is the only use of the hypothesis.
2. (eq:wall) Polymath §8.4 claim (a) is the enlarged barrier box [X, X+1]×[0.2, 1]×[0, 0.2], with X = X₀ − ½ and X₀ = 6·10¹⁰ + 83952. It is stated for H_t/B_t ≠ 0 and established by Rouché with f_t and an error allowance of 1.25e-3. The wording matches what the paper uses. Polymath used Platt's 3.06·10¹⁰ for its hypothesis (i). The paper uses full RH up to (X+1)/2 = 3.0000042·10¹⁰, which Platt also covers.
3. Lemma head (finite region, |x| ≤ X, t ≤ 1/5): no zero reaches y = 1/5. A zero cannot enter through x = X (wall) or through x = 0 (H_t(iy) > 0). The first contact is at height exactly 1/5, and there h' ≤ −1/h < 0. Multiple contacts are excluded by backward vertical Hermite splitting. Prop 3.1(ii) gives the roots √2(t−t₀)^{1/2}λ_j, which are imaginary for t < t₀. Sound.
4. Comparison (old form, p = 4): K_h ≥ (h/p)K_p when v ≤ h and p² ≥ 5h². The numerator is (B−A)[D² + (A+B−6r)D + (A−r)(B−r)], and A+B−6r ≥ 6(A−r) ≥ 0. Sound.
5. Lemma approx (3 ≤ y ≤ 6, x ≥ 6·10¹⁰): structure F = H/B ≈ P_N(s_B) + γP_N(s_A).
   - Summands: Polymath Prop 6.1, valid for every real σ and T > 10.
   - Remainders: Polymath Prop 6.2. The source (TeX lines 684–720) says "σ real, T' > 0, K positive integer". For σ < 0 the remainder needs K+σ ≥ 2. Prop 6.3 itself is restricted to 0 ≤ σ ≤ 1, which is why the paper rederives the σ ∈ [−3, 4] case. The rederivation follows Polymath's proof line by line.
   - Checked inequalities: the E(U) exponent bound ((U²+6|U|+2/3)/4 ≤ U²+1, because the discriminant is −4); the Gaussian moments E[e^{bU+cU²}] = D^{−1/2}exp(…); and every scalar (evidence/partD.out).
   - Normaliser facts: Re α(w) > L/2, since Re 1/(2w) > 0, Re 1/(w−1) > 0 for y > 1, and |w| > x/2. For |γ|: log|M₀(v)/M₀(w)| = −(y/2)L + O(y³/x²). The heat factor has real part ≈ −tπy/(16T) < 0. So the paper's e^{1e-8·y} is safe; Polymath's own e^{0.02y} is only needed for small x.
6. Prop two-density (p = 61/20, r₀ = 1/20). The disk lies in 3 ≤ y ≤ 3.1. Cauchy bounds: η₁ = 2e-6/r₀ and η₂ = 2·2e-6/r₀² are correct.
   - δ_j: terms n ≤ N₁ use u = σ_p + 5.57t⁻. The terms N₁ < n ≤ N₀ and n > N₀ are disjoint and are covered by n^{−u₈} through an integral tail.
   - ρ_j: |λ_n| ≤ L+2 < 2log N₀ + 3, and monotone decrease in N₀ (exponent ≤ −1.025, log factor < 0.16).
   - Formulas: (F'/F)' = s_B''ϱ + (s_B')²ϱ' + ε''/(1+ε) − (ε'/(1+ε))², so G* is right. −Im(B'/B) = ½Re β(w) ≥ L/4 − 1e-6.
7. Lemma kernel-derivative: proved (algebra re-derived by hand; Φ, K_p and ∂_pK_p re-expanded in sympy inside partA.py).
8. Theorem derivative-step:
   - k_i formula: the own-pair values K_p(0,h) = 2p/(B−q) and ∂_pK_p(0,h) = (2w−4B)/w² with w = B−q.
   - Barrier: the perturbed Q_ε gives a strict first contact. A highest zero is attained for t > 0 because Polymath Thm 1.5(i) (Theorem "Distribution of zeroes", (i): x ≥ e^{C/t} ⇒ y = 0) leaves finitely many non-real zeros.
   - The contact zero has q ≤ q_i by the strip theorem. That is the certificate's r-range, and k_i is evaluated at q_i (the bracket decreases in q).
9. Final: t₂₄ = 23/200 + 463/100000 = 11963/100000, q₂₄ = 1/25, so Λ ≤ t₂₄ + q₂₄/2 = 0.13963.

## Numeric certificates (evidence/)
- partA.py (system python3 + sympy, about 13 s): exact fractions from the printed tables give ε₀, C^cert, G^cert, k_i and the E₂₀ gates. Bernstein: Φ → Ψ by compactification, 8 pieces in U, 25 coefficients each. All 4800 coefficients exceed 13/10000. The minimum equals the printed fraction exactly: row 21, piece 4, (a,b) = (3,4).
- partB.py (python-flint 0.9.0, Arb at 128 bits, about 17 s, 33 MB): recomputes Z_*, ψ_*, ψ_{2,*}, δ_j (sum to 65536 plus the integral tail at u₈) and ρ_j for all 24 rows. Same for the triangle table (U_Z, U_J) and the Euler table (U_Z, U_ψ, U_0, U_1). Every printed entry is a strict upper bound; the entries are clearly upward roundings of the true values.
- partD.py: every scalar inequality in §§9–9.2, the threshold brackets, the Gaussian checks, the tail allowances, and the gates and margins of the older tables.
- sanity.py (mpmath, 30 digits, prime powers up to 10⁶): exact t = 0 values. At x = X it gives S = 5.339 and ∂_pS = 0.147. X sits near the worst case for S because Polymath chose it where the Euler product is large.

## Polymath (91) and the barrier code (outside the paper's own argument)
- (91) as printed divides the time term by 2n. That is wrong, and the paper's footnote says so correctly.
- Code at lines 911–976: num = ⌈D_z⌉ points per edge, spacing 1/(num−1), step Δt = (min|f| − ½)/D_t, abort if min|f| < 1.
- The mesh term D_z/(2(num−1)) is in (½, 0.54] for the mesh sizes reported (11076 down to 56 points). So the forward bound min|f| − D_tΔt − mesh term is about −0.04…0 at the very end of each step.
- The node at the next step covers that end backwards: min|f| ≥ 1, minus a mesh term ≤ 0.54, minus D_t·δ ≤ 0.037. The barrier claim therefore stands, except possibly the tail of the last step near t = 0.2, which the paper does not use (it needs t ≤ 0.1197).
- Not replayed: the code's derivative bounds. For example, generate_ddzbound integrates instead of summing over n plus an n = 1 term. This is a published input, cited as such by the paper.

## Verdict inputs
- No false or unproved step inside the paper.
- Premises: the 7/8 statement, plus Polymath (claim (a), Thm 1.5(i), Props 3.1/6.1/6.2, eqs (37)–(41)) and Platt (2017). All equation and section numbers the paper cites match arXiv v2.
- Our REPORT5 (0.148409) and REPORT6 (0.147374) rest on the extra premises A1–A3 and give a larger bound, so neither is a record.
