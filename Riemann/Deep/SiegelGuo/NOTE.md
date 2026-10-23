# NOTE: audit of Guo, "An explicit elementary proof of no Landau–Siegel zeros" (arXiv 2610.30476v1)
Source: https://arxiv.org/e-print/2610.30476 (Elementary_Frobenius_Siegel_v19.tex, 1387 lines). Line numbers below are TeX lines.
Claim (l.90-98): (1-β) log q > 1/100 for every real zero β of L(s,χ), χ primitive quadratic, conductor q ≥ 3 (appendix: 1/64).

## Chain and step verdicts (all checked line by line; no false or unjustified step found)
Preliminaries (l.252-407)
- (prelim:numerical) l.281-297: atanh-series enclosures with the stated (x,m) reproduce every bound exactly (check_certificates.py).
- Chebyshev ψ ≤ x log 4 (l.329-337): binomial argument valid for ψ (each p-adic digit term of binom(2m,m), binom(2m+1,m) is 0 or 1, and is 1 for every prime power in the window).
- M(x) ≤ log x + log 4 (l.338-342): nM(n) ≤ log n! + ϑ(n). Correct.
- Σ_{odd p} log p/p² < 30607/88200 (l.343-357): f'' = (6 log t − 5)/t⁴ > 0, Hermite–Hadamard on residue classes 11, 13 mod 6. Exact.
- (prelim:von-mangoldt-bounds) l.370-374 and the interval bound (l.387-406): A(x) ≥ log x − 1 uses log(1+1/n) ≤ 1/n; W(x) ≤ V, W(y) ≥ log2/4 + log3/9 (y ≥ 9). Loss = 2 − 7/675 exact.
Arithmetic part (l.410-845)
- Lemma arith:counting (l.431-491): image mod p^ν has ≤ Π p^{ν−h_i} classes (coordinates in the adapted basis mod p^{ν−h_i}); CRT; box of side d0/(n w_i) injects (sup-norm of a difference < d0). Correct.
- ‖a+bt‖ submultiplicative (l.505-517): |ac+dbe| + √|d||ae+bc| ≤ (|a|+√|d||b|)(|c|+√|d||e|). Correct.
- Parameters (l.519-533): m² > p0 (prime not a square), 2(m−1) < p0. Correct.
- Lemma arith:injectivity (l.557-638): f_{u,v} irreducible mod p0 since (d/p0) = −1, v ≢ 0; distinct (v+v' < p0); f_{u,v} | h in Z[ξ] (remainder (au+b)+avt = 0); ξ^{p0} ≡ 2u−ξ mod (p0, f_{u,v}) by Fermat + Euler; d/dξ f(ξ,ξ^{p0}) ≡ f_z mod p0 iterates; k! invertible (k < s0 < p0), gcd(f, f') = 1; degree 2m(m−1)s0 ≥ m² s0 > p0 s0 > deg g. Correct.
- Triangular basis (l.654-667): unitriangular in r, stays in the box (h ≤ ⌊r/p⌋, s+h ≤ s0−1). Correct.
- (arith:frobenius), (arith:divisibility) (l.669-692): ρ(t^p) = d^{(p−1)/2} t exactly; every Leibniz factor ∂^{j}(z^p − w) has remainder ≡ 0 mod p (j = 0 by Frobenius, 1 ≤ j ≤ p by p!/(p−j)!, j > p zero). Correct.
- Coordinate bound (l.695-704): ‖u±vt‖ ≤ (m−1)(1+√|d|) < b0; r!/(r−j)! < p0^j; b0 ≤ p0 from p0 ≥ 4q. Correct.
- (arith:H-upper) (l.705-734): Σ(r+s) = p0 s0(p0+s0−2)/2; λ = (s0−1)/(p0−1) < 1/30; coefficient of λ = (3/2)log p0 − (1/2)log q − log 2 > 0; (29/30)·(1/2) + 1/15 = 11/20. Correct.
- (arith:min-sum), (arith:floor-sum), (arith:floor-square), (arith:H-lower) (l.736-772): identities re-derived ((k−s0)(k−s0+1) ≥ 0; pa(a−1)/2+ab formula); the k² loss uses Σ_odd log p/p² < 7/20. Correct.
- Remainder (l.775-838): every term re-computed exactly; R < 599/100 < 6 (tightest margin of the paper, 0.01, rigorous). Small cases p0 < 10⁴, p0 < 4q, no odd inert prime: correct.
Real-variable part (l.846-1024)
- |A_u(v)| ≤ min(r, q−r) ≤ q/2 (l.862-868), (real:tail) (l.870-877): correct.
- a = 1*χ ≥ 0, multiplicative, a(m²) ≥ 1 (l.879-900): correct.
- Prop real:hyperbola (l.902-1014): hyperbola split correct (d > u ⇒ m < u); first line ≤ (q/2)t^{1/2−β} uses L(β)=0 and t/d ≥ √t ≥ 1; Euler remainder ∈ [−u^{−β}, u^{−β}]; main terms (t/m)^ε m^{−β} = t^ε/m; E0 ≤ (1+½log t)t^{1/2−β}; pairing: E1 = u^ε Σ_{m>u} χ(m)g(m), g = (1/v)∫_u^v w^{−β}, g ≤ ½u^{−β} for β ≥ ½ via (√x−2)² ≥ 0, unimodal ⇒ ∫|g'| ≤ u^{−β}; |E1| ≤ (q/2)t^{1/2−β}. Uniform in q and β; no 1/ε. Correct.
- Lemma real:large-prime (l.1041-1134): (n,p) ↦ np injective, a(np) = a(n)a(p); Abel summation main term = Lκ exactly (sympy); |ε log t − 1| ≤ 1 from ε log(xy) ≤ 1; ∫_y^∞ B/t = 2B(y)+2y^{−1/2} (sympy); L(1)/ε = y^{−ε}(H_β(y)−E(y)) ≤ H_1(y)+B(y) is an UPPER bound only (no lower bound for L(1,χ), no sign assumption, nothing circular); inert-lower uses p > y ≥ q ⇒ χ(p) = ±1. Correct.
Final proof (l.1136-1221) and appendix (l.1232-1278): ε log(xy) = 76δ, β > ½, y ≥ 729; H_1(y) ≥ Σ_{m≤27} m^{−2} > 8/5; f_k' as printed is CORRECT (sympy; our earlier scout report of a wrong derivative at l.1164-1168 was a false positive); κ = 2870(log q)², all certificates exact. Appendix (q ≥ 7, x = q^46, y = q^4) and q = 3,4,5 positivity blocks: correct.

## Brief's watch list
Circular use of β < 1 − c/log q: none (only β ≥ ½, derived from δ ≤ 1/100). Non-uniform errors: none (every bound explicit in q, t; none has 1/ε). Brun–Titchmarsh/sieve: not used. Lower bound for L(1,χ): not used. Hyperbola main-term signs: re-derived, correct.

## Data tests (scripts in this folder; outputs *.out)
- check_certificates.py: 56/56 exact checks (all displayed rationals, log enclosures, f_k', two integrals, Abel main term, kernel g).
- check_prelim.py (x ≤ 10⁷): max ψ(n)/n = 1.039; max M(x) − log x = −0.347; min M(x)−M(y)−log(x/y) = −0.194 (claim > −2); V = 0.7554.
- check_arith.py: 1820 fundamental discriminants |D| ≤ 3000, x ≤ 10⁶: no violation; largest S₋ − (11/20)log x − (29/60)log q = −0.566 (bound 6). Exact H at p0 ≥ max(10⁴, 4q), 8 characters: H ≤ counting bound, displayed lower/upper forms bracket H.
- check_inject.py: injectivity analog (p0 ≤ 97, s0 up to p0−1): 251/251 full rank mod 2³¹−1; divisibility: 774 basis vectors exact, 0 failures; Frobenius 760/760; negative control (split p) fails 12/12, so the test is sensitive.
- check_real.py (9 characters, β ∈ {.5,.6,.75,.9,.99,.999}, t ≤ 3·10⁵; proof run without L(β)=0, extra term L(β)G(t)): max |E|/bound = 0.112; each of the three pieces ≤ 0.411 of its bound; hyperbola identity error ≤ 3·10⁻¹⁴; Lemma (126 (x,y) pairs): |R|/bound ≤ 0.029, J − rhs ≤ −7.2.
- check_model.py (Siegel pattern: all odd primes ≤ p0 inert, q = 3): exceeds the counting bound from p0 ≈ 10⁶ (9.55 > 8.80), true χ₋₃ stays at 4.84: the counting step is what excludes the Siegel pattern, and it holds for real characters only via the Frobenius divisibility.

## Context
Same Frobenius mechanism as OpenAI (1 Oct 2026) and Kerr 2610.21099 (audited correct here, Riemann/Deep/Siegel). The arithmetic part here is self-contained (counting lemma + Stepanov-type injectivity, van Dijk style) and was re-verified independently. If correct, c = 1/100 (1/64 for the appendix) beats the stated all-conductor constants of Banks (0.003) and Kimura (5.8·10⁻⁴). Not our result: nothing for DISCOVERIES.md or REFUTATIONS.md.
