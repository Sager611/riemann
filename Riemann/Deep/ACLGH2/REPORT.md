BARRIER: the onset grows at least like the angle-weighted number of coupled directions, n* ≥ N_S(κ♮) ≥ Σ_{κ_n>κ♮} κ_n² ≈ 2∏_{p∈S} E_θ e^{δ_p(θ)} = exp((½ + o(1)) log² λ²), faster than any power of |S| or λ and slower than exponential. The reason: the angle κ♮ at which the low-frequency weight falls below τ does not depend on S (κ♮ ≈ τ/A_R, A_20 ≈ 8, A_5 ≈ 2, flat for |S| = 6…26), while the pair has a coupled direction at every joint near-resonance τ log p ≡ 0 (mod 2π), with κ² ≈ (area of its coupling pocket)/2π. Theorem condpos then needs n₀ ≥ 2σ²F_S and |M| ≥ 2λ²σ²F_S, so the polynomial sizes of Remark condpos(b) and Prop. reduction (τ = λ^{−A}, dim V' = O(Aσ² log λ)) fail. This rests on a semiclassical count, checked in bands |τ| ≤ 1200 (|S| ≤ 6) and ≤ 589 (|S| ≤ 16) and at single pockets up to τ = 2·10⁴ (|S| = 6). Theorem condpos at fixed λ is untouched.

# ACLGH2: onset of the localisation (H2ii) of Andrews–Connes–La Gatta (arXiv 2611.13129) against the number of primes (aclgh2, 2026-11-11)

Source: TeX `prolate-weil-fine15.tex`, (H2ii) l. 1255–1261, Theorem condpos l. 1318–1330, Remark condpos(b) l. 1335, strategy and Prop.
reduction l. 1364–1380. Models and identities: [DERIVATION.md](DERIVATION.md) (builds on [../ACLGFals/DERIVATION.md](../ACLGFals/DERIVATION.md)).
Notation: S = {∞} ∪ {p ≤ P}, t = 1 unless stated, R = 5, 20, τ = 10⁻², 10⁻³; onset n*(S, R, τ) = min{n₀ : T(n₀, R) ≤ τ}.

## 1. Model
- Cut model (aclgfals's Nyström kernel 2tC_S(t²uw), Laurent series cut at n ≤ nmax): code/semi_pair2.py.
- Band-limited model (new): exact symbol u_∞∏u_p on |τ| ≤ τ_1 = 0.75τ_max, unitary Laurent operator on a grid Δx = π/τ_max, Γ_t =
  Hankel matrix (c_{i+l+1}); an exact pair of projections; no cut of the p-adic Laurent series. code/laurent.py, laurent_count.py.
- Onset identity (DERIVATION §3): T(n₀, R) = W(κ_{n₀}, R), W(κ, R) = Σ_{κ_n<κ} κ_n m_n(R). So n* = N_S(κ♮): the low-frequency weight W
  fixes an angle κ♮(S, R, τ), and the counting function N_S(κ) = #{κ_n > κ} fixes its index; N_S(κ♮) ≥ Σ_{κ_n>κ♮} κ_n² as κ_n ≤ 1.
- Semiclassical count (DERIVATION §4): N_sc = (1/2π)∫(D_t)_+ dτ, D_t = 2 log t − 2ϑ' + Σ_p δ_p(τ log p) (the delay = phase derivative of
  the symbol); by Kronecker–Weyl, N_sc = 2t²F_S, F_S = ∏_p E_θ e^{δ_p(θ)}, log F_S = Σ_{p≤x}(log p)²/p + O(1) = ½ log² x + O(log x).
  Each coupling pocket {D > 0} of area A carries one direction per sign of τ with κ² ≈ A/2π (thin pockets, measured) or ≈ A/2π directions
  with κ ≈ 1 (fat ones): Σκ² over coupled directions ≈ N_sc.

## 2. Validation
| check | result |
|---|---|
| S = {∞}, κ_0…κ_4 = 1 − δ_n(1) = .999971376, .979484735, .524085896, .058976589, .002732329 | cut model: rel. error ≤ 9·10⁻⁸; band-limited: ≤ 3.7·10⁻² (τ_max 100), ≤ 2.3·10⁻³ (400), O(Δx²) |
| S = {∞,2}, nmax 128: semi_pair2 against aclgfals sp_p2_n128 | κ to 4·10⁻⁶; window sums equal to 3 digits; tabulated against direct kernel: κ to 7·10⁻¹¹ |
| S = {∞,2,3}, nmax 256: against aclgfals sp_p23_n256 | window sums equal to 3 digits (R = 20: 1.21, .963, .179 on (10,20], (20,50], (50,100]) |
| S = {∞,2}: band-limited (τ_max 400) against cut (nmax 1024), windows (10,20], (20,50] | R = 20: .0241/.00227 against .0237/.00215; R = 5: .00256/.000612 against .00248/.000578 |
| W(κ, R) and κ♮ against τ_max 100/200/400, depth 20/30/45, closure b₁ .6/.75/.85 | S = {∞,2,3,5,7}, R = 20: W(.01) = .0553/.0554/.0554 (depth), .0554/.0550/.0551 (τ_max); κ♮ changes ≤ 3% |
| N(κ > ½) in the band against N_sc(τ_1) | Table 3.4, within 20% for |S| = 4…16 |

## 3. Results
### 3.1 S = {∞, 2}: a low-frequency direction at κ = 0.023 (corrects aclgfals §2.4)
| model | κ | m(20) | m(5) | index in the model |
|---|---|---|---|---|
| band-limited τ_max 100 / 200 / 400 (depth 25) | .0232 / .0236 / .0232 | .90 / .80 / .90 | .62 / .55 / .62 | 21 / 37 / 75 |
| band-limited τ_max 200: b₁ .6, b₁ .85, depth 40 | .0234, .0233, .0236 | .90, .90, .80 | .62, .62, .55 | 33, 41, 37 |
| Laguerre (aclgfals hk_p2, exact symbol) M 1000 / 1500 | .0210 / .0218 | .87 / .86 | .59 / .59 | 21 / 23 |
| cut nmax 256 / 512 / 1024 | .046 / .037 / .032 (Aitken .022) | .75 / .89 / .81 | .18 / .32 / .43 | 129 / 257 / 261 |
Three models agree. Its true index is N(.0232) = (C/.0232)², with C = lim κ_n√n = .88 ± .05 (band-limited: .74–.98 for n = 8…70; cut
1024: .84–.96 for n = 16…128), i.e. 1.3–1.6·10³. So T(n₀, 20) ≥ .021 and T(n₀, 5) ≥ .014 for n₀ < 1.3·10³. The weight below it is
.0004 (R = 20) and .0003 (R = 5): the next low-frequency direction is κ = 5·10⁻⁴, m(20) = .86, at index ~3·10⁶. Hence
n*({∞,2}, R, τ) ≈ 1.4·10³ for both R and both τ, against ~10–100 read off the cut model. aclgfals's "T ≈ 0.006 R^0.8/n₀ (20 ≤ n₀ ≤ 100)"
is 30× (R = 20) to 60× (R = 5) too low at n₀ = 100: the cut model put this direction near its cut and it was read as an artifact.

### 3.2 The cut model fails for ≥ 3 primes
κ_0 at nmax 128/256/512: S = {∞,2,3,5}: .743/.822/.880; {∞,2,3,5,7}: .644/.703/.729; {∞,…,11}: .631/.687 (nmax 256).
Band-limited: κ_0 = 1.00000 and κ_1 = .99879, .99936, .99960. The top direction sits at mean depth 6.1, 7.3, 8.3 in x = log u, beyond
log nmax. Window sums of the cut model change by factors 2–10 per doubling of nmax (logs/p235*_n*, p2357*, p235711*). The cut model
cannot give error bars there. Every multi-prime number below comes from the band-limited model.

### 3.3 Angle threshold κ♮(S, R, τ) (band-limited, t = 1; τ_max 200 for |S| ≤ 6, else 100 with checks at 200)
| \|S\| (P) | R=5: τ=10⁻² | R=5: 10⁻³ | R=20: 10⁻² | R=20: 10⁻³ | W(.01, 20) | W(.001, 20) | Tr(M_χ N), R = 5 / 20 |
|---|---|---|---|---|---|---|---|
| 1 (∞) | .0589 | .0027 | .0589 | .0027 | .0023 | .0001 | 2.44 / 2.53 |
| 2 (2) | .0232 | .0232 | .0232 | .0232 | .0004 | .0004 | 2.47 / 4.95 |
| 3 (3) | .0141 | .0141 | .0141 | .0141 | .0001 | .0001 | 3.33 / 7.19 |
| 4 (5) | .0211 | .0211 | .0115 | .0101 | .0006 | .0001 | 3.69 / 8.87 |
| 5 (7) | .0189 | .0034 | .0028 | .0015 | .055 | .0008 | 4.29 / 9.84 |
| 6 (11) | .0068 | .0013 | .00125 | .00029 | .069 | .0074 | 4.45 / 10.67 |
| 7 (13) | .0035 | .00046 | .00144 | .00018 | .097 | .0068 | 4.76 / 11.16 |
| 10 (23) | .0061 | .00043 | .00124 | .00009 | .088 | .0083 | 5.87 / 12.16 |
| 16 (47) | .0055 | .00046 | .00123 | .00011 | .094 | .0080 | 7.64 / 14.95 |
| 26 (97) | .0042 | .00045 | .00112 | .00009 | .100 | .0071 | 9.29 / 17.00 |
All 21 sets with |S| = 6…26 lie in the ranges R = 20: κ♮(10⁻²) ∈ [.00093, .00146], κ♮(10⁻³) ∈ [.00009, .00029]; R = 5: [.0029, .0068],
[.00027, .00131]. So W(κ, R) ≈ A_R κ near 0, uniformly in S, with A_20 ≈ 8 (7–11) and A_5 ≈ 2 (1–4). In angle form (H2ii) holds,
polynomially in R and uniformly in S: Tr(M_χR N^{(κ)}) ≈ A_R κ (A_20 from κ♮(10⁻³): 3–11). t = 1.5, R = 20, τ = 10⁻² (10⁻³), |S| = 1…6: .0285 (.0020), .083 (.0052),
.102 (.0037), .0178 (.0100), .0031 (.0020), .0018 (.0004): the same picture.

### 3.4 Counting function: resonance directions
In the band-limited model (S = {∞,2,3,5,7}, τ_max 800) the 21 directions with κ > ½ sit at |τ| ≈ 0 (κ = .99, .98, .63) and, in couples,
at τ ≈ 23, 45, 54, 74, 109, 172, 281, 371, 480, with quartile width ≤ 1 (code/centroid.py). At τ = 280.9: τ log p/2π = 30.99, 49.11, 71.95,
87.00 for p = 2, 3, 5, 7. These are coupled directions at joint near-resonances. Their number grows with the band as the semiclassical
count predicts:
| \|S\| | N(κ > ½) / N_sc(τ_1) at τ_1 = 300 | 600 (589) | 885 | 2F_S = N_sc(∞) | horizon 2πe^{M_S} |
|---|---|---|---|---|---|
| 4 | 13 / 13.1 | 15 / 16.1 | 15 / 16.9 | 25.0 | 4.9·10⁴ |
| 5 | 17 / 17.5 | 21 / 22.1 | 25 / 24.9 | 51.8 | 5.2·10⁵ |
| 6 | 25 / 21.6 | 31 / 28.7 | 31 / 32.5 | 96.9 | 4.1·10⁶ |
| 7 | 20 / 25.1 | 28 / 34.5 | – | 174 | 2.9·10⁷ |
| 10 | 39 / 34.6 | 59 / 49.3 | – | 769 | 5.4·10⁹ |
| 16 | 40 / 50.1 | 66 / 73.9 | – | 6.8·10³ | 2.7·10¹³ |
At τ_1 = 1200 (depth 12): |S| = 4: 15 against 18.1, |S| = 5: 25 against 27.1. At t = 1.5, |S| = 5, τ_1 = 600: 48 against 39.5
(2t²F_S = 117). Complete counts: S = {∞,2}: 4 against 4.9; S = {∞,2,3}: 7 against 10.5 (semiclassics high by 1.5× at |S| ≤ 3).
Beyond |τ| ≈ 2000 every pocket is thin. S = {∞,2,…,11} has 970 pockets on 0 ≤ τ ≤ 5·10⁶ (past the horizon 4.1·10⁶), only one with
A > 2π (at τ ≈ 0); median A = .13, 90% .84; Σ A/2π = 48.7 (× 2 signs = 97.3 = 2F_S). Pockets for |S| = 3, 4, 5: 69, 188, 439 (Σ A/2π
= 5.6, 12.3, 25.7). Shifted-band model (code/window_count.py: modulation by e^{iτ_c x} commutes with P; window τ_c ± 75):
| τ_c (one pocket) | A/2π | largest κ | √(A/2π) | Σκ² in window | N_sc in window |
|---|---|---|---|---|---|
| 2447.6 | .428 | .631 | .654 | .756 | .587 |
| 3099.8 | .444 | .647 | .666 | .764 | .645 |
| 16751.6 | .265 | .480 | .515 | .300 | .265 |
| 19851.5 | .210 | .425 | .458 | .275 | .244 |
| 5·10⁴, 6·10⁵, 1.5·10⁶ (no pocket) | 0 | .18, .053, .024 | 0 | .092, .007, .003 | .010, 0, 0 |
So κ ≈ .95√(A/2π) and Σκ² ≈ N_sc + a ladder background ≤ .09. Centred at τ_c = 0 the window model reproduces the real model's counts.

### 3.5 Growth: polynomial or exponential
Lower bound n* ≥ Σ_{κ_n>κ♮} κ_n² ≈ 2F_S (S = {∞} ∪ {p ≤ x}; code/arith_factor.py, quadrature for every p ≤ 10⁶):
| x | 11 | 13 | 23 | 47 | 97 | 199 | 997 | 9973 | 99991 | 999983 |
|---|---|---|---|---|---|---|---|---|---|---|
| \|S\| | 6 | 7 | 10 | 16 | 26 | 47 | 169 | 1230 | 9593 | 78499 |
| 2F_S | 97 | 174 | 769 | 6.8·10³ | 9.0·10⁴ | 3.2·10⁶ | 3.4·10¹⁰ | 3.1·10¹⁸ | 6.5·10²⁸ | 2.9·10⁴¹ |
| log F_S − ½ log²x | +1.0 | +1.2 | +1.0 | +.7 | +.2 | +.3 | −.3 | −.5 | −.6 | −.7 |
Fits for |S| ≤ 26. Exponential, log F = a + b|S|: b = .40 per prime, rms .64, residuals concave (−1.04 … +.85 … −.81). Power,
log F = a + b log|S|: b = 3.9, rms .45, residuals convex (+.82 … −.43 … +.72). Local rates: d log F/d|S| falls .87 → .26 (|S| = 26) →
.0003 (78499); d log F/d log x rises 2 → 3.6 (x = 97) → 13.5 (x = 10⁶). So the growth is quasi-polynomial, F_S = exp(½ log²x + O(log x)),
i.e. exp(2 log²λ + O(log λ)) for x = λ².
Onset estimates for |S| ≤ 6:
| \|S\| | κ♮: R=5 (10⁻², 10⁻³), R=20 (10⁻², 10⁻³) | in-band index of κ♮ (τ_max 1180) | 2F_S ≈ Σκ² | 2 × pockets | ladder estimate of n* = N_S(κ♮) |
|---|---|---|---|---|---|
| 2 | .0232 (all four) | ≥ 209 | 4.7 | – | 1.4·10³ (1.3–1.6·10³) |
| 3 | .0141 (all four) | ≥ 497 | 11 | 138 | 1.1·10⁴ |
| 4 | .0211, .0211, .0115, .0101 | ≥ 537, 538, 722, 746 | 25 | 376 | 5.8·10³, 5.8·10³, 2.8·10⁴, 3.9·10⁴ |
| 5 | .0189, .0034, .0028, .0015 | ≥ 593, 1112, 1153, 1301 | 52 | 878 | 9.5·10³, 8.6·10⁵, 1.4·10⁶, 6.8·10⁶ |
| 6 | .0068, .0013, .0012, .0003 | ≥ 904, 1443, 1454, 1836 | 97 | 1940 | 1.7·10⁵, 1.4·10⁷, 1.8·10⁷, 6·10⁸ |
"2 × pockets" counts one direction per pocket and sign of τ; all have κ ≈ √(A/2π) > κ♮ except pockets with A/2π < κ♮² ≈ 10⁻⁶.
The ladder estimate is the independent-rung count Σ_{|γ(n)|/√n > κ} 2n (code/counting.py). It matches S = {∞,2} (63, 255 against 63–67,
253 at κ = .1, .05) but undercounts multi-prime S at κ = .1 by 1.4–8× (|S| = 3…6), so treat it as an order of magnitude. Fitted to n*(R = 20, τ = 10⁻²):
exponential, slope 2.4 in log n* per prime (×11, rms .66); power, exponent 8.3 (rms 1.10). That exponential look on |S| ≤ 6 comes from κ♮ falling 20×. Beyond
|S| = 6 κ♮ is flat, so the onset grows like N_S at a fixed angle, which is at least the quasi-polynomial 2F_S.

## 4. What it means for ACLG's route
- (H2ii) at fixed S: consistent. T(n₀, R) = W(κ_{n₀}, R) ≈ A_R κ_{n₀} → 0, polynomial in R; at S = {∞,2} it is a staircase set by
  isolated low-frequency directions.
- Theorem condpos needs ε_σ(n₀) ≤ log σ/(C₀λ), so τ ~ 1/λ, κ♮ ~ 1/(A_R λ), and n₀ ≥ N_{S,σ}(κ♮) ≥ Σ_{κ>κ♮}κ² ≈ 2σ²F_S. Also
  |M| = #{κ_m(λσ) > κ_{n₀}(σ)} ≥ 2λ²σ²F_S, and the codimension 2 + 2(n₀ + 1)|M| ≥ 8λ²σ⁴F_S² = exp((1 + o(1)) log² λ²).
- Remark condpos(b) says the sizes are "both polynomial in λ under a polynomial decay of the angles". The angles do decay polynomially for
  each S (κ_n ~ C_S n^{−1/2}), but C_S² ≥ 2F_S grows quasi-polynomially. Prop. reduction's bookkeeping assumes polynomial dimensions:
  α, β = O(τ^{1/2}λ^{C'}) with τ = λ^{−A}, dim V' = O(Aσ² log λ). With dim V ≥ exp(log² λ²) it needs τ ≤ exp(−c log² λ) (A ~ log λ),
  and the block V' has superpolynomial dimension. Theorem condpos itself, for each λ, is unaffected.
- Mechanism (new): coupled directions of the semilocal pair at joint near-resonances of (log p)_{p∈S}, one per coupling pocket, with
  κ² ≈ A/2π. Their angle-weighted number is the Cramér sum 2t²∏_p E e^{δ_p}, spread up to the horizon 2πt²exp(Σ_p 2 log p/(√p − 1)).

## 5. Novelty
python3 .scratch/lst.py over 2026-10/11 in 13 categories (math.NT, CA, CV, SP, PR, FA, CO, math-ph, OA, NA, DS, MP, GM): no title on
semilocal angles, resonances or onset (hits: 2610.32082 multiband prolate matrices, unrelated). Zenodo "semilocal": 22701239 and 22701808
(Dollard-renormalized semilocal prolate–Jacobi, S = {∞, p} and finite S, constructions), 21906620 (certified positivity c = 14, 17, 19),
23442026 (CCM determinant criteria): none counts angles or treats (H2ii). arXiv search, three queries: only CCsel 2006.13771,
CC 2310.18423 and ACLG 2611.13129; ACLG expects polynomial sizes (l. 1335) and names the angle counting function as open (l. 1260). Our
files: no NOVELTY.md mentions resonances; ExploreHankel (NOGO explore13) is about Nyman–Beurling Hankels. Known ingredients: Landau/Weyl
phase-space counting for Hankel–Toeplitz operators, Kronecker–Weyl equidistribution of (τ log p), Mertens sums. Their combination into
N_sc = 2t²F_S and the barrier for ACLG's sizes were not found. Candidates for the orchestrator: NOGO (barrier, N13); DISCOVERIES after a
reviewer verdict; the correction of aclgfals §2.4 (§3.1).

## 6. Caveats
- Σκ² ≈ 2F_S is semiclassical. It is checked against exact counts in bands |τ| ≤ 1200 (|S| ≤ 6) and ≤ 589 (|S| ≤ 16), within 20%
  (|S| ≥ 4; 1.5× high at |S| = 3), and at four single thin pockets of |S| = 6 up to τ = 2·10⁴ (κ ≈ .95√(A/2π)). The bulk of 2F_S lies
  far beyond for large S (|S| = 16: 74 of 6.8·10³ below |τ| = 589; horizon 2.7·10¹³), where pockets are thinner still; the rule
  κ² ≈ A/2π is not checked there.
- Equidistribution of (τ log p) mod 2π is used up to the horizon; the exact pocket scan for |S| = 6 (Σ A/2π = 97.3 against 2F_S = 96.9)
  supports it there.
- Computed at t = 1 (and t = 1.5); Theorem condpos uses σ > 1 and scale λσ, where semiclassics gives 2t²F_S. W and κ♮ come from directions
  with Mellin mass at |τ| ≤ 40: the band-limited model resolves them, with no dependence on τ_max from 100 to 400.
- The point estimates of n* (10³–10⁹) are heuristic; only the lower bounds (in-band index, Σκ² ≈ 2F_S, 2 × pockets) are claimed.

## 7. Process note
Five early runs broke the 300 MB cap (/usr/bin/time, decimal MB): 302, 451, 496 and 668 MB (p235_n1024, p2357_n1024, p235711_n512,
p235711_n1024, each ≤ 2.6 s before the 1-s poller killed them) and 300.4 MB (p235711_n256, finished). The cause was semi_pair2.py's
kernel table, built in up to 65536 × nterms chunks; it was fixed to ≤ 2M entries per chunk. p2_n1024 (before the fix) peaked at
298.7 MB. All later runs stayed ≤ 245 MB (poller 0.2 s, kill at 270 MB). One
heuristic counting run (N_ovl in §3.5's ladder comparison, 77 s) ran outside the monitor: peak not measured, estimated ≤ 250 MB from its
~10⁶-entry dictionary. One compute process at a time throughout; single-threaded.

RH relevance: N13 (sizes of the conditional route of ACLG: (H2ii), the counting function of the angles, Theorem condpos, Prop. reduction)

## Files
- DERIVATION.md (models, onset identity, semiclassical identity), this REPORT.md.
- code/: semi_pair2.py (cut model; copy of aclgfals's semi_pair.py with tabulated kernel and scale t), hankel.py (copy, helpers),
  laurent.py (band-limited model: κ_n, m_n(R), depth profile), laurent_count.py (all eigenvalues: counting function), centroid.py (Mellin
  location of large angles), window_count.py (shifted band at one pocket), pockets.py (pocket scan), semicl.py (N_sc(T) by integration),
  arith_factor.py (F_S for p ≤ 10⁶), counting.py (ladder estimate),
  windows.py, profile.py, lowfreq.py, natsum.py, ctable.py (tables), run*.sh (monitored runners), batch_*.sh.
- logs/: v_* (validation), p*_n* (cut model), L_* (band-limited; *_s1.5 at t = 1.5; L_P* up to P = 97), C_* (counting), CEN_p2357_t800
  (centroids), W_*, W2_* (single pockets), PK_* (pocket scans), SC* (semiclassics), AF.json (F_S), *.time (wall time, /usr/bin/time
  peak RSS, polled peak), *.killed.* (cap breaches).
