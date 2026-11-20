BARRIER: value distribution theory makes Levinson's count an exact a-point count of g = (ζ′/ζ)/(χ′/χ) (every value with Re a > ½ is taken N⁺ + O(1) times right of the line, every value with Re a < ½ (|a| not tiny) N⁺ + N̄₀ + O(1) times), but in this hyperbolic half-strip its second main theorem is the first main theorem (Littlewood, Levinson's own step) plus the log-derivative lemma, which bounds zero counts only from below and, at T = 10⁴, with errors larger than the counts (best SMT lower bound −112 against true values −9 to 770, ζ and Davenport–Heilbronn alike) (node N8).

# explore54 — value distribution theory (Nevanlinna FMT/SMT, log-derivative lemma, Ahlfors) → N8, with N13, N9 side checks

## Transfer
- Tool: Nevanlinna's first main theorem (FMT), second main theorem (SMT) and lemma on the logarithmic derivative (LDL), in Littlewood's rectangle form; Ahlfors' length–area ratio L/S; Milloux/Hayman/Clunie are SMT/LDL corollaries (same structure).
- Node N8, Levinson–Conrey branch: κ ≥ 1 − 2N_V⁺/N, N_V⁺ = #{zeros of V_b = F + (b/L)F′, σ > ½}, bounded by Littlewood + Jensen (explore45: the concavity step is 58% of the loss at T = 10⁴).
- Bridge: the Levinson family is the a-point family of ONE meromorphic function. For F(s) = χ(s)F(1−s), F(s̄) = F(s)‾ (F = ζ or DH):
  g := (F′/F)/(χ′/χ); on σ = ½, 2 Re F′/F = χ′/χ, so g = ½ − iη exactly, η := Im(F′/F)/|χ′/χ|;
  V_c := F′ − (½ + c)(χ′/χ)F = (χ′/χ)F·(g − ½ − c); Levinson's V_b ≈ (b/L)V_c with c = 1/b − ½ (χ′/χ ≈ −L; b < 2 ⟺ c > 0).

## S1 (a-point identity; argument principle = Nevanlinna's 'zeros of f′ against zeros of f' bookkeeping)
For a window [T₁, T₂], Re a ≠ ½ and |a| > max |g| on the far side σ = 3 (≤ |ζ′/ζ(3)|/|χ′/χ| ≈ 0.03 for ζ at T = 10⁴) (counts in σ > ½, with multiplicity; N̄₀ = distinct zeros on the line, N⁺ = zeros of F with σ > ½):
  #{g = a} = N⁺ + E (Re a > ½),   #{g = a} = N⁺ + N̄₀ + E (Re a < ½),   E = O(1) edge terms (O(log T) in general).
Proof: on the line Re(g − a) = ½ − Re a has fixed sign while −η runs once over ℝ between consecutive line zeros (g has a simple pole with negative residue 1/(χ′/χ) at each); a gap and the indentation round the next zero contribute 0 or 2π to the winding; the far side σ = 3 has g − a ≈ −a. a = 0 (zeros of F′; Speiser/Levinson–Montgomery): subtract the far-side winding of F′/F (ζ: (log 2/2π)H).
Corollaries: (i) Levinson's first step is lossless: N − N₀^mult = 2N_{V_b}⁺ + O(log T) for 0 < b < 2; (ii) line zeros are a difference of two a-point counts, N̄₀ = #{g = a} − #{g = a′} (Re a < ½ < Re a′); (iii) reflection V_c(1−s) = −V_{−c}(s)/χ(s) (checked to 10⁻³¹, refl.py), and V_c has N zeros in all (c ≠ ±½), so the 'good' count #{V_{−c} right} = N − #{V_c right} (simple zeros): a lower bound for good a-points IS Levinson's upper bound for bad ones (no new freedom); (iv) at c = ½ the reflection reads V_½(s) = F′ − (χ′/χ)F = −χ(s)F′(1−s): Levinson's function at b ≈ 1 (optimal linear Q: b = 1.032, κ = 0.3474) is the reflected derivative, so his bad count is #{ζ′ zeros left of ½} (Speiser) and his Littlewood-weight loss (explore45: 0.106 of 0.249 at T = 10⁴) is the weighted count of ζ′ zeros with ½ < β′ < ½ + R/L, reflected: N8's Levinson branch and N13's Speiser node are one count.

## Tests of S1 (fz.py: Dirichlet series + Euler–Maclaurin, rel. error vs mpmath ≤ 5.5·10⁻¹⁰ (ζ), 2·10⁻¹¹ (DH); integer windings; c1.py, c2.py, c4.py)
| window | F | N̄₀ (sign changes) | N_tot | N⁺ | #{g = a} right, Re a > ½ | #{g = a} right, Re a < ½ |
|---|---|---|---|---|---|---|
| [10, 240.3] | ζ | 102 | 102 | 0 | 0 (6 real a) | 102 (c = −0.49: 76, edge) |
| [10, 240.3] | DH | 152 | 160 | 4 | 4 | 157 (= 4 + 152 + 1 edge) |
| [10⁴, 1.4·10⁴] | ζ | 4808 | 4808 | 0 | 0 (c = 2, 0.5, 0.2; b = 1.03, 1.5) | 4807 (c = −0.2, −0.4; b = 3, 6) |
| [10⁴, 1.1·10⁴] | DH | 1301 | 1437 | 68 | 68 (same 5 values) | 1369 = 68 + 1301 |
| [10⁴, 10500] | ζ | 590 | 590 | 0 | 0 (a = 0.8+0.1i, 2+2i) | 589–590 (a = 0.45, 0.2+0.5i, −1−i, 0.1−0.3i); a = 0: 534 (pred. 534.84) |
| [10⁴, 10500] | DH | 650 | 716 | 33 | 33 | 683–684 = 33 + 650 |
DH attains Levinson's inequality: 1301 = 1437 − 2·68. By S1, ExploreSymm's 275 DH V-zeros right of ½ on [10⁴, 1.4·10⁴] are DH's off-line zeros right of ½ (± O(1)). The c = 0.5 rows are V_½ = −χ(s)F′(1−s): DH′ has 68 zeros left of ½ in [10⁴, 1.1·10⁴], its 68 left off-line zeros (Levinson–Montgomery), ζ′ none (Speiser). Values on Re a = ½ are taken on the line (once per gap) and are ill-conditioned for winding (excluded). DH, a = 0: 616, i.e. the far-side winding of DH′/DH at σ = 3 is 67, not (log 2/2π)H = 55 (its 4^{−s} term competes). Peaks ≤ 196 MB.

## S2 (Cauchy law of η on the line; KNOWN: arXiv 2609.15862)
Littlewood for V_c/F on [½, 3] × [T₁, T₂] gives, with H = T₂ − T₁, W(h) := Σ over zeros of h with σ > ½ of (β − ½):
  Dlog(c) := (1/H)∫ log|c + iη| dt − log(c + ½) = (2π/H)[W(V_c) − W(F)] + O(1/L²) + edges/H   (c > 0).
Under RH both W vanish, so η is Cauchy(0, ½) in log-potential for every zero statistics (Boole's identity: deterministic, no √N fluctuation); an isolated off-line pair at depth δ contributes −(δ − u₊)·2π/H, u₊ = (√(1 + A²δ²) − 1)/A, A = c|χ′/χ| (local model 1/(u−δ) + 1/(u+δ) + A = 0), so off-line zeros concentrate η near 0.
| quantity (c2.py) | ζ, [10⁴, 1.4·10⁴] | DH, [10⁴, 1.1·10⁴] |
|---|---|---|
| KS distance to Cauchy(0, ½) | 1.6·10⁻⁴ | 1.24·10⁻² |
| P(\|η\| < ½) (Cauchy: ½) | 0.4999 | 0.5231 |
| Dlog(0.05), Dlog(0.5), Dlog(2), Dlog(4) | +2.4, +1.4, +0.7, +0.4 (·10⁻⁴) | −0.0664, −0.0418, −0.0182, −0.0103 |
| (2π/H)[W(V_c) − W(F)] at c = 0.5, 2 (W from Littlewood integrals) | 0, 0 | −0.04176, −0.01814 |
DH: W(F) = 11.600 (68 off-line zeros, mean β − ½ = 0.171), W(V_2) = 8.712, W(V_0.5) = 4.954, W(V_0.2) = 2.721. The identity closes to 10⁻⁴. This is 2609.15862's mechanism (Cauchy limit iff the right-of-line distances sum to o(T)); a re-encoding of Littlewood's lemma, not new.

## S3 (the transfer proper: second main theorem in the critical-scale half-strip)
Littlewood form on R_X = [½, ½ + X/L] × [T₁, T₂], values a₁..a_q (δ = min |a_i − a_j|), scale λ, I = ∫ over σ = ½, J = ∫ over σ = ½ + X/L, I⁺ = ∫ log⁺:
  2π[Σ_j W(g = a_j) − W(g′ = 0) − (q − 2)W(g = ∞)] ≥ RHS := Σ_j I⁺(g − a_j) − I⁺(λg′) − I(log⁺Σ_j|λg′/(g − a_j)|) − (q−1)H log⁺(2/δ) − Σ_j J(log|g − a_j|) + J(log|λg′|),
W now weighted by min(β − ½, X/L); poles of g = zeros of F (simple), bad values Re a_j > ½ count only off-line zeros (S1).
Structure lemma (exact algebra, every q, a_j, λ, X and every meromorphic g, so also V_b, VΨ, ζ itself): inserting Littlewood's identities 2πW(g − a_j) = I(log|g − a_j|) − J(log|g − a_j|) + edges and 2πW(g′) = I(log|λg′|) − J(log|λg′|) + edges, all zero counts cancel (Σ_j W(a_j) − W(g′=0) − (q−2)W(∞) = Σ_j W(g − a_j) − W(g′)), and LHS − RHS = I[log⁺(1/|λg′|) + log⁺Σ_j|λg′/(g − a_j)| + (q−1)log⁺(2/δ) − Σ_j log⁺(1/|g − a_j|)] ≥ 0 holds pointwise for arbitrary numbers (w, w′) = (g, λg′): no analyticity. So the SMT is FMT (Littlewood) + arithmetic; its only analytic content in Nevanlinna theory is the LDL, an upper bound for I(log⁺|g′/(g − a)|), which by Littlewood for g′/(g − a) is an upper bound for W(g′ = 0) − W(g = a) (+ pole terms): a LOWER bound for a-point counts. Levinson needs an UPPER bound for N⁺; in this theory only Jensen on the FMT supplies one (Levinson's own step). Bad values enter the SMT only on the large side: it is an RH-disproof inequality, never a proof one.
Numbers (c3.py, window [10⁴, 1.1·10⁴], H = 1000, N ≈ 1181 (ζ), 1437 (DH); grid: 21 values (9 bad, 7 good, mixed), q ∈ {2, 3, 4, 6}, ≤ 60 subsets each, λL ∈ {¼, 1, 4}, X ∈ {½, 1, 2, 4}; best RHS/2π):
| values | ζ: best lower bound | ζ: true Σ_j [W(a_j) − W(∞)] | DH: best lower bound | DH: true |
|---|---|---|---|---|
| bad, q = 2, 3, 4, 6 | −196, −482, −822, −1393 | −0.01 (all) | −205, −493, −822, −1394 | −3.4, −4.2, −7.0, −9.1 |
| good, q = 2, 3, 4, 6 | −116, −410, −567, −1500 | 134, 309, 419, 769 | −112, −400, −563, −1493 | 144, 221, 420, 607 |
| mixed, q = 2, 3, 4, 6 | −250, −586, −948, −1715 | 76, 76, 131, 330 | −245, −574, −929, −1698 | 83, 139, 138, 338 |
Smallest case [14, 240] (102 ζ zeros, 160 DH zeros; logs c3_small_*): best bounds −26 to −388 for both, true values −2.2 to 172. Every bound is negative by 112–1715 at T = 10⁴ (slack ≍ (q − 1)H/2π·O(1), the size of the counts): the SMT cannot detect DH's 68 off-line zeros (W(F) = 11.6) and gives nothing for ζ. RHS/2π is the usable bound: the dropped nonnegative terms (q − 2)W(F) (the unknown itself) and W(g′ = 0) = W(g′) + 2W(F) (critical points of g; W(g′) = 79.7, 158.5, 306.9, 545.5 for ζ at X = ½, 1, 2, 4, DH 83.4, 167.4, 332.6, 598.1, c3w.py: about one critical point per zero right of the line even under RH) are unknown to any unconditional argument, and adding them back the selected configurations stay negative (e.g. ζ good q = 2: −116 + 80 = −37). Sampled grid, not a certified optimum; the structure lemma above is the certified part. Peaks 91 MB (ζ), 123 MB (DH).
Ahlfors (no exhaustion): in the RH model g₀ = ½(1 − coth(z/2)), z = L(s − ½), each period 2π of the half-strip maps once onto {Re w < ½}: covering number S = (1 + 1/√5)/2 = 0.7236, boundary length L = π·2/√5 = 2.810 per period (sphere of area π), L/S = 3.88 for every height, so Ahlfors' (q − 2)S ≤ Σ n̄ + hL is void and the omission of the open half-plane {Re w > ½} (= RH) is permitted.

## Global regimes (why no other domain helps; b12.py, dhz.py)
- Disc (classical theory of ξ): moving every zero horizontally anywhere in 0 ≤ β ≤ 1 changes N(r, 0, ξ) = Σ log⁺(r/|ρ|) by at most Σ_{γ>0} γ⁻² = 0.023105 for every r (each |log|ρ̃| − log|ρ|| ≤ 1/(2γ²); Σ_{γ>0} 1/(¼ + γ²) = ½(2 + γ_E − log 4π) = 0.023095709 with RH to 3·10¹², plus ¼Σγ⁻⁴ < 10⁻⁵), and by ≤ 2.96·10⁻¹² if only zeros above 3·10¹² move (Σ_{γ>T} γ⁻² ≤ (log(T/2π) + 1)/(πT)); the SMT for ξ with values 0, a, ∞ is slack by ≍ N(r, a) ≍ T(r, ξ) ~ (r/π) log r. RH-blind.
- Half-plane σ > ½ (Nevanlinna's half-plane characteristic, Goldberg–Ostrovskii Ch. I): ζ has BOUNDED characteristic, S(∞) = A(∞) = (2/π)∫₁^∞ log⁺|ζ(½+it)| t⁻² dt = 0.0296 (t ≤ 3000) + ≤ 0.0008 (tail, Hiary's |ζ(½+it)| ≤ 0.618 t^{1/6} log t); semicircle term → 0; DH: 0.127. For bounded characteristic the SMT and defect relations are void (error O(log r) exceeds the characteristic), and the FMT for the value 0 is Carleman's formula, i.e. Balazard–Saias–Yor (RH ⟺ ∫ log|ζ(½+it)|/(¼ + t²) dt = 0). Off-line Carleman count C(∞, 0) = 2Σ_{β>½}(β − ½)/γ²: ≤ 5.9·10⁻¹² for ζ; ≥ 1.269·10⁻⁴ for DH (first four off-line zeros 0.80852 + 85.69935i, 0.65083 + 114.16334i, 0.57436 + 166.47931i, 0.72426 + 176.70246i, refound by findroot). Also blind to Lindelöf: A(∞) < ∞ for every μ(½) < ∞.
- Critical-scale strip: the only regime with counts of order N; there S3's structure lemma applies.

## Side nodes
- N13 (Speiser): S1 at a = 0 is Levinson–Montgomery (#ζ′ zeros left of ½ = #ζ zeros left of ½ + O(log T)); it holds verbatim for DH (c4: DH a = 0 count 616 = 33 + 650 − 67, 67 = far-side winding of DH′/DH at σ = 3), so a quantitative Speiser separates ζ from DH only through N⁺ itself.
- N9: Gundersen/Hayman bounds |ζ^{(k)}/ζ(s)| ≤ |s|^{kε} off an exceptional set are weaker than Landau's ζ′/ζ ≪ log² t and blind to growth exponents; Clunie's lemma needs an algebraic differential equation (ζ satisfies none: Hilbert, Ostrowski, Voronin); the half-plane characteristic is finite for every μ(½). No transfer.

## Barrier (unconditional; what is certified)
For N8 (Levinson branch) value distribution theory supplies (a) identities — S1, the reflection V_c(1−s) = −V_{−c}(s)/χ(s), the Cauchy identity S2, Carleman/BSY — all re-encodings of Littlewood's lemma; (b) inequalities — the SMT, which by the structure lemma is Littlewood + a pointwise arithmetic inequality valid for every pair (w, w′), so its zero content beyond Littlewood is the LDL, a LOWER bound for a-point counts (Milloux, Hayman's alternative and Clunie's lemma are SMT/LDL corollaries of the same shape); and Ahlfors' covering inequality, void since L/S = 3.88 per period. Hence every upper bound for N⁺ obtainable here is an FMT bound (Littlewood + Jensen on a base line, i.e. Levinson's method for some F); moreover the SMT terms that must be bounded below (Σ_j I⁺(g − a_j), J(log|λg′|)) get only upper bounds from Levinson's mean-square inputs. Certified: S1 (proof + exact integer windings, ζ and DH), reflection (10⁻³¹), structure lemma (algebra, all q, a_j, λ, X, any meromorphic g), B1/B2 constants. Sampled: S3's slack table, S2's statistics.
Input-class ceiling (explore45, certified there): no inequality valid for all analytic functions plus Levinson's 42 line second moments (leading coefficient 1, bandwidth 1.61) proves ν < 0.200 at T = 10⁴, θ = ½. The branch itself is far below the record: best Levinson–Conrey constant 5/12 against the TEN record 0.6736.

## Adversarial (Davenport–Heilbronn)
Every statement holds for DH with its own N⁺ (68 in [10⁴, 1.1·10⁴]) and for any F with F(s) = χ(s)F(1−s), real coefficients: S1 and S2 tell ζ from DH only through N⁺ itself (0 against 68; KS 1.6·10⁻⁴ against 1.2·10⁻²), and the SMT is vacuous for both. No input of the theory sees the Euler product: no RH information.

## Novelty
- S2 known: arXiv 2609.15862 (Sept 2026; also cited in NOGO R386).
- S1: a = 0 is Levinson–Montgomery 1974 (Speiser 1934); the all-a version and the equality form of Levinson's first step are a routine argument-principle extension, not found as stated: axq 'Levinson-Montgomery' (4 hits), title 'logarithmic derivative' ∧ zeta (14), 'roots of the equation' ∧ zeta, au:Garunkstis (13), Speiser ∧ zeta, Nevanlinna ∧ 'Riemann zeta', 'second main theorem' ∧ zeta (2401.14481 Borel lemma; 1407.4545), Boole ∧ zeta; monthly listings 2026-10/11 and Zenodo (2 phrase queries): nothing. Routine; no DISCOVERIES entry.
- Structure lemma: Nevanlinna's own proof read with Littlewood's identities; folklore-level. Corollary (iv), V_½ = −χ(s)ζ′(1−s), is a two-line consequence of the functional equation; folklore-level.

## Where it breaks / what remains
- Breaks: the half-strip is hyperbolic (no defect relation; errors O(H) like the counts) and the analytic input (LDL) points the wrong way. RH is the statement that g omits {Re w > ½} in the right half-strip, which no value-distribution law forbids or forces.
- Left for N8: the TEN/pair-correlation branch (positivity, outside this field); inside Levinson only explore45's zero-atomicity gap (≤ 17% of the bound at T = 10⁴).
- The SMT is an RH-disproof direction (lower bounds for bad a-points); at T = 10⁴ it does not even detect DH's off-line zeros.
- Pointer (from (iv), not from this field): at b = 1 Levinson's Littlewood-weight loss is the near-line ζ′ mass G = Σ_{½<β′<½+R/L}(R/L − (β′ − ½)); an unconditional lower bound for G would cut that loss directly. The Levinson–Montgomery total Σ(β′ − ½) = (T/2π) log log T + O(T) is too coarse (it gives G ≥ (N/L)(R − log log T + O(1)) < 0); a usable bound needs near-line ζ′ zeros from small gaps (pair correlation), i.e. the TEN input class.

## Suggested NOGO line
- explore54 (Orca, exploration: value distribution theory — Nevanlinna FMT/SMT, log-derivative lemma, Ahlfors covering surfaces, Milloux/Hayman/Clunie → N8, N13, N9) BARRIER: Levinson's family is the a-point family of g = (ζ′/ζ)/(χ′/χ) (values with Re a > ½ taken N⁺ + O(1) times right of the line, others N⁺ + N̄₀; Levinson's first step is lossless, DH 1301 = 1437 − 2·68), but in the hyperbolic half-strip the SMT is Littlewood plus a pointwise arithmetic inequality, its only analytic input (the log-derivative lemma) bounds zero counts from below, and at T = 10⁴ its best bound is ≤ −112 against true values −9 to 770 for ζ and DH alike; disc and half-plane theories are RH-blind (zero moves change N(r, 0) by ≤ 0.0231; ζ's half-plane characteristic is 0.030); the Cauchy law of Z′/Z is known (2609.15862) [ExploreValue](Riemann/Deep/ExploreValue/REPORT.md)

RH relevance: N8 (Levinson–Conrey branch: barrier, unconditional); N13 (Speiser = Levinson–Montgomery, no ζ/DH separation beyond N⁺); N9 (no transfer).
