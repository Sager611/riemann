SOUND: Re s > B_* − 1/160000 = 0.87446205769353… is zero-free for every finite-order Hecke L-function over every number field, so also for ζ(s) and Dirichlet L-functions, conditional on Liu's two unrefereed 2026 preprints (Hexagon 2610.00655v1, which supplies β_*(F) ≤ B_* = 0.874468307693530 and the comparison machinery; arXiv 2610.26835v2, which supplies theta, M1/M2 including the prime-free M2 at arbitrary lengths, detector, count and Mellin continuation). No false step in the new parts; the endpoint is already beaten by Liu's own later Hexagon claims.

# Chang–Wang, arXiv 2611.00924v1: audit (node N7, 2026-11-05)

Chang and Wang take Liu's Hexagon boundary B_* and lower it by 1/160000. Their new tool is a stopping rule. When primes are removed from a returned product and all of them are gone, they bound the rest with the prime-free moment theorem, which allows any lengths. Liu's version had required a negative length defect at that point. With this rule a slightly longer extra prime product is allowed near the critical length, which gains 1/20000 in the first mean. Every finite inequality in the paper checks out, mostly by interval proofs (verify.py: 77 checks, 0 failed). The analytic steps follow Liu's framework, and Liu proves the same norm theorem himself in his v2 of 3 Nov. The result stands or falls with Liu's unrefereed inputs.

Sources, read as TeX: CW `main.tex` (e-print, 668 lines; "l."), [H] Liu, Hexagon 2610.00655v1 `number_field_new_v1.tex` (29 Oct; "H l."), its v2 (3 Nov 03:16 UTC) and v3 (3 Nov 17:40 UTC; "H2 l."), [Liu] arXiv 2610.26835v2 `number_fields_v4.tex` ("Liu l."). Kazhdan–Patterson and Milne were not downloaded (peer-reviewed or textbook).

## 1. Inputs and their status

| # | Input as used (CW line) | Source (locator) | Status |
|---|---|---|---|
| I1 | General-field arithmetic: S, generators, Θ, theta coefficients a₀, twisted multiplicativity, reciprocity, full local Gauss identities, orbit Poisson (l.75–97) | [Liu] Part III, Thm vi:prop:transfer (Liu l.11906), Part II §§ Poisson/compact dual (Liu l.3476ff.) | Unrefereed arXiv. Our audit Riemann/Deep/LiuHecke: no refutation |
| I2 | Metaplectic normalization (l.97) | Kazhdan–Patterson, Publ. IHES 59 (1984), 62 (1985) | Peer-reviewed |
| I3 | M1 and M2: z = 0 at any bounded lengths (principal rows excluded); z > 0 with n₁+n₂+6κz ≤ M for 13/18 ≤ κ ≤ 1 under β_* ≤ (1+κ)/2 (l.97, 298–301) | [Liu] Thm v:thm:M2-assembled (Liu l.7575–7583), Thm v:thm:marked-inverse-new (Liu l.5657) | Unrefereed. κ below 3/4 audited in Riemann/Deep/Kappa (no refutation). The z = 0 clause is marked "proved in Part II" in our LiuHecke audit |
| I4 | Detector, class count R_k, same-object contour comparison, both AFE arms, direct estimate B(ℓ) = 11/12 − ℓ/4 (l.99, 407–412, 554–558) | [Liu] §14, Rem. v:rem:refined-admissible (Liu l.10112); [H] Prop jdc:count (H l.3472), Lemmas jwa:profiles/contours (H l.3142, 3212), Lemma jca:additive (H l.3683) | Unrefereed |
| I5 | Prior bound β_*(F) ≤ B_* < 0.8745. It licenses the cap κ = 749/1000 and the count parameter k_* = 2B_* − 1 (l.101–104) | [H] Thm jcg:theorem (H l.4286), Thm nf:main (H l.145) | Unrefereed Hexagon preprint, revised twice on 3 Nov, not audited by us |
| I6 | Theta/Mellin identity with H_η holomorphic and nonzero on Re s > 437/500; principal savings l_y/20, h/600, ω_P (l.106–110, 615–618) | [Liu] eq. v-principal-buffered-savings (Liu l.9442) and the continuation section | Unrefereed |
| I7 | [H] machinery reused inside the new norm: valuation classes and completed lengths (Thm zpb:theorem, Cor. cnr:homogeneous), small-formal estimate (Prop jtf:terminal, H l.2412), benchmark (Lemma jdr:benchmark, H l.2543), two-transform return (Prop jdr:return, H l.2677), amplitude selection (Lemma jem:amplitude, H l.2089), envelope (Lemma jca:taper, H l.3786; Prop jca:realize, H l.3876), small and lowest amplitudes (Prop jsb:comparison, Lemma jfm:floor) | [H] | Unrefereed |
| I8 | Abelian factorizations for the descent (l.639–656) | Milne, CFT notes | Standard |

**Conditional on:** I1, I3–I7, all of them Liu's, over every F = K(ζ₁₂). I5 is a full theorem of [H] taken as a premise. It is not circular, since CW only bootstrap from it (l.104).

**Applies to ζ:** yes. Take K = Q and F = Q(ζ₁₂). ζ is a factor of ζ_F (Thm 1.1, l.66–68; descent l.639–656). The result is not limited to Hecke L-functions of special fields.

## 2. Findings on the new parts

1. **Endpoint arithmetic (l.41–64, 439–452).** r_*, ℓ_*, b_*, k_*, B_* agree with [H]'s 20 printed digits. The derivative definition of b_* (l.56) equals [H]'s closed form. B = 0.87446205769353045955. The root box, the Taylor coefficients at .45057, the three derivative enclosures (rigorous intervals), k(r) ≡ 5/6 − ℓ(r)/2, n_b(r,k(r)) ≡ 0 and 3D²{ℓn_e+n_0} ≡ −1944r(6r−5)²p(r) are all exact. A_Q = 2329.39 > 2200. (verify.py A1–A15.)
2. **Exhaustion rule (Lemma 4.2, l.285–303).** The removal arithmetic b ≤ min{z_rem, D₊/K+ω}, κb ≤ D₊/6+κω is exact, because κ/K = 1/6. At exhaustion with a positive defect the lemma calls M2 with z = 0. [Liu] states that case for arbitrary bounded plain lengths, with only principal rows excluded (Liu l.7580–7581), and the children lie outside Θ ⊃ {1}. The call is to a preceding input, so it is not circular. [H] v1 had instead forced D_free ≤ −η < 0 at exhaustion, through the coupled condition W < 7(1−2r)/54 (H l.3010–3015, 3098). Dropping that is the whole new idea. Liu makes the same move in [H] v2, Lemma en:exhaustion (H2 l.7105–7128: "No condition D ≤ pz′ is required"). No flaw.
3. **Local normalization (l.244–263, 311–321).** The D_free identity (l.261–263), the R_loc simplification and R_loc + D₀/6 (l.313–319) are exact identities (C2). Nonpositivity needs only nonnegative allocation variables, η_K ≥ 0, z_C ≥ 0 and w_f ≤ μ_com. It is [H] Cor. cnr:homogeneous with 9/2 replaced by K, and it equals H2 l.7270.
4. **Other branches of Thm 4.1 (l.264–267, 323–345).** M* ≥ 1−2r_o−6W and m* ≥ 1−2r_o−7W hold. The small-formal bound is −1/2+3d₀/K+7/200 ≤ −93/200+90/749 = −0.3448 < −1/3, and it only has to be negative. It rests on [H] Prop jtf:terminal's identity, which we checked exactly and which does not use [H]'s η > 0. The origin bound −c_Kd₀+5W < 0 and the exceptional bound −c_Kd₀+4W < 0 hold because c_K = 1747/2247 > 5/7. Restoration: 3/2−2m+K(z+w) ≤ 5/2−4m+Kw < 1. As a cross-check, CW's domain (l.231–234) lies inside the domain of Liu's later Thm en:variable-rectangle-norm (H2 l.7049–7064) with short side t = 1/4 and p = 6κ = 4.494 ∈ [112/25, 9/2].
5. **Cap κ = 749/1000 (l.104, 220): minor point CW do not address.** [Liu]'s refined comparison is stated for κ ≥ 7497/10000 (Liu l.10115), and CW use 0.749 on the numerator side without comment. Both uses still hold. The crude supply is (21/100)(6κ−1/8) − 11/12 = 8.2·10⁻⁴ > 0; Liu's own threshold is κ > 2263/3024 = 0.74835 (H2 Rem. kb:scalar-corr(b), l.5901–5915). Also 6κℓ > 3/4. (C17.)
6. **Gain ν = 1/20000 (Prop 5.1, l.347–404).** The margins 40523/6741000 and 15343/500000 are exact; their minima sit at m = .43 and m = .41. The cover margin is 0.0042 > 1/250 (l.367). z_K < 41/1000, W_new < 1/40, and a(z_K) ≥ 0.13175 > 13/100 (minimum at δ = 7/20, x = 2/5, m = .41). The F′(s) formula (l.387–390) is exact and F′ > 0, so G_K ≥ G_0 since 1/K = 500/2247 > 2/9. G_new < 1/100, and the net gain is 13/200000 − 9·10⁻⁶ = 5.6·10⁻⁵ > 1/20000. The bookkeeping matches [H]: the fourth moment saves 2gz₀+4γ (H l.3713), the initial-product share sits in E through θ₀(1−2m)/2 = gz₀/2, and a(z₀(m)) ≥ a_□ = a(1/18) (checked). So G₀ dominates the source envelope, and only the +1/1000 width counts as new (l.391–393). No flaw.
7. **Scalar comparison (§6, l.406–550).** Every printed inequality holds. Interval branch-and-bound at the real parameter ranges (b ∈ [.09175,.09185], ℓ ∈ [ℓ_*,.16883], k = k_*) confirms:
   - |E_b| < 1/5 and 5/16 ≤ E_ℓ ≤ 5/4 for k ∈ [2/3,3/4];
   - Lemma 6.1: Q_x > 0 and the denominator is < 0, so Q ≤ 0 with equality only at (r_*,1/2); S ≥ 71/300;
   - complementary classes E < −53/400000 (certified max −1.3250·10⁻⁴);
   - full box E < 779/400000; short sums E + hΔ_s < −1/100 (certified −0.01001).
   The closed forms of E(δ,1/2) and E(δ,2/5), the y-expansions and the concavity data are exact. One cosmetic slip: the −69911/6480000 in (6.6) (l.547) equals 1/500 − .79·1049/64800, so it uses 1/500 instead of the sharper 779/400000. Harmless.
8. **Global comparison (Prop 7.1, l.563–604).**
   - Q_ℓ < 593/500 (certified max 1.18599), so Q(ℓ_*+1/40000) < 593/20000000. The true maximum is 2.835·10⁻⁵, at the contact class, where Q_ℓ = 1.1338.
   - Band exponent ≤ −9.85·10⁻⁶; the true value is −1.16·10⁻⁵. λ = 0.577491, m₀ = 0.422509, λ < 11/19.
   - The envelope claims of l.576–578 (slopes 7/225 and 49/432, excess ≥ 7/150000 outside [.41,.43], 7δ/2160 above 9/20) are stated without proof. They are true for [H]'s envelope (Lemma jca:taper, H l.3786–3860): certified excess ≥ 4.69·10⁻⁵, true minimum 8.97·10⁻⁵ at m = 3/8, δ = 7/20, x = 1/2. [H]'s identities for Γ_L − Γ_N and for the outer reserve are exact.
   - The coarse-box conditions (l.589–594) hold, as do ℓ/(h+1/2000) > 21/100 and ℓ/h < 2/9.
   - The case d ≤ h − 1/80 (l.582–584) is correct but redundant, since the exponent increases with d.
   (D1–D12.)
9. **ℓ outside [H]'s certified window.** ℓ = ℓ_* + 1/40000 = 0.1688184 lies outside [H]'s verified interval .16876 ≤ ℓ ≤ .16880 (H l.4156). So CW cannot cite [H]'s complementary, short-sum and geometric margins. They re-prove them up to ℓ = .16883 (l.487–550, 588–596), and those proofs check.
10. **Continuation and descent (§§8–9, l.606–656).** Same structure as [Liu]. The saving s = min(r_np, r_P)/32 is more conservative than Liu's /16. The descent through E = MF and the factorizations of ζ_E over F and over M is standard.
11. **False steps: none found.** The weakest links are the inputs I1 and I3–I7, and the sketch-level reuse of [H]'s return (l.278–283) at κ = 0.749 without the coupled margin. Liu writes that extension out himself in [H] v2, Part II.

## 3. N7 context (orchestrator: TREE.md and R279 need fixing)
- CW's B_* = 0.874468307693530 is Liu's Hexagon 2610.00655v1 (29 Oct), not arXiv 2610.26835 (0.874886647040). TREE.md N7, R279 and the brief conflate the two.
- Later unrefereed claims by Liu, all on Hexagon, none citing CW:
  - 2610.00655v2 (3 Nov): B′ = 0.874221250 without exhaustion; B_ex = 0.874072848 with the same exhaustion idea as CW;
  - 2610.00655v3 (3 Nov): B_tw = 0.874058504;
  - **2611.00082v1 (4 Nov): B₂ = (777 − √1081)/852 = 0.873381967**, from an ordinary product estimate with condition n₁+n₂+2z ≤ M.
  The arXiv math.NT listing for November 2026 (140 entries) has no other half-plane paper.
- CW's endpoint (arXiv, 31 Oct) is therefore not a record. As mathematics it is consistent.
- Suggested N7 text: "Chang–Wang 2611.00924: B = 0.874462058, sound given Liu's inputs ([audit](Riemann/Deep/ChangWang/REPORT.md)); superseded by Liu, Hexagon 2610.00655v2/v3 and 2611.00082 (B₂ = 0.873381967, unaudited)."

## 4. Reproduction
`python3 verify.py` (stdlib + mpmath): 77 checks, 0 failed, about 6 min, 27 MB, single thread. Exact identities use Fractions. Range claims use mpmath.iv branch-and-bound, which is rigorous up to mpmath's outward rounding. The "info" lines are non-rigorous grid values.

Sources: [arXiv 2611.00924](https://arxiv.org/abs/2611.00924) · [Hexagon 2610.00655](https://hexagonmath.org/2610.00655v1) (v2, v3) · [arXiv 2610.26835v2](https://arxiv.org/abs/2610.26835v2) · [Hexagon 2611.00082](https://hexagonmath.org/2611.00082v1)
