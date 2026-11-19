# nbaudit FULLPROOF: the N5 chain with source line numbers

Sources (TeX in .scratch/nbaudit/tex, deleted after): BD = Báez-Duarte math/0202141; Bu = Burnol math/0103058; BCF = 1211.5191; K = 2611.03924 riemann.tex;
Z2 = Zenodo 21898513 manuscript/rh_nyman_beurling_paper.tex; NB = Riemann/Deep/NBRate (REPORT.md, code/s2_rig.py = s2).

## 1. Definition (holds)
- BD l.44, 56-66: H = L_2(0,∞), ρ(x) = x − [x], ρ_a(x) = ρ(1/(ax)), χ = indicator of (0,1]; RH ⟺ χ ∈ closure span{ρ_a : a ∈ ℕ}.
  So d_N² = inf_c ‖χ − Σ_{k≤N} c_k ρ_k‖²_{L²(0,∞)}, c free (no constraint): NB's V1. ✓
- Bu l.168-176, 186-190: same K = L²(]0,∞[), χ = 1_(0,1], D(λ) = dist(χ, span{ρ(θ/t): λ ≤ θ ≤ 1}). Since {1/k : k ≤ N} ⊂ [1/N, 1], d_N ≥ D(1/N).
- BCF l.64-83 quote BBLS's liminf d_N² log N ≥ Σ 1/|ρ|² and Burnol's Σ m_ρ²/|ρ|², and Σ_{Re ρ=1/2} m(ρ)/|ρ|² = 2 + γ − log 4π under RH (l.79).
- The two comparison sources compute the CONSTRAINED distance (Σ c_k/k = 0, approximants supported on (0,1]) = NB's V2 ≥ V1:
  K l.749-751 (G_k, k = 2..n, H²(C_{1/2}) = L²(0,1) by l.161); Z2 l.44-47, 101: ‖1 − Σ_{n=2}^N c_n ρ_{1/n}‖_{L²(0,1)}, ρ_a = {a/x} − a{1/x}.

## 2. Gram, target, norm (holds)
- G_jk = ⟨ρ_j, ρ_k⟩ = ∫_0^∞ {1/(jx)}{1/(kx)} dx = ∫_0^∞ {u/j}{u/k} du/u² (u = 1/x). Homogeneity: x → x/d gives G_{dj,dk} = G_jk/d.
- b_k = ⟨χ, ρ_k⟩ = ∫_1^∞ {u/k} du/u² = (1/k)[log k + Σ_{j≥1}(log(1+1/j) − 1/(j+1))] = (log k + 1 − γ)/k. ‖χ‖² = 1.
- Vasyunin (coprime h,k; quoted NB REPORT l.6-7, s2 l.84): G_hk = [(log 2π − γ)(h+k) + (k−h) log(h/k) − π(V(h,k)+V(k,h))]/(2hk),
  V(h,k) = Σ_{m=1}^{k−1} {mh/k} cot(πm/k). Not in any TeX we hold (K l.789 cites Vasyunin only); checked independently in §5 (direct exact integration).
- G_11 = ∫_0^1 1 + ∫_1^∞ {u}²/u² = log 2π − γ ✓ (formula with V ≡ 0).
- s2's sum: pairing m ↔ k−m, V(h,k) = (1/k)Σ_{m≤(k−1)/2} (2(mh mod k) − k) cot(πm/k) (cot(π/2) = 0 drops m = k/2) ✓ (s2 l.28-55);
  V(k,h) = V(k mod h, h) ✓ (s2 l.83); every (i,j) = (d·a, d·b), d = gcd, gcd(a,b) = 1, visited exactly once (s2 l.73-92) ✓.

## 3. Linear algebra (holds)
- Q(c) = ‖χ − Σ c_kρ_k‖² = 1 − 2bᵀc + cᵀGc. G ≻ 0 ⇒ minimiser c* = G⁻¹b, d_N² = 1 − bᵀG_N⁻¹b, and for ANY c: Q(c) − d_N² = rᵀG_N⁻¹r, r = G_Nc − b_N.
  Hence d_N² ∈ [Q(c) − ‖r‖²/λ_min(G_N), Q(c)] (s2 l.4-5, 177-180) ✓.
- λ_min(G_N) ≥ λ_min(G_1280): G_N is a principal submatrix (Cauchy interlacing) ✓. λ_min(G) ≥ μ − ‖E‖_2 ≥ μ − ‖E‖_F for G − μI = L̂L̂ᵀ + E (Weyl) ✓.
- c for each N: Cholesky of G_1280 truncated (leading block of L = Cholesky of leading block) ✓ (s2 l.132-133, 177). Any c is valid; c only affects width.

## 4. Code audit of s2_rig.py (holds)
- Cot sums (l.28-55): C_m = round(2¹⁴⁰ cot(πm/k)) from 220-bit Arb balls; |C_m 2⁻¹⁴⁰ − cot| ≤ rad + 2⁻¹⁴¹ = e_m; R = 2(ma mod k) − k, |R| < k < 2¹¹;
  26-bit limbs; R @ climb exact (|partial sums| < 2¹¹·2²⁶·640 < 2⁴⁷ < 2⁵³, l.46-47); V = S/(k 2¹⁴⁰) ± h·max e_m (l.41) ✓.
- g (l.73-84) in Arb with V as balls; Gt = round(2¹²⁰ mid), ε_G = max(rad·2¹²⁰ + ½)/d + ½ ulps over d (l.85-92) ✓; limbs < 2²⁰ ✓; b likewise (l.98-103) ✓; Part 2 = l.57-104.
- certify (l.140-165): ĉ = round(2⁶⁰c); w = G̃ĉ exact (20-bit limbs, < 2⁵¹; int64 sums < 2⁵⁴); Q exact rational; |Q − Q̃| ≤ ε_G S1² + 2ε_b S1 (l.160);
  |r_j − r̃_j| ≤ ε_G S1 + ε_b (l.163); r2 ≥ ‖r‖² with (1+10⁻¹²), (1+10⁻⁹) margins covering float conversion, summation (≤ 1.4·10⁻¹³) and the division r2/λ (l.164, 179) ✓.
- lo = Q − errQ − r2/λ, hi = Q + errQ exact Fractions (l.180); sign by Arb: (ball·log N − C) < 0 is True only if certain (l.182-185) ✓. N = 1..NMAX, no gaps ✓.
- λ certificate (l.110-127): E := G − μI − L̂L̂ᵀ exact; L̂L̂ᵀ = Σ_ab P_ab 2^{17(a+b)−100} exactly (P_ab integer < 2⁴⁷; Z = Z0 + 2¹⁷Z1 + 2³⁴Z2 exact for negative Z too);
  |E − Ê| ≤ ε_G + u|G̃| + γ₈ Σ_ab|P_ab|2^… + u(|G_f|+μ) + u(|G_f|+μ+|S|) ≤ ε_G + (3u+O(u²))(|G_f|+μ) + (9u+O(u²))‖Y_i‖‖Y_j‖ < U (12u) ✓;
  Σ_ab |P_ab,ij| 2^… ≤ Σ_k Y_ik Y_jk ≤ ‖Y_i‖‖Y_j‖ (Cauchy–Schwarz; Y exact: integers < 2⁵² times 2⁻⁵⁰); λ_min(G) ≥ μ − ‖E‖_F (Weyl) ✓. Reviewer: CORRECT.
- Width (not in the code's printed summary): "max width hi-lo over N: 2.22e-16" (l.208) is the float width of nextafter-rounded bounds; the exact width
  was printed only at 18 N (l.193-195). Exact tracking (a2/c1, 3 added lines): max 3.138·10⁻²³ at N = 1068 (NMAX 1280), 4.108·10⁻²⁴ at N = 622 (NMAX 640).

## 5. Independent recomputation (holds)
- b2 proves the Gram entries pairwise without Vasyunin: on [A,B] with ⌊u/h⌋ = P, ⌊u/k⌋ = Q: ∫ (u/h−P)(u/k−Q)u⁻² du = (B−A)/(hk) − (Q/h+P/k) log(B/A) + PQ(1/A − 1/B);
  tail Σ_{j≥J} ∫_0^L f(v)(jL+v)⁻² dv = Σ_n (n+1)(−1)ⁿ m_n ζ(n+2, J)/L^{n+2}, m_n = ∫_0^L f vⁿ dv exact rationals (f ≤ 1: terms ≤ J^{−n}). 40 digits.
- b1 (Vasyunin, 60 digits, one Cholesky, d_N² = 1 − Σ_{i≤N} y_i² for all N ≤ 100) lands inside every certified interval; b3 (double) reproduces all 16 runs.

## 6. Constant and comparisons
- C = 2 + γ − log 4π = 0.0461914179 (BCF l.79: Σ_{Re ρ=½} m(ρ)/|ρ|² under RH). Burnol (Bu l.202-205): liminf D(λ)² log(1/λ) ≥ Σ m_ρ²/|ρ|² unconditionally (only on-line zeros contribute);
  d_N ≥ D(1/N) since {1/k : k ≤ N} ⊂ [1/N, 1], so liminf d_N² log N ≥ C under RH; BCF l.75-83: conjectured d_N² ~ C/log N. Sign changes at finite N contradict nothing (liminf).
- V2 (constrained, Σ c_k/k = 0 ⟺ support in (0,1]): Gram A = TᵀGT, A_jk = G_jk − G_j1/k − G_1k/j + G_11/(jk), b′_k = b_k − b_1/k = log k/k (= |⟨E, G_k⟩| of K l.784).
  Check: A_22 = log 2/4 = K's (1/4)Σ_m (m mod 2)²/(m(m+1)) ✓. Values: REPORT comparison table.
