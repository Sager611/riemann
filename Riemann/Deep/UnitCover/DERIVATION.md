# Covering every unit group (ℤ/qℤ)^×, q ≤ Q, by subset products of O(log Q) small primes — derivation

Every step and every numeric constant used by `cert.py` (section numbers match its comments). Notation: Q ≥ 3; q a modulus,
2 ≤ q ≤ Q; G = G_q = (ℤ/qℤ)^×, H = |G| = φ(q); Ĝ its character group (Ĝ ≅ G); χ, ψ characters, 1 the principal one;
G[d] = {g : g^d = 1}; ω(q) the number of prime factors; τ_D(q) the number of ordered factorizations q = n₁⋯n_D (τ = τ₂);
L := max(log Q, log 9). "P covers q" means {∏_{p∈F} p mod q : F ⊆ P, p ∤ q for all p ∈ F} = G_q.

## 0. Inputs

- **(H)** Bebel (QuasiRiemannTracker, Lean-kernel-checked): ζ(s) and every Dirichlet L(s, χ) are zero-free on Re s > B₂ + 2.2·10⁻²⁵⁴,
  B₂ = (777 − √1081)/852 = 0.87338196661320751993…. We use **B = 0.8733819666132076** (certified B > B₂ + 10⁻¹⁷),
  b = 1 − B = 0.1266180333867924, **E := 1/(1 − B) = 7.89776916646**. Consequence used: for primitive χ every nontrivial zero has
  b ≤ β ≤ B (functional equation), no exceptional zero.
- **explprog** (Riemann/Deep/ExplProg, DERIVATION §1–§3, imported unchanged by cert.py): Lemma EF (explicit formula for weights
  supported in [t₀, t₁] ⊂ (0, ∞)), Lemma P (pairing ρ ↦ 1 − ρ̄), Lemma Z (zero sums from the refereed zero count of
  Bennett–Martin–O'Bryant–Rechnitzer, Math. Comp. 90 (2021) 1455–1482, Thm 1.1), Lemma Zζ (Σ_ρ 1/γ² ≤ 0.046371186 over ζ zeros,
  from Σ Re 1/ρ = 1 + γ_E/2 − ½ log 4π), Lemma K1 (kernel decay), (1.1) (trivial zeros), the triangle kernel of §3.1:
  g on [−7, 0] rising from 0 at −7 to 1 at −21/10 and falling to 0 at 0, w(t) = g(log t)t^{−B}, supp w = [e^{−7}, 1],
  **W(1) = 2.424245725, W(0) = 115.9708945, ‖w‖_∞ ≤ 38.854381**, c_∞ = z_a/(2W(1)) = 0.453311.
- **(RS62)** Rosser–Schoenfeld 1962, Thm 9: θ(x) < 1.01624x (x > 0).
- Elementary: φ(q) ≤ q − 1 (q ≥ 2); Σ_{q≤Q} τ_D(q) ≤ Q(1 + log Q)^{D−1}; D^{ω(q)} ≤ τ_D(q).
- Not used: GRH, any zero-density theorem, any unrefereed zero count, Linnik's theorem.

## 1. The random construction (Theorem 1)

Fix y (chosen in §2) and the probability law ν on primes p ∈ (e^{−7}y, y]: **ν(p) = log p · w(p/y)/M**, M := Σ_p log p · w(p/y).
Draw p₁, p₂, …, p_N independently from ν; F_j := σ(p₁, …, p_j). For a modulus q define g_j = g_j(q) ∈ G_q by
g_j = p_j mod q, except **g_j = 1 ("no-op") if p_j | q or p_j ∈ {p₁, …, p_{j−1}}**. Put S₀ = {1}, S_j = S_{j−1} ∪ g_j S_{j−1}
(so S_j = S_{j−1} at a no-op), T_j = G ∖ S_j. Let P := {p₁, …, p_N} (distinct values), |P| ≤ N.

**Lemma 1.1.** S_N(q) ⊆ {subset products of P ∖ {p : p | q}} (mod q). *Proof.* An element of S_N is ∏_{j∈F} g_j for some F ⊆ {1..N};
dropping no-ops (g_j = 1) leaves indices whose primes are pairwise distinct (a repeat is a no-op) and coprime to q. ∎
So if S_N(q) = G_q for every 3 ≤ q ≤ Q, then P covers every 2 ≤ q ≤ Q (q = 2: G trivial). Hughes's form (q coprime to ∏P) is the
special case where no p ∈ P divides q.

## 2. Lemma D: character means of the law (cert.py sec 2)

**Lemma D.** Let 3 ≤ q ≤ Q, ψ a nonprincipal character mod q, 1 ≤ j ≤ N. Then Re E[ψ(g_j) | F_{j−1}] ≤ δ(y, Q, N) with

  δ := num/den,  num = (y^B Z(L) + y^{1/2} Z_L(L))/2 + W(0)/(1 − e^{−2Y}) + PP(y) + 2‖w‖_∞((N − 1) log y + L),
  den = yW(1) − (y^B z₁ + y^{1/2} z₁^L) − W(0)/(1 − e^{−2Y}) − PP(y) (> 0),

Y = log y − 7, (n₀, k₁, k₂), (m₀, k₁^L, k₂^L) from Lemma K1 at Y, Z = Z_G(L) and Z_L = Z_{G_L}(L) from Lemma Z, z₁ = 0.0231856·k₂ and
z₁^L = 0.0231856·k₂^L (Lemma Zζ: 0.0231856 = 0.046371186/2), PP(y) = ‖w‖_∞·1.01624(√y + y^{1/3} log y/log 2).

*Proof.* Let ψ* (conductor q* ≥ 3, q* | q) induce ψ; ψ(p) = ψ*(p) for p ∤ q. Let A_j = {p₁..p_{j−1}} ∪ {p : p | q}. Given F_{j−1}, g_j = 1
on A_j and g_j = p mod q otherwise, so E[ψ(g_j)|F_{j−1}] = Σ_{p∉A_j} ν(p)ψ*(p) + ν(A_j), and since |ψ*(p)| ≤ 1,
Re E[ψ(g_j)|F_{j−1}] ≤ |Σ_p ν(p)ψ*(p)| + 2ν(A_j). (i) ν(A_j) ≤ (‖w‖_∞/M)((j − 1) log y + Σ_{p|q} log p) ≤ ‖w‖_∞((N − 1) log y + L)/M,
since ν(p) ≤ ‖w‖_∞ log y/M and Σ_{p|q} log p ≤ log q ≤ L. (ii) M·Σ_p ν(p)ψ*(p) = Σ_n Λ(n)ψ*(n)w(n/y) − Σ_{p^k, k≥2} log p ψ*(p^k)w(p^k/y);
the second sum is ≤ ‖w‖_∞(ψ(y) − θ(y)) ≤ ‖w‖_∞·1.01624(√y + y^{1/3} log y/log 2) = PP(y) in modulus (ψ − θ = Σ_{k≥2} θ(y^{1/k}), at most
log y/log 2 terms, RS62). By Lemma EF for ψ* (no main term, ye^{−7} > 1) and (1.1),
|Σ_n Λψ*w(n/y)| ≤ Σ_ρ |y^ρ W(ρ)| + W(0)/(1 − e^{−2Y}), and by Lemmas P, K1, Z (conductor q* ≤ e^L, L ≥ log 9)
Σ_ρ |y^ρW(ρ)| ≤ (y^B Z(L) + y^{1/2}Z_L(L))/2. (iii) M ≥ Σ_n Λ(n)w(n/y) − PP(y) and, by Lemma EF for ζ, (1.1) and Lemma Zζ,
Σ_n Λ(n)w(n/y) ≥ yW(1) − (y^B z₁ + y^{1/2}z₁^L) − W(0)/(1 − e^{−2Y}). Combine (i)–(iii). ∎

**Monotonicity.** Write δ = n(y)/d(y) with n = num/y^B, d = den/y^B:
n = Z/2 + y^{1/2−B}Z_L/2 + (W(0)/(1 − e^{−2Y}) + PP)/y^B + 2‖w‖_∞((N − 1) log y + L)/y^B,  d = y^b W(1) − z₁ − y^{1/2−B}z₁^L − (W(0)/(1 − e^{−2Y}) + PP)/y^B.
The k's of Lemma K1 decrease in Y, so G, G_L decrease pointwise and Lemma Z applied to the larger envelope stays valid; y^{1/2−B},
y^{1/3−B} log y (log y > 1/(B − 1/3)) and log y/y^B (log y > 1/B) decrease. Hence n is non-increasing and d non-decreasing in y, so
δ is non-increasing in y; δ is non-decreasing in N. No monotonicity in L is needed: for every Q with log Q ≤ L₁ all conductors are
≤ e^{L₁}, so Lemma Z evaluated at L₁ applies directly, and Σ_{p|q} log p ≤ L₁; den does not depend on L. Hence on a cell
L ∈ [L₀, L₁] with y = (c'L)^E ≥ (c'L₀)^E: δ ≤ num(y₀, L₁, N)/den(y₀), y₀ = (c'L₀)^E, with den(y₀) > 0 certified at that corner
(cert.py delta_direct with L0; reviewer remark, 6 Nov).

**2.3 Uniform version (cert.py delta_cells, delta_tail).** Put **y := (c'L)^E**. Claim: with **c' = 820, δ₀ = 1/1000**, and
N − 1 ≤ 4L + 200: **δ ≤ 1/1000 for every Q ≥ 3**. Proof: (a) 664 cells [L₀, L₀ + max(1/100, L₀/200)] cover [log 9, 60.275]; on each
the certified cell bound is ≤ 0.0009837. (b) For L ≥ L_m = 60.275: δ ≤ δ₀ ⟺ num + δ₀(y^Bz₁ + y^{1/2}z₁^L + triv + PP) ≤ δ₀yW(1);
dividing by y^B L (y^b = c'L) it suffices that
  Z(L)/(2L) + y^{1/2−B}Z_L(L)/(2L) + (1 + δ₀)[W(0)/((1 − e^{−2Y})y^B) + ‖w‖_∞1.01624(y^{1/2−B} + y^{1/3−B} log y/log 2)]/L
  + 2‖w‖_∞((4 + 200/L) log y + 1)/y^B + δ₀(z₁ + y^{1/2−B}z₁^L)/L ≤ δ₀c'W(1).
Every term is non-increasing in L ≥ L_m when the k's are taken at Y_m = Y(L_m) (k's decrease in Y; Z(L)/L by explprog's majorant
Z_over_L_majorant; y^{1/3−B} log y and log y/y^B decrease for log y > 2 and > 1/B; (4 + 200/L) decreases), so evaluating at L_m
suffices: certified margin > 0.344. The same with **c' = 820000, δ₀ = 10⁻⁶** (worst cell ≤ 9.826·10⁻⁷, tail margin > 0.345).
(The least c' with δ ≤ 10⁻³ is about 0.80/δ₀ near L = 3 and tends to c_∞/δ₀ = 453.3; the table of §5 uses the least c' per row.)

## 3. Lemma C: covering a fixed modulus (cert.py sec 3)

Fix q, H = φ(q) ≥ 2. Assume Re E[ψ(g_j)|F_{j−1}] ≤ δ for all nonprincipal ψ and all j ≤ N, with **0 ≤ δ < 1** (Lemma D gives
δ > 0 for the finite horizon N; δ ≥ 0 keeps every weight π + δ(1 − π) and β nonnegative, which the over-count (3.1) and the hole
stage need — reviewer remark, 6 Nov).

**3.1 Fourier facts.** μ_j := law of ∏_{i≤j} g_i^{ε_i} with independent fair bits ε_i (given g); supp μ_j = S_j and
μ̂_j(χ) := Σ_x μ_j(x)χ(x) = ∏_{i≤j}(1 + χ(g_i))/2. Put c_χ(g) := |(1 + χ(g))/2|² = (2 + χ(g) + χ̄(g))/4 ∈ [0, 1] and the collision
**C_r := Σ_{χ≠1}|μ̂_r(χ)|² = H‖μ_r‖₂² − 1** (Parseval), so |μ̂_r(χ)|² = ∏_{j≤r} c_χ(g_j). If |T_r| = h then (Cauchy–Schwarz on the support)
‖μ_r‖₂² ≥ 1/(H − h), so C_r ≥ h/(H − h). Hence, with **t := α/(1 − α)**, |T_r| ≥ αH ⇒ C_r ≥ t.

**3.2 The k-th moment (k = 1, 2, 3, 4).** For a k-tuple a = (χ₁, …, χ_k) of nonprincipal characters let
R(a) = {e ∈ {−1,0,1}^k : ∏ χ_i^{e_i} = 1} and **π(a) := Σ_{e∈R(a)} 2^{−k−|e|}** (|e| = #{i : e_i ≠ 0}). Expanding,
∏_i c_{χ_i}(g) = Σ_{e∈{−1,0,1}^k} 2^{−k−|e|} ∏_i χ_i^{e_i}(g) (weights sum to 1, the sum is real), so
E[∏_i c_{χ_i}(g_j)|F_{j−1}] ≤ π(a) + δ(1 − π(a)). By the tower property (factors in [0, 1], F_j-measurable),
**E[C_r^k] = Σ_a E∏_{j≤r}∏_i c_{χ_i}(g_j) ≤ Σ_a (π(a) + δ(1 − π(a)))^r**, and by Markov **Pr(|T_r| ≥ αH) ≤ E[C_r^k]/t^k.**
(k = 1 is Hughes's collision step, arXiv 2610.04678 l.149–158, with no bad characters.)

**3.3 Closed relation sets.** R(a) is *closed*: R = span_ℤ(R) ∩ {−1,0,1}^k and span_ℤ(R) contains no unit vector (R ⊆ Λ(a) :=
{e ∈ ℤ^k : ∏χ_i^{e_i} = 1}, Λ(a) ∌ unit vectors as χ_i ≠ 1, and span_ℤ(R) ∩ cube ⊆ Λ(a) ∩ cube = R). For a closed R with
lattice span_ℤ(R) of rank ρ and elementary divisors d₁ | … | d_ρ, the tuples with ∏χ_i^{e_i} = 1 for all e ∈ R form the kernel of
Ĝ^k → Ĝ^R, of size **H^{k−ρ}∏_i|G[d_i]|** (Smith normal form; |Ĝ[d]| = |G[d]|). Hence
  Σ_a (π(a) + δ(1 − π(a)))^r ≤ Σ_{R closed} H^{k−ρ_R}∏_i|G[d_i]| · (π_R + δ(1 − π_R))^r.                                    (3.1)
*Enumeration* (cert.py closed_sets): breadth-first from R₀ = {0}; from each found R and each cube vector v ∉ R form
span_ℤ(R ∪ {v}) (Hermite normal form), discard it if it contains a unit vector, else add its closure. Every closed R is reached
(add its elements one at a time; each intermediate closure lies inside span_ℤ(R) ∩ cube = R). Elementary divisors from the
determinantal divisors (gcd of j×j minors). Result: **1, 4, 48, 2340 closed sets for k = 1, 2, 3, 4** (57 classes for k = 4 by
(π, k − ρ, #{d_i > 1}, ∏d_i); max π = 1/2, attained only by the real "all equal" tuples).
*Independent checks*: enumcheck.py recomputes the closed sets for k = 2, 3 by brute force over all symmetric subsets of the cube
(4 and 48, identical) and meets only enumerated sets for k = 4 among 32403 random and structured tuples in 12 groups with element
orders up to 840 (691 distinct sets met); classcheck.py verifies, on all tuples of 118, 66, 29 small groups (k = 2, 3, 4), that R(a) is enumerated, π
matches and the class count is at most the Smith-form kernel.
*Hand classification* (cross-check, also verified by classcheck.py):
k = 2: π = 1/4 for (H − 1)² generic pairs; 3/8 for the 2(H − |G[2]|) pairs ψ ∈ {χ, χ̄}, χ non-real; 1/2 for the |G[2]| − 1 pairs χ = ψ real.
k = 3 (i ~ j means χ_j ∈ {χ_i, χ̄_i}): no i ~ j and no support-3 relation: 1/8, ≤ (H−1)³; no i ~ j, one ± pair of support-3 relations:
5/32, ≤ 4(H−1)²; no i ~ j, two non-opposite support-3 relations (their sum and difference give χ_i² = 1 for some i): ≤ 1/4,
≤ 12(|G[2]|−1)H; exactly one pair, non-real, χ₃ ∉ {χ₁^{±2}}: 3/16, ≤ 6(H−1)²; same with χ₃ = χ₁^{±2}: ≤ 1/4, ≤ 12(H−1); exactly one
pair, real: 1/4, ≤ 3(|G[2]|−1)(H−1); all three related, χ₁ real: 1/2, ≤ |G[2]|−1; non-real of order 3: 11/32, ≤ 4(|G[3]|−1);
order ≥ 4: 5/16, ≤ 4(H − |G[2]|).

**3.4 Sums over moduli.** For a class with m = k − ρ, n = #{d_i > 1}, D = ∏d_i:
Σ_{3≤q≤Q} φ(q)^m ∏|G_q[d_i]| ≤ **Q^{m+1}/(m+1)** (n = 0, m ≥ 1; φ(q) ≤ q − 1 and Σ_{j≤Q−1} j^m ≤ Q^{m+1}/(m+1)), **≤ Q** (n = m = 0),
**≤ 2^n Q^{m+1}(1 + log Q)^{D−1}** (n ≥ 1), using |G[d]| ≤ 2d^{ω(q)} (for odd p^a the factor is gcd(d, φ(p^a)) ≤ d; for 2^a it is
gcd(2, d)·gcd(2^{a−2}, d) ≤ 2d), d₁^{ω}⋯d_n^{ω} = D^{ω(q)} ≤ τ_D(q) (τ_D(p^a) ≥ D), and Σ_{q≤Q}τ_D(q) ≤ Q(Σ_{n≤Q}1/n)^{D−1} ≤ Q(1 + log Q)^{D−1}.
(For k ≤ 3 also |G[2]| ≤ τ(q) and |G[3]| ≤ 3^{ω(q)} ≤ τ₃(q); verified by classcheck.py on all q ≤ 120.)

**3.5 Hole stage and union.** For j ≥ r, T_{j+1} = T_j ∩ g_{j+1}T_j and E[|T ∩ gT| | F_j] = (1/H)Σ_χ |F_T(χ)|² Re E[χ(g)|F_j]
≤ |T|²/H + δ(H|T| − |T|²)/H (F_T(χ) = Σ_{x∈T}χ(x), Σ_{χ≠1}|F_T|² = H|T| − |T|²; Hughes l.160–170). On {|T_r| < αH}: |T_j| < αH for
j ≥ r, so E[|T_{j+1}| | F_j] ≤ β|T_j| with **β := δ + (1 − δ)α**, and E[|T_{r+s}|·1{|T_r|<αH}] ≤ αHβ^s. So for N = r + s:
  **Σ_{3≤q≤Q} Pr(S_N(q) ≠ G_q) ≤ t^{−k} Σ_classes c·Σ_q(count)·(π + δ(1 − π))^r + α(Q²/2)β^s =: U(Q; k, δ, α, r, s).**   (3.2)
If U < 1, some outcome covers every 3 ≤ q ≤ Q (§1), with |P| ≤ r + s primes, all ≤ y.

## 4. Theorem 1 for every Q ≥ 3 (cert.py sec 4, cover_section)

With ρ_g := 2^{−k}(1 − δ) + δ (generic class), **a₁ := (k + 1)/log(1/ρ_g), b₁ := 2/log(1/β)**, r = ⌈a₁ log Q + a₀⌉, s = ⌈b₁ log Q + b₀⌉:
each class term of (3.2) is ≤ c·coef·Q^e(1 + log Q)^f ρ^{a₁ log Q + a₀} = c·coef·(1 + log Q)^f Q^{−ε}ρ^{a₀}, ε := a₁ log(1/ρ) − e
(ε = 0 for the generic class, whose term is ρ_g^{a₀}/(k + 1); ε > 0 certified for every other class), and
sup_{log Q ≥ log 3}(1 + L)^f e^{−εL} ≤ (f/ε)^f e^{ε−f} (or its value at L = log 3 when f/ε − 1 ≤ log 3); the hole term is ≤ (α/2)β^{b₀}.
Certified sup_Q U < 1 (δ = 1/1000, α = 1/10⁴; y = (820L)^E; N ≤ 4L + 200 holds so Lemma D applies). All values below are upper
ends of the certified balls, rounded up (cert.py ceil_str):

| k | a₁ ≤ | a₀ | b₁ ≤ | b₀ | sup U < | **#P ≤ C log Q + c** | η = C − 1/log 2 ≤ |
|---|---|---|---|---|---|---|---|
| 1 | 2.889557 | 13 | 0.293577 | 0 | 0.62 | 3.1832 log Q + 15 | 1.7405 |
| 2 | 2.168729 | 29 | 0.293577 | 0 | 0.47 | 2.4624 log Q + 31 | 1.0197 |
| 3 | 1.930068 | 42 | 0.293577 | 0 | 0.73 | 2.2237 log Q + 44 | 0.7810 |
| 4 | 1.813106 | 56 | 0.293577 | 0 | 0.55 | **2.1067 log Q + 58** | **0.6640** |
| 4, Q ≥ 10¹⁰⁰ | 1.813106 | 13 | 0.293577 | 0 | 0.54 | 2.1067 log Q + 15 | 0.6640 |

(C ≥ a₁ + b₁ (≤ 2.106682 for k = 4), c = a₀ + b₀ + 2 for the two ceilings.) With **δ = 10⁻⁶, α = 10⁻⁷, y = (820000L)^E**:
k = 3: a₁ ≤ 1.923600, a₀ = 72, b₁ ≤ 0.145771, sup U < 0.66, #P ≤ 2.0694 log Q + 74 (η ≤ 0.6267); k = 4: a₁ ≤ 1.803379, a₀ = 96,
sup U < 0.48, **#P ≤ 1.9492 log Q + 98** (η ≤ 0.5065). As δ, α → 0 (c' → ∞) the k = 4 slope tends to 5/(4 log 2) = 1.8034.

## 5. Table (cert.py table_section): least certified N = r + s per Q, best α on a grid, δ = 1/1000

| Q | log₂ Q | k=1 | k=2 | k=3 | k=4 | **N** | c' (least, 2 digits, rounded up) | max P ≤ (c' log Q)^E ≤ |
|---|---|---|---|---|---|---|---|---|
| 10³ | 10.0 | 28 | 27 | 30 | 33 | **27** | 750 | 10^29.34 |
| 10⁶ | 19.9 | 52 | 44 | 47 | 50 | **44** | 670 | 10^31.33 |
| 10¹² | 39.9 | 97 | 79 | 77 | 82 | **77** | 600 | 10^33.33 |
| 10³⁰ | 99.7 | 231 | 181 | 165 | 164 | **164** | 530 | 10^36.05 |
| 10¹⁰⁰ | 332.2 | 745 | 579 | 525 | 498 | **498** | 490 | 10^39.91 |
| 10¹⁰⁰⁰ | 3321.9 | 7336 | 5676 | 5127 | 4858 | **4858** | 460 | 10^47.59 |

δ = 10⁻⁶: Q = 10¹²: 77 (c' = 600000, 10^57.02); 10¹⁰⁰: 470 (c' = 490000, 10^63.60); 10¹⁰⁰⁰: 4510 (c' = 460000, 10^71.28).
Each row: r, s, α printed in logs/cert_output.txt, U < 1 certified, and Lemma D re-checked at that Q with the largest N of the row.

## 6. Small Q by exact computation (cert.py sec 5)

Exact subset-product closure mod q for every 2 ≤ q ≤ 2000 (S ← S ∪ pS over the primes p ∤ q, in order; stop when |S| = φ(q)).
The **first 17 primes 2, 3, …, 59 cover every q ≤ 2000**; the first 15 (≤ 47) every q ≤ 1000, 13 (≤ 41) every q ≤ 300,
11 (≤ 31) every q ≤ 100. The 16 least primes > 2000 (2003, …, 2113) cover every q ≤ 2000; 14 primes in (1000, 1091] every q ≤ 1000.
Hardest moduli: 1995, 1955, 1885, 1785 (first primes, m = 17), 1999 (primes > 2000, m = 16). Adversarial moduli rechecked
with the printed sets: 3, 4, 6, 8, 30, 210, 1024, 1155 = 3·5·7·11, 1365, 1001, 1980, 1993, 1997, 1999, 2000. Counting bound:
covering q needs ≥ log₂ φ(q) primes (10.96 for q = 1999).

## 7. Certificate conventions
python-flint Arb at 200 bits; every proof inequality is a certified ball sign (pos), else the run aborts. Floats only choose
r, s, α, a₀, b₀, c' and cell boundaries; the chosen values are then certified. Roundings: all upper bounds are ball upper ends;
input constants as stated (1.01624, BMOR's 0.22737 and 1.567 inside explprog) or rounded up (B). Integers (closed sets, Smith
divisors, closures) are exact. Printed statement constants are ball upper ends rounded up (ceil_str). One command, about 50 s,
35 MB peak (measured), one thread:
`uv run --no-project --with python-flint --with mpmath python3 cert.py all` → logs/cert_output.txt.

## 8. Adversarial checklist
- Moduli with many small prime factors (primorials, 2^k, products of primes of P): only |G[d]|, ω(q) and Σ_{p|q} log p ≤ log q enter
  (§2 (i), §3.4); all q ≤ 2000 are checked exactly (§6). Primes of P dividing q are no-ops (§1), not removed from P.
- q ≤ 4: H = 2, a single real character; included in (3.2) (q = 2 is trivial).
- q prime: G cyclic; |G[2]| = 2; nothing special. 2-power q: (ℤ/2^aℤ)^× = ℤ/2 × ℤ/2^{a−2} handled by |G[d]| ≤ 2d^{ω(q)}.
- Repeated draws: no-ops, priced in δ (§2 (i)). Imprimitive characters: §2 (ψ(p) = ψ*(p) for p ∤ q).
- Real characters (the subgroup obstruction: all p in an index-2 subgroup) give the π = 1/2 classes; they are not the bottleneck.
- No step uses GRH, a zero-density estimate or an unproved zero count; the only unrefereed input is (H) (Lean-checked).

## 9. What limits the constants
- Union over ≈ Q moduli: the collision stage needs (1 + 1/k) log₂ Q primes (k-th moment of C_r); k ≤ 4 because the closed relation
  sets are enumerated (2340 for k = 4). **Open (generic dominance):** if, for every k, every closed R ≠ {0} of rank ρ had
  log₂(1/π_R) > (k − ρ + 1)k/(k + 1), the slope would tend to 1/log 2 + 2/log(1/β) as k → ∞, i.e. (1/log 2 + η) log Q + c_η with
  primes ≤ (c'_η log Q)^E for every η > 0. Least slack log₂(1/π_R) − (k − ρ + 1)k/(k + 1) at δ = 0: 0.0817 bits (k = 2, π = 3/8),
  0.1650 (k = 3, π = 3/16, k − ρ = 2), 0.2000 (k = 4, π = 1/2, the real all-equal tuples, slack 1/(k + 1) for every k); the
  all-equal non-real tuples (π = C(2k,k)/4^k) have slack log₂(4^k/C(2k,k)) − 2k/(k + 1) > 0 for all k ≥ 2. Not proved for k ≥ 5.
- Hole stage: each step kills holes only at rate δ, so it costs b₁ log Q = 2 log Q/log(1/β) (b₁ ≤ 0.293577 at δ = 10⁻³, α = 10⁻⁴; ≤ 0.145771 at δ = 10⁻⁶, α = 10⁻⁷).
- δ ≈ (0.45…0.80)/c': the kernel constant c_∞ = 0.4533 of explprog (72% of it from zeros below height 5/7, BMOR).
- Monte Carlo (mc.py, not part of the proof): random primes from (10⁶, 2·10⁶] cover every q ≤ 100, 300, 1000 with median N = 10, 12, 14
  (log₂ Q = 6.6, 8.2, 10.0); the certified N is about twice the observed value at Q = 10³; the certified slope is ≤ 1.4603 log₂ Q (k = 4, δ = 10⁻³).
