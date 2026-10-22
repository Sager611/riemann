# STLemP: Lemma P lifted to ring class groups; family mean squares of L(½, Sym^m f ⊗ θ_χ); new thresholds for Frobenius conditions in ideal classes

Date 2026-10-22 (arXiv listings fetched 2026-10-21). Brief: `.scratch/orca/stlemp-brief.md`. Checks: [checks/](checks/) (§8; every script under a hard alarm, every process
< 10 MB; the scratch directory `.scratch/stlemp/` with the arXiv listings was deleted). References: [C] = [Combos](../Combos/NOTE.md),
[D2] = [Density2](../Density2/NOTE.md), [R] = [RCApps](../RCApps/NOTE.md), [R2] = [RCApps2](../RCApps2/NOTE.md), [R3] = [RCApps3](../RCApps3/NOTE.md),
[R4] = [RCApps4](../RCApps4/NOTE.md), [R5] = [RCApps5](../RCApps5/NOTE.md) (Lemma P, Corollary P, Theorem 1), [RS] = [RSFam](../RSFam/NOTE.md),
[CS] = [ChebST](../ChebST/NOTE.md) (with [FIX](../ChebST/FIX.md)); HT, PT, BTZ and inputs (J1)–(J6) as in [RS] §0 and [CS] header;
[HB] = Heath-Brown, *A mean value estimate for real character sums*, Acta Arith. 72 (1995), Theorem 1 (the quadratic large sieve).
Notation of [RS] §0 and [CS] §0, §3 throughout (f a non-CM newform of level N_f, λ_m, L_m(s,χ) = L(s, Sym^m f ⊗ θ_χ), 𝔇(X), Ĝ, Ĝ_c,
N_m(σ;t,χ), N_D^{c,m}(σ,T), e(σ) := (1−σ)/(3−2σ)); notation of [R3] §0 and [R5] §2 for orders (𝒟, Δ, f(Δ), G_Δ, T_n, λ_χ(n), Φ,
Q′_Δ, 𝒬′). Paper proofs, not refereed, no Lean. Unconditional, no 2026 half-planes; the preprint PT is used only where marked (the
degree-3 threshold), as in [CS].

## Verdicts

| Item | Result | Status |
|---|---|---|
| task 1 | the binding term of every threshold (§0). None of [RS], [CS] uses [R3] Theorem 1, so the N³X^{−1} term that Lemma P removed does not occur there; the losses a family count can remove are the *per-D* second moments (I₂) and the *per-D* class I sieve below and just above the conductor | §0 |
| tool 1 | Theorem 1′: [R5] Theorem 1 for every sequence whose quadratic twists are controlled (class 𝔅: (B1), (B2)); main-term budget √X′·𝔙E^θ. Contains [R3]'s detector sequences, the Sym^m f-twisted detector sequences and the approximate-functional-equation pieces of L(s, Sym^m f ⊗ θ_χ), m = 1, 2 (Lemma W). Not for arbitrary coefficients ((B2) fails for b ≥ 0 at q = 1) | new |
| tool 2 | Lemma L (ring-class lift): for b vanishing on multiples of a prime p₀, Σ_{D∈𝔇(X)}Σ_{χ∈Ĝ(D)}\|Σ_n b(n)λ_χ(n)\|² ≤ C(p₀²X)^{1/2}log·𝒬′_{p₀²X}(b): class group characters of D lift to the ring class group of conductor p₀, and Lemma P makes all other orders of discriminant ≍ p₀²X nonnegative. Corollary L: the family sieve beyond the conductor, ≪ X^{O(δ)}max{X, N}^{3/2}N^{1−2σ} for class-𝔅 sequences (per-D sum: X·N·N^{1−2σ}) | new; first bound beyond the conductor of [R3] §6, [R5] §10 (lossy by (N/X)^{3/2}) |
| mean squares | Theorem MS: Σ_{D∈𝔇(X)}Σ_{χ∈Ĝ(D)}\|L(½+δ+iv, f⊗θ_χ)\|² ≪ X^{3/2+60δ} (Lindelöf on average for the double family; per-D sum X^{2}) and Σ\|L(½+δ+iv, Sym²f⊗θ_χ)\|² ≪ X^{9/4+60δ} (per-D sum X^{5/2}) | new |
| density, m = 1 | Theorem 2₁: Σ_{D∈𝔇(X)}Σ_{χ∈Ĝ(D)}N_1(σ;t,χ) ≪ X^{2−σ+40δ}, ½+3δ ≤ σ ≤ 1: the per-D density hypothesis for f ⊗ θ_χ, on average over D. Supersedes [RS] Theorem C (per D) and, for σ < 0.87, [RS] Proposition F (X^{11e(σ)}) | new |
| density, m = 2 | Proposition F₂: Σ_DΣ_χN_2(σ;t,χ) ≪ X^{15(5−4σ)/(16(2−σ))+Cδ}, ¾ ≤ σ ≤ 13/14 (per-D mollifiers of length Y^{1/4}, class I by Lemma P below X and by the lift above). [RS] Proposition F: X^{12e(σ)}; at σ = 5/6: 25/18 = 1.389 against 3/2 | new |
| one-sided | Theorem D′, degree-one minorants (λ_f(p) > −α with α > 0, λ_f(p) < α, …): every A > 1 (was 9/4) | **new record**; 1 is the GRH limit of the method |
| sign | Theorem D′, degree-two minorants (λ_f(p) > 0, < 0, ≠ 0, \|λ_f(p)\| ≶ β, λ_f(p) > α for α < ½): every A > A₂* := (23 − 3√5)/8 = 2.036474… (was 3) | **new record** |
| Sato–Tate in classes | [CS] Theorems 1, 2 with 𝒜(0) = 𝒜(1) = 1, 𝒜(2) = A₂*, 𝒜(3) = 3 (PT; 10/3 without PT), 𝒜(M) unchanged for M ≥ 4 (11/3 for 4 ≤ M ≤ 8, (M+3)/3 for M ≥ 8) | M ≤ 3 new; 𝒜(3) = 3 comes from re-optimizing [CS] Proposition F′ (free Y), not from Lemma P |
| discrepancy | [CS] Theorem 2: 10/(3A−3) for A > 11/3, unchanged (it is governed by M ≥ 8) | unchanged |
| D > 0 | the same thresholds for narrow classes (Lemma L⁺, [R5] Lemma P⁺) | §6.5 |
| limits | M ≥ 4: class I at N = X^{11/4}, where PT (X^{3/4}N) meets HT (X^{7/2}); the lift costs N^{3/2} = X^{33/8} there. M ≥ 8: the second-moment slope; the lift gives I₂ ≪ X^{3(m+1)/4}, worse than X^{(m+3)/2} for m ≥ 3 | NOGO for this method (§7) |

## 0. What sets each threshold (task 1)

Write η := 1 − σ and give exponents of X (X^{O(δ)} dropped). For the top twist m and exponent A, zeros with β below
1 − 1/(2A) (or below ½ + 4δ when a family bound reaches down to σ = ½) are counted per D trivially (h log zeros per window, each
weighted by x^{β−1} ≤ |D|^{−A(1−β)}); above, Markov over D needs Σ_{D∈𝔇(X)}N_D^{c,m}(σ, T_X) ≤ X^{1+Aη−ε′} (the *bar*).

| M | old 𝒜(M) | binding point | binding terms | where they come from | does Lemma P touch it? |
|---|---|---|---|---|---|
| 1 | 9/4 ([RS]) | σ = 5/6 | [RS] Theorem C per D, (\|D\|³)^{e(σ)}, against [RS] Proposition F, X^{11e(σ)} | Theorem C: per-D sieve and per-D Lemma B; Proposition F: HT's threshold Q²\|𝒮\| = X^{7/2} times I₂ = Σ_D(per-D Lemma B) = X^{2} | yes: I₂ is a family-sieve problem at length \|D\| = conductor of θ_χ. Lemma P with the lift gives I₂ ≪ X^{3/2} (Theorem MS); then the per-D detector of [CS] §3 with all class I lengths below X^{1−4δ} gives the density hypothesis on average (Theorem 2₁), and Theorem C is not needed |
| 2 | 3 ([RS]) | σ = 5/6 | trivial per-D count, Theorem C, Proposition F (X^{12e}) and the bar coincide at X^{3/2} | Proposition F = (X^{7/2}·I₂)^{2e}, I₂ = X^{5/2} | partly: the AFE length X^{3/2} is beyond the conductor, so Lemma P alone does not reach I₂; the lift does (I₂ ≪ X^{9/4}), and it also handles the class I lengths X ≤ N ≤ X^{15/8} (Proposition F₂) |
| 3 | 10/3 ([CS], route F) | σ = 5/6 | route F, X^{13e} | [CS] route H′ was evaluated at the Y that balances class II against a *diagonal* class I, which drags in N = X^{11/4} | no; with Y as small as class II allows, route H′ (PT) gives 3 (§6.3) |
| 4–8 | 11/3 ([CS]) | σ = 1 − 1/(2A) | class II forces Y > X^{11/4}; at N = X^{11/4} the best of per-D (XN), HT (N + X^{7/2}), PT (NX^{3/4} + X^{5/2}) is X^{7/2}, giving X^{3/4+(11/2)η}, the bar at A = 11/3 | lengths X^{7/4} beyond the conductor | no: Lemma P stops at the conductor and the lift costs N^{3/2} = X^{33/8} > X^{7/2} at N = X^{11/4} (§7) |
| ≥ 8 | (M+3)/3 ([CS]) | σ = 1 − 1/(2A) | class II slope ⅓ per degree from I₂ = X^{(m+3)/2} | AFE length X^{(m+1)/2} ≥ X^{2} | no: the lift gives X^{3(m+1)/4} ≥ X^{(m+3)/2} for m ≥ 3 |
| discrepancy | 10/(3A−3) ([CS] Thm 2) | | 𝒜(2n) = (2n+3)/3, 2n ≥ 8 | the row above | no |

Rankin–Selberg conductor growth enters only through the AFE length |D|^{(m+1)/2} (rows 2, ≥ 8); the HT large-sieve constant enters
only route F; Proposition F's family count |𝒮| ≍ X^{3/2} is used at its true size everywhere; PT enters rows 3–8.

## 1. The family sieve for sequences with controlled quadratic twists

For X′ ≥ 32 and b supported on admissible n (odd, squarefree, every prime factor ≥ 19) put, as in [R5] Corollary P,

  𝒬′_{X′}(b) := Σ_{Δ∈𝒟}Φ(|Δ|/X′)Q′_Δ(b) = ½Σ_{n,n′}b(n)b̄(n′)Σ_{c|(n,n′)}𝒩_c(nn′/c²; X′).

**Definition (class 𝔅).** Let 0 < δ ≤ 1/100, X′ ≥ 32, N = 2^j, ½ − δ ≤ σ ≤ 1, 𝔙 ≥ 1, 0 ≤ θ ≤ ½, C_b ≥ 1. A sequence b is in
𝔅 = 𝔅(N, σ, 𝔙, θ; X′) if supp b ⊂ {admissible n ∈ [N/16, 4N]} and
- (B1) |b(n)| ≤ C_b log(2N)·d(n)^{C_b}n^{−σ};
- (B2) for every admissible κ ≤ 4N, every odd squarefree q ≥ 1 and every character ψ mod 8, with ψ̃ := ψ·(·/q),
  |Σ_{m admissible, (m,κ)=1} b(κm)ψ̃(m)| ≤ C_b X′^{3δ}N^{1/2−σ}κ^{−1/2}d(κ)^{C_b}𝔙^{1/2}·𝓛_b(q,ψ),
  where 𝓛_b(q,ψ) ≥ 0 and Σ_{q≤Q, q odd squarefree}𝓛_b(q,ψ)² ≤ C_b X′^{3δ}Q^{1+θ} for every 1 ≤ Q ≤ X′.

(B1) gives Σ_n d(n)⁴|b(n)|² ≤ C(δ, C_b)X′^{δ}N^{1−2σ} for N ≤ X′. [R3] Lemma 9 (with Jutila's mean value) says exactly that
[R3]'s detector sequences b_{j,t′}, ∂_{t′}b_{j,t′} lie in 𝔅(N, σ, V, 0). Arbitrary sequences do not: for b ≥ 0, (B2) at κ = q = 1,
ψ = 1 asks |Σ_m b(m)| ≤ N^{1/2−σ}𝔙^{1/2}·O(X′^{O(δ)}), false for 𝔙 small; this is the pole term √X′|Σb|² of [R5] Remark (3), the
reason for the reviewer's restriction. Class 𝔅 keeps the restriction: (B2) is a hypothesis to be verified for each sequence (Lemma W).

**Theorem 1′.** Let b ∈ 𝔅(N, σ, 𝔙, θ; X′), X′ ≥ X₀(δ, C_b) and N ≤ X′^{1−4δ}. Then, with E := max{1, 23N/√X′},

  𝒬′_{X′}(b) ≤ C(δ, C_b) X′^{20δ}(X′ + √X′·𝔙·E^θ)·N^{1−2σ}.

*Proof.* This is the proof of [R5] Theorem 1 at scale X′, i.e. [R3] §3.4–§3.8 at w = 1 (Corollary P: 𝒬′_{X′}(b) = ½Σ_cQ_c, Q_c the
w = 1 terms, c odd squarefree, β_c = b(c·) supported on admissible n₂ ∈ [L/64, L], L = 4N/c). That proof uses b only as follows.
(a) Lemma 6 of [R3] (Poisson in T) needs (c, M) = 1, M ≤ X′² and X′ ≥ 64c√M X′^δ for M = n₂n₂′: the first because n = cn₂ is
squarefree, the others because c√M = √(nn′) ≤ 4N ≤ 4X′^{1−4δ} ([R5] Theorem 1, (i)). (b) Lemma 7 (local densities): M = n₂n₂′ is
odd and cube-free. (c) §3.5: S(2^ate√v) = 0 unless 2^ate ≤ E₀ := 23N/(c√X′); if E₀ < 1 the set 𝔈 is empty and the y ≠ 0 part
vanishes. (d) y = 0 (Lemma 5(a)): ≤ 4X′Σ|β_c|². (e) The secondary term §3.6 uses only [HB] and Σ_n d(n)⁴|b(n)|², and (§3.6 (iv)
of [R5]) E₀ + L/κ ≤ 2L/κ, which holds as N ≤ 4cX′/529; it is ≪ X′^{8δ}c^{−1}X′Σ_{c|n}d(n)⁴|b(n)|²; σ enters nowhere. (f) The main
term §3.7 uses b only through the bound for B_{κ,q}(1) = ±Σ_{m:(m,κ)=1}b(cκm)ψ(κm)(m/q) given by Lemma 9, and through the mean value
of 𝓛 over q ≤ E₀. By (B2) (with cκ in place of κ), |B_{κ,q}(1)| ≤ C_bX′^{3δ}N^{1/2−σ}(cκ)^{−1/2}d(cκ)^{C_b}𝔙^{1/2}𝓛_b(q,ψ), which is
Lemma 9's bound with V replaced by 𝔙 and d(·)² by d(·)^{C_b}; and by partial summation over dyadic ranges,
Σ_{q≤E₀}q^{−1}σ_{−1}(q)𝓛_b(q,ψ)² ≤ C(δ)X′^{4δ}E^θ (σ_{−1}(q) ≤ C log log q; E₀ ≤ 23N ≤ X′). Every other step of §3.7 (the weights (2^ate)^{−1},
|ϑ| ≤ 2^{ω(T₄)}T₄T₅, (κκ′)^{−1/2} = (T₂T₃)^{−1/2}(T₄T₅)^{−1}, and the sums over a, t, ψ, T₂, T₃, T₄) is unchanged; the larger divisor
powers cost X′^{O(δ)} through Σ_{T}d(T)^{C}T^{−1} ≪ (log X′)^{C′}. So |𝔐_c| ≪ X′^{13δ}N^{1−2σ}𝔙E^θc^{−1} and the main part of Q_c
is ≪ X′^{13δ}c^{−2}√X′𝔙E^θN^{1−2σ}. Summing (d), (e), (f) over c ≤ 4N gives Σ_c|Q_c| ≤ C(δ, C_b)X′^{15δ}(X′ + √X′𝔙E^θ)N^{1−2σ}. □

For [R3]'s detector sequences (𝔙 = V ≤ X′^{1/2}, θ = 0) this is [R5] Theorem 1.

## 2. The ring-class lift

**Lemma L.** Let X ≥ 32, let p₀ be 1 or an odd prime, X′ := p₀²X, and let b be supported on admissible n with p₀ ∤ n. Then

  Σ_{D∈𝔇(X)}Σ_{χ∈Ĝ(D)}|Σ_n b(n)λ_χ(n)|² ≤ C X′^{1/2}log X′·𝒬′_{X′}(b).

*Proof.* For p₀ = 1 this is [R5] Corollary P. Let p₀ be an odd prime, D ∈ 𝔇(X) and Δ := p₀²D: Δ ≡ D (mod 8) is a discriminant of
conductor f(Δ) = p₀, |Δ| > 4. In the realisation of [R3] §0, G_Δ = I_K(p₀)/P_{K,Z}(p₀) and G_D = I_K/P_K ≅ I_K(p₀)/P_K(p₀), so
π: G_Δ → G_D, [𝔞]_Δ ↦ [𝔞], is a surjective homomorphism, and χ ↦ χ∘π embeds Ĝ(D) into Ĝ_Δ. For n prime to p₀ every O_K-ideal of norm
n is prime to p₀, so λ_{χ∘π}(n) = (χ∘π)(T_n) = Σ_{N𝔞=n}χ(𝔞) = λ_χ(n). Hence, by [R3] Lemma 2(i),
  Σ_{χ∈Ĝ(D)}|Σ_n b(n)λ_χ(n)|² ≤ Σ_{ψ∈Ĝ_Δ}|Σ_n b(n)λ_ψ(n)|² = h(Δ)Σ_{C∈G_Δ}|S_C^{(Δ)}(b)|².
By Lemma P (f = p₀, so e ∈ {1, p₀}), Q′_Δ(b) = Σ_{C∈G_Δ}|S_C^{(Δ)}(b)|² + (w_D/2)Σ_{C∈G_D}|S_C^{(D)}(b_{p₀})|², and
b_{p₀}(m) = b(p₀m)1[(m, p₀) = 1] = 0. So the left side of the lemma is ≤ Σ_{D∈𝔇(X)}h(p₀²D)Q′_{p₀²D}(b). The map D ↦ p₀²D is injective,
Φ(|p₀²D|/X′) = Φ(|D|/X) = 1 on 𝔇(X), and every term Φ(|Δ|/X′)Q′_Δ(b) of 𝒬′_{X′}(b) is ≥ 0 by Lemma P; so the sum is
≤ max_{D∈𝔇(X)}h(p₀²D)·𝒬′_{X′}(b). Finally h(p₀²D) = h(D)(p₀ − (D/p₀)) ≤ (p₀ + 1)C X^{1/2}log X ≤ 2C X′^{1/2}log X′ ([C] (F3)). □

Check (`liftcheck.py`, §8): 43 pairs (D, p₀), random b on admissible n ≤ 400 prime to p₀; h(p₀²D)Σ_C|S_C^{(p₀²D)}|² ≥ h(D)Σ_C|S_C^{(D)}|² in
all cases (smallest ratio 1.98), Σ_C S_C equal on both sides, h(p₀²D) = h(D)(p₀ − (D/p₀)).

**Corollary L (the family sieve beyond the conductor).** Let X ≥ X₀(δ, C_b), N ≥ 1, and choose p₀ := 1 if N ≤ X^{1−4δ}, and
otherwise a prime p₀ ∈ [P, 2P], P := (N^{1+5δ}/X)^{1/2}, p₀ > 576, p₀ ∤ N_f. Put X′ := p₀²X (so N ≤ X′^{1−4δ}, and
X′ ≤ 4max{X, N^{1+5δ}}). If b ∈ 𝔅(N, σ, 𝔙, θ; X′) and p₀ ∤ n on supp b, then
  Σ_{D∈𝔇(X)}Σ_{χ∈Ĝ(D)}|Σ_n b(n)λ_χ(n)|² ≤ C(δ, C_b)X′^{1/2+21δ}(X′ + √X′𝔙E^θ)N^{1−2σ}.
In particular for 𝔙E^θ ≤ X′^{1/2}: ≪ X^{O(δ)}max{X, N}^{3/2}N^{1−2σ}, against X(4N + X^{1/2})N^{1−2σ} for the sum of the per-D sieves.

*Proof.* N ≤ X′^{1−4δ}: for N > X^{1−4δ}, X′ ≥ N^{1+5δ} and (1+5δ)(1−4δ) ≥ 1. Lemma L and Theorem 1′. □

**Lemma L⁺ (D > 0).** For fundamental D > 0, Δ := p₀²D, the narrow ring class group G⁺_Δ surjects onto G⁺_D, and [R4] Lemma 2⁺(i), (ii)
give Σ_{χ∈Ĝ⁺(D)}|Σ_nb(n)λ_χ(n)|² ≤ ŵ(0)^{−1}h⁺(Δ)R⁺_ΔQ_Δ(b), h⁺(Δ)R⁺_Δ ≤ C√Δ log(3Δ)log log(9Δ), and Q_Δ(b) = Q′_Δ(b) for b vanishing on
multiples of p₀ ([R5] Lemma P⁺: the block e = p₀ vanishes). With [R5] Lemma P⁺ (positivity of every non-square order):
Σ_{D∈𝔇⁺(X)}Σ_{χ∈Ĝ⁺(D)}|Σ b λ_χ|² ≤ CX′^{1/2}log²X′·Σ_{Δ>0 non-square}Φ(Δ/X′)Q′_Δ(b), and [R5] Theorem 1⁺ at w = 1 (whose main term is
[R3] §3.7 verbatim, [R4] Theorem 1⁺) gives Theorem 1′ for this sum; square Δ contribute ≪ X′^{δ}N‖b‖² ≤ X′‖b‖² as N ≤ X′^{1−4δ}.

## 3. Quadratic twists of Sym^m f–sequences (m = 1, 2)

For squarefree n, λ_m(n) := Π_{p|n}U_m(x_p) (λ₁ = λ_f, λ₂(p) = λ_f(p)² − 1, |λ_m(n)| ≤ (m+1)^{ω(n)}); P_m := 64(m+1)². For a finite
set 𝒫 of primes > P_m with Σ_{p∈𝒫}p^{−1/2} ≤ 1, call n **(m,𝒫)-admissible** if n is squarefree and every p | n has p > P_m, p ∤ N_f,
p ∉ 𝒫. Such n are admissible in [R3]'s sense. For a Dirichlet character ψ̃, L(w, Sym^m f ⊗ ψ̃) is entire for m = 1, 2 and every ψ̃,
the principal one included (twists of the cuspidal f, Sym²f by characters; f has no CM); L_p(w, Sym^m f ⊗ ψ̃)^{−1} =
Π_{i≤m+1}(1 − α_{m,i}(p)ψ̃(p)p^{−w}) with |α_{m,i}(p)| ≤ 1.

**Lemma H (mean values of quadratic twists).** Let π ∈ {f, Sym²f}, d := 2 resp. 3, 0 < δ ≤ ¼, ψ a character mod 8, u ∈ R, Q ≥ 1. Then
  Σ_{q≤Q, q odd squarefree}|L(½+δ+iu, π ⊗ ψ(·)(·/q))|² ≤ C(δ, f)(Q(|u|+3))^{d/2+δ}.

*Proof.* Let χ := ψ(·)(·/q) and χ* the primitive character inducing it (conductor q·c, c | 8, since (·/q) is primitive for odd
squarefree q); L(s, π⊗χ) = L(s, π⊗χ*)·L₂(s, π⊗χ*)^{−1} with |L₂^{−1}| ≤ 2^d on Re s ≥ 0. π⊗χ* is cuspidal on GL_d with conductor
𝔮 ≍_f q^d and a gamma factor that depends only on χ*(−1); split the q by χ*(−1) = ±1. The approximate functional equation
(Iwaniec–Kowalski Thm 5.3) at s = ½+δ+iu gives L(s, π⊗χ*) = S₁ + ε(s)S₂, S₁ = Σ_nλ_π(n)χ*(n)n^{−s}V_s(n/√𝔮),
S₂ = Σ_nλ_π̄(n)χ̄*(n)n^{s−1}V_{1−s}(n/√𝔮), |ε(s)| ≤ C(𝔮(|u|+3)^d)^{−δ}, (y d/dy)^kV_s(y) ≪_{k,B}(1 + y/(|u|+3)^{d/2})^{−B}. Cut n at
P₀ := (Q(|u|+3))^{d/2}Q^{δ/4} (the tail is ≪ Q^{−10}), split into dyadic blocks n ~ P with the partition φ, and separate the
q-dependence of V by Mellin: φ(v)V_s(Pv/√𝔮) = (2π)^{−1}∫F̂(τ)v^{−iτ}dτ with |F̂(τ)| ≤ C_A(1+|τ|)^{−A} uniformly in P and q. So the
block equals (2π)^{−1}∫F̂(τ)B_τ(q)dτ, B_τ(q) := Σ_nλ_π(n)χ*(n)n^{−s−iτ}φ₂(n/P) (φ₂ = 1 on supp φ), and |block|² ≤ C∫(1+|τ|)^{−A+2}|B_τ(q)|²dτ.
χ*(n) = 0 for even n when c = 8 or 4; in general split n = 2^an′ (the 2-part has ≤ log P values and (2/q) depends on q mod 8, so split
q mod 8 once more). For odd n = n₁n₂² (n₁ squarefree): χ*(n) = ψ*(n₁)(n₁/q)·1[(n₂, q) = 1] (ψ* is quadratic). Cauchy–Schwarz over n₂
with weights n₂^{−1−δ} and dropping the indicator 1[(n₂,q) = 1] ≤ 1 outside the inner sum: |B_τ(q)|² ≤ ζ(1+δ)Σ_{n₂}n₂^{−1−3δ}
|Σ_{n₁ odd squarefree}c_{n₂}(n₁)(n₁/q)|² with |c_{n₂}(n₁)| ≤ d_d(n₁)d_d(n₂²)n₁^{−½−δ}. [HB] for each n₂:
Σ_{q≤Q}|Σ_{n₁}c_{n₂}(n₁)(n₁/q)|² ≤ C(δ)(QP)^{δ/4}(Q + P/n₂²)Σ_{n₁}|c_{n₂}(n₁)|² ≤ C(δ)(QP)^{δ/4}(Q + P)d_d(n₂²)². Summing over n₂ and the
O(log P₀) blocks: Σ_q|S₁|² ≤ C(δ)(QP₀)^{δ/4}(Q + P₀) ≤ C(δ, f)(Q(|u|+3))^{d/2+δ} (d ≥ 2, so P₀ ≥ Q). S₂ is the same with coefficients
of size n^{−½+δ}, i.e. an extra factor ≤ P₀^{2δ}, cancelled by |ε(s)|² ≤ C(𝔮(|u|+3)^d)^{−2δ}. □

**Lemma W (the sequences are in 𝔅).** Let m ∈ {1, 2}, 0 < δ ≤ 1/1000, X ≤ X′ ≤ X^{3}, |t′| ≤ X′^{δ}, 𝒫 as above, and C_b := C(δ, f, m).
- (i) (twisted detector) Let 2 ≤ V ≤ Y, N = 2^j ≤ 4Y log²Y, ½ + 2δ ≤ σ ≤ 1, W(v) := e^{−vN/Y}φ₂(v), and
  b(n) := 1_{(m,𝒫)-adm}(n)a_V(n)λ_m(n)n^{−σ−it′}W(n/N), or ∂_{t′} of it (a_V(n) = Σ_{k|n,k≤V}μ(k); these are the class I sequences
  b_{j,t′}, ∂_{t′}b_{j,t′} of [CS] §3 with (m,𝒫)-admissibility). Then b ∈ 𝔅(N, σ, V, (m−1)/2; X′).
- (ii) (AFE pieces) Let s₀ ≤ X′ be an integer, P ≥ 1, |σ − ½| ≤ δ, w := σ + it′, N := 2^{⌈log₂(P/s₀)⌉}, and
  b(k) := 1_{(m,𝒫)-adm}(k)1[(k, s₀) = 1]λ_m(k)k^{−w}φ₂(s₀k/P). Then b ∈ 𝔅(N, σ, 1, (m−1)/2; X′).

*Proof.* (B1) is clear (|a_V(n)| ≤ d(n), |λ_m(n)| ≤ d(n)^{2}; ∂_{t′} multiplies by −i log n). (B2), case (i): let κ ≤ 4N be admissible,
q odd squarefree, ψ mod 8, ψ̃ := ψ(·/q), s := σ + it′. Terms with κ not (m,𝒫)-admissible vanish. For (m,𝒫)-admissible m coprime to κ,
λ_m(κm) = λ_m(κ)λ_m(m) and a_V(κm) = Σ_{k₁|κ}Σ_{k₂|m, k₁k₂≤V}μ(k₁)μ(k₂). By Mellin inversion,
  Σ_m b(κm)ψ̃(m) = λ_m(κ)κ^{−s}(2πi)^{−1}∫_{(2)}𝒟_κ(s+u)(N/κ)^uŴ(u)du,  Ŵ(u) := ∫_0^∞W(v)v^{u−1}dv,
  𝒟_κ(z) := Σ_{k₁|κ}μ(k₁)Σ_{k₂≤V/k₁, adm, (k₂,κ)=1}μ(k₂)λ_m(k₂)ψ̃(k₂)k₂^{−z}·Σ_{m′ adm, (m′, κk₂)=1}λ_m(m′)ψ̃(m′)m′^{−z}.
The last sum is L(z, Sym^m f ⊗ ψ̃)·Π_{p non-adm or p|κk₂}L_p(z)^{−1}·Π_{p adm, p∤κk₂}h_p(z), h_p := (1 + λ_m(p)ψ̃(p)p^{−z})L_p(z)^{−1} =
1 + Σ_{k≥2}c_kp^{−kz}, |c_k| ≤ 4^{m+1}. On Re z = ½ + δ: |L_p(z)^{−1}| ≤ 2^{m+1}, the product over the non-admissible p is ≤ C(f, m)
e^{(m+1)Σ_{𝒫}p^{−1/2}} ≤ C(f, m), the product over p | κk₂ is ≤ d(κk₂)^{m+1}, and Π_{p adm, p∤r}|h_p| ≤ C(δ, m)2^{ω(r)} (|h_p − 1| ≤
2·4^{m+1}p^{−1−2δ} ≤ ½ for p > P_m). Hence on Re z = ½ + δ,
  |𝒟_κ(z)| ≤ C(δ, f, m)d(κ)^{m+3}Σ_{k₂≤V}d(k₂)^{m+2}k₂^{−½−δ}|L(z, Sym^m f⊗ψ̃)| ≤ C(δ, f, m)d(κ)^{m+3}V^{½}(log 2V)^{2^{m+2}}|L(z, Sym^m f⊗ψ̃)|.
𝒟_κ and Ŵ are entire (supp W ⊂ [¼, 4]), and (v d/dv)^kW ≪_k (5 log²Y)^k on supp φ₂ (N/Y ≤ 4 log²Y), so |Ŵ(u)| ≤ C_A(log 2Y)^{2A}
(1+|Im u|)^{−A} for |Re u| ≤ 2. Move the line to Re(s+u) = ½ + δ (Re u = ½ + δ − σ ∈ [−½, −δ]); no pole is crossed (also for ψ̃
principal: L(z, Sym^m f) is entire). On the new line |κ^{−s}(N/κ)^u| = N^{½+δ−σ}κ^{−½−δ}. Thus (B2) holds with 𝔙 = V and
  𝓛(q,ψ) := ∫_R|L(½+δ+i(t′+r), Sym^m f ⊗ ψ̃)|(1+|r|)^{−8}dr,
the factors N^δ, (log 2Y)^{C}, |λ_m(κ)| ≤ d(κ)² being absorbed in X′^{3δ}d(κ)^{C_b}. By Cauchy–Schwarz in r and Lemma H (d = m+1),
Σ_{q≤Q}𝓛(q,ψ)² ≤ C∫Σ_q|L(½+δ+i(t′+r))|²(1+|r|)^{−8}dr ≤ C(δ, f)Q^{1+(m−1)/2}·Q^{δ}(|t′|+3)² ≤ C_bX′^{3δ}Q^{1+(m−1)/2}
for Q ≤ X′ (|t′| ≤ X′^δ). The ∂_{t′} case: W^♯(v) :=
−i log(Nv)W(v) has the same Mellin bounds. Case (ii): the same computation with V = 1 (no k₁, k₂ sums), W(v) := φ₂(s₀Nv/P)
(smooth, supported in v ∈ [1/16, 4], derivative bounds uniform), the extra coprimality (m′, s₀) = 1 costing d(s₀)^{m+1} ≤ X′^{δ}, and
the line moved from Re u = 2 to Re u = ½ + δ − σ ∈ [0, 2δ] (no pole). □

## 4. Family mean squares

**Theorem MS.** Let m ∈ {1, 2}, 0 < δ ≤ 1/1000, X ≥ X₀(δ, f), |v| ≤ X^{δ}. Then
  Σ_{D∈𝔇(X)}Σ_{χ∈Ĝ(D)}|L_m(½+δ+iv, χ)|² ≤ C(δ, f)X^{3(m+1)/4+60δ}.
(Summing [CS] Lemma F1(d) over D gives X^{(m+3)/2+ε}: X² for m = 1, X^{5/2} for m = 2. For m = 1 the new bound is the size of the
family, i.e. Lindelöf on average for {f ⊗ θ_χ : |D| ~ X, χ ∈ Ĝ(D)}. The same proof on Re s = ½ (σ = ½ for both sums in Lemma W(ii),
|ε_m| = 1) gives Σ_{D,χ}|L_m(½+iv,χ)|² ≪_ε X^{3(m+1)/4+ε}, δ being arbitrary.)

*Proof.* Fix a prime p₀ ∈ [X^{(m−1)/4+6δ}, 2X^{(m−1)/4+6δ}] with p₀ > P_m, p₀ ∤ N_f, and put X′ := p₀²X ∈ [X^{(m+1)/2+12δ},
4X^{(m+1)/2+12δ}], 𝒫 := {p₀}, 𝔅𝔓 := {p ≤ P_m} ∪ {p | N_f} ∪ {p₀}. Let D ∈ 𝔇(X), s′ := ½ + δ + iv, q_m = q_m(D) ≍ |D|^{m+1} the
conductor of L_m(·, χ) (independent of χ, [RS] (R2)). By [RS] (R3), L_m(s′,χ) = S₁(χ) + ε_m(χ,s′)S₂(χ) with
S₁ = Σ_na_χ(n)n^{−s′}V_{s′}(n/√q_m), S₂ = Σ_nā_χ(n)n^{s′−1}V_{1−s′}(n/√q_m), a_χ(n) := Σ_{N𝔞=n}λ_m(𝔞)χ(𝔞), |a_χ(n)| ≤ d_{2m+2}(n),
and |ε_m(χ,s′)| ≤ C(q_m(|v|+3)^{2m+2})^{−δ}.
*Factorisation.* Each n ≥ 1 is uniquely n = s·k with k (m,𝒫)-admissible, (k, s) = 1, and s ∈ 𝒮 := {s : p | s ⟹ p ∈ 𝔅𝔓 or p² | s}
(k := the product of the p ∥ n with p ∉ 𝔅𝔓). a_χ is multiplicative, and a_χ(k) = λ_m(k)λ_χ(k) for (m,𝒫)-admissible k ([CS] (F0)).
*Truncation and separation.* The terms n > P₀ := X^{(m+1)/2+4δ} contribute ≪ X^{−10} to S₁ (decay of V_{s′}; √q_m(|v|+3)^{m+1} ≤
CX^{(m+1)/2+3δ}). For a dyadic P ≤ 2P₀ let S₁^{(P)} be the part with weight φ(n/P). Since (v d/dv)^k[φ(v)V_{s′}(Pv/√q_m)] ≪_k 1
uniformly in P, D, Mellin inversion in v gives φ(n/P)V_{s′}(n/√q_m) = (2π)^{−1}∫F̂_{P,D}(τ)(n/P)^{−iτ}φ₂(n/P)dτ with
|F̂_{P,D}(τ)| ≤ C_A(1+|τ|)^{−A}. Hence |S₁^{(P)}(χ)|² ≤ C_A∫(1+|τ|)^{−A+2}|B_{P,τ}(χ)|²dτ, where
  B_{P,τ}(χ) := Σ_na_χ(n)n^{−s′−iτ}φ₂(n/P) = Σ_{s∈𝒮}a_χ(s)s^{−s′−iτ}Σ_kb_{s,P,τ}(k)λ_χ(k),
  b_{s,P,τ}(k) := 1_{(m,𝒫)-adm}(k)1[(k,s)=1]λ_m(k)k^{−s′−iτ}φ₂(sk/P)
(independent of D and χ). By Cauchy–Schwarz with |a_χ(s)s^{−s′}| ≤ d_{2m+2}(s)s^{−½},
  |B_{P,τ}(χ)|² ≤ 𝔖·Σ_{s∈𝒮, s≤4P}d_{2m+2}(s)s^{−½}|Σ_kb_{s,P,τ}(k)λ_χ(k)|²,  𝔖 := Σ_{s∈𝒮, s≤4P}d_{2m+2}(s)s^{−½} ≤ C(f,m)(log X)^{C(m)},
since the primes of 𝔅𝔓 give Π_{p∈𝔅𝔓}(1 − p^{−½})^{−2m−2} ≤ C(f, m) (the factor of p₀ is ≤ 2) and the squarefull part gives
Π_{p≤4P}(1 + C(m)/p).
*The family sieve.* For |τ| ≤ X^{δ}: b_{s,P,τ} ∈ 𝔅(N, ½+δ, 1, (m−1)/2; X′) with N = 2^{⌈log₂(P/s)⌉} ≤ 8P₀ (Lemma W(ii), t′ = v + τ), it
vanishes on multiples of p₀, and N ≤ 8X^{(m+1)/2+4δ} ≤ X′^{1−4δ}. Lemma L and Theorem 1′ give
  Σ_{D∈𝔇(X)}Σ_{χ∈Ĝ(D)}|Σ_kb_{s,P,τ}(k)λ_χ(k)|² ≤ C(δ, f)X′^{1/2+21δ}(X′ + √X′E^{(m−1)/2})N^{−2δ} ≤ C(δ, f)X′^{3/2+21δ},
because E ≤ 23X′^{1/2}, so √X′E^{(m−1)/2} ≤ 5X′^{½+(m−1)/4} ≤ 5X′ for m ≤ 3. For |τ| > X^δ use, for each D, [C] Lemma 2:
Σ_χ|Σ_kbλ_χ|² ≤ C log²(3X)(4N + X^{1/2})Σ_kτ(k)|b(k)|² ≤ CX^{(m+1)/2+6δ}; with A := 20/δ the weight (1+|τ|)^{−A+2} makes this ≪ X^{−5}.
*Collecting.* Σ_{D,χ}|S₁(χ)|² ≤ (log 4P₀)Σ_PΣ_{D,χ}|S₁^{(P)}|² + O(X^{−5}) ≤ C(δ, f)(log X)^{C}X′^{3/2+21δ}. S₂ is treated identically
(ā_χ = a_{χ̄} since λ_m is real, and χ̄ runs over Ĝ(D)), with k^{−(½−δ)−iv} in place of k^{−s′}: Lemma W(ii) with σ = ½ − δ gives the
extra factor N^{4δ} ≤ (8P₀)^{4δ} ≤ CX^{(2m+2)δ+20δ²}, compensated by |ε_m|² ≤ C(q_m(|v|+3)^{2m+2})^{−2δ} ≤ CX^{−(2m+2)δ}·C. Since
X′^{3/2+21δ} ≤ 8X^{3(m+1)/4+18δ}·(4X^{3/2+12δ})^{21δ} ≤ 16X^{3(m+1)/4+50δ}, the theorem follows. □

## 5. Family densities for m = 1, 2

**The detector.** [CS] Lemmas F1, F2 for m ∈ {1, 2} with "admissible" replaced by "(m,𝒫)-admissible" and δ″ := δ: the primes of 𝒫 move
into the finite product E of Lemma F1(a), which stays entire and bounded by C(f,m)e^{(2m+2)Σ_𝒫p^{−1/2}} on Re s ≥ ½ (the only range
used: the contour of Lemma F2 and class II lie on Re s ≥ ½ + δ); Lemma F1(b), (c), (d)
and Lemma F2 are unchanged (F1(c) needs V ≤ |D|^{1/2}). For each D ∈ 𝔇(X), χ ∈ 𝒳(D) (the χ ∈ Ĝ(D) with a zero ρ_χ of L_m(·,χ) in
[σ,1] × [t−1, t+1]) is in class I or class II. As in [CS] Proposition F′ (proof):
- *Class II* (Cauchy–Schwarz over the pairs (D, χ) jointly, not per D):
  Σ_D|𝒳_II(D)| ≤ C(δ)Y^{½+δ−σ}∫e^{−|v−t|}(Σ_{D,χ}|L_m(½+δ+iv,χ)|²)^{½}(Σ_{D,χ}|M̃(½+δ+iv,χ)|²)^{½}dv,
  with Σ_{D,χ}|M̃|² ≤ C(δ)X^{3/2}(log X)^{C} (Lemma F1(c), V ≤ X^{1/2}/2), Theorem MS for |v| ≤ X^δ and [CS] Lemma F1(d) (per D,
  ≤ C|D|^{(m+1)/2+δ}(|v|+3)^{m+2}) with e^{−|v−t|} ≤ e^{−|v|/2} for |v| > X^δ ≥ 4|t|. Hence
  Σ_D|𝒳_II(D)| ≤ C(δ, f)X^{3(m+3)/8+31δ}Y^{½+δ−σ}  (X^{3/2} for m = 1, X^{15/8} for m = 2).
- *Class I*: Σ_D|𝒳_I(D)| ≤ CJ²Σ_j∫(1+|τ|)^{−A}∫_{I_τ}Σ_{D,χ}(|P_j(χ,t′)|² + |∂_{t′}P_j(χ,t′)|²)dt′dτ over the blocks
  V/2 < N = 2^j < 2Z, P_j(χ,t′) = Σ_nb_{j,t′}(n)λ_χ(n); for |t′| > X^{δ} (only |τ| > X^δ/2) the per-D sieve and the decay of the
  weight give ≪ X^{−5}.

**Theorem 2₁ (m = 1: the density hypothesis of each D, on average).** For 0 < δ ≤ 1/1000, X ≥ X₀(δ, f), ½ + 3δ ≤ σ ≤ 1 and
|t| ≤ X^δ/4,
  Σ_{D∈𝔇(X)}Σ_{χ∈Ĝ(D)}N_1(σ; t, χ) ≤ C(δ, f)X^{2−σ+40δ}.

*Proof.* 𝒫 := ∅, V := X^{1/2}/2, Y := X^{1−5δ}, Z := 2Y log²Y ≤ X^{1−4δ}/4. Class II: X^{3/2+31δ}Y^{½+δ−σ} ≤ X^{2−σ+35δ}, as
(1−5δ)(½+δ−σ) ≤ ½ − σ + 4δ. Class I: every block has N ≤ 2Z ≤ X^{1−4δ}, so Lemma L with p₀ = 1, Theorem 1′ and Lemma W(i)
(𝔙 = V ≤ X^{1/2}, θ = 0, X′ = X) give Σ_{D,χ}|P_j|² ≤ C(δ, f)X^{1/2+21δ}(X + √X·V)N^{1−2σ} ≤ C(δ, f)X^{3/2+21δ}N^{1−2σ}, and the same
for ∂_{t′}P_j. Summing over N > V/2: Σ_D|𝒳_I(D)| ≤ C(δ, f)X^{3/2+22δ}(V/2)^{1−2σ} ≤ C(δ, f)X^{2−σ+23δ}. Finally
Σ_χN_1(σ;t,χ) ≤ C log(XT)|𝒳(D)| ([RS] (R4)). □

For each D the family {f ⊗ θ_χ : χ ∈ Ĝ(D)} has h ≍ |D|^{1/2±ε} members; its density hypothesis is |D|^{1−σ}, i.e. X^{2−σ} over 𝔇(X).
[RS] Theorem C gives each D only |D|^{3e(σ)} (trivial below σ = ¾), [RS] Proposition F gives X^{11e(σ)}; Theorem 2₁ is better than
Proposition F for σ < √14/2 − 1 = 0.8708 (where 11e(σ) = 2 − σ) and than Theorem C summed over D for every σ < 1.

**Proposition F₂ (m = 2).** For 0 < δ ≤ 1/1000, X ≥ X₀(δ, f), ¾ ≤ σ ≤ 13/14 and |t| ≤ X^δ/4,
  Σ_{D∈𝔇(X)}Σ_{χ∈Ĝ(D)}N_2(σ; t, χ) ≤ C(δ, f)X^{F₂(σ)+Cδ},  F₂(σ) := 15(5−4σ)/(16(2−σ)).
(At σ = 5/6: F₂ = 25/18 = 1.3889, against 3/2 for [RS] Proposition F, Theorem C and the trivial count.)

*Proof.* η := 1 − σ ∈ [1/14, ¼], y := 15/(8(1+η)) ∈ [3/2, 15/8], v := y/4 ∈ [3/8, 15/32]; V := X^{v}, Y := X^{y}, Z := 2Y log²Y.
For each block N_j = 2^j ∈ (X^{1−4δ}, 2Z] choose a prime p₀(j) as in Corollary L (P_j = (N_j^{1+5δ}/X)^{1/2} ≥ X^{δ/4} as
(1−4δ)(1+5δ) − 1 ≥ δ/2), and let 𝒫 be the set of these ≤ 2 log Z primes (Σ_𝒫p^{−1/2} ≤ 2 log Z·X^{−δ/8} ≤ 1). Class II: C(δ, f)X^{15/8+31δ}Y^{½+δ−σ} = C(δ, f)X^{c₂+O(δ)},
c₂ := 15/8 − y(½ − η). Class I, with ν := log N/log X:
- blocks with N ≤ X^{1−4δ}: Lemma L (p₀ = 1), Theorem 1′, Lemma W(i) (𝔙 = V, θ = ½): ≤ C(δ, f)X^{1/2+21δ}(X + √X·V·E^{1/2})N^{1−2σ},
  exponent ½ + max{1, ¼ + v + ν/2} + ν(2η−1) + O(δ) (E = 1 when 23N < √X). On ν ∈ [v, 1] its maximum is R(v) := 3/2 − v(1−2η):
  the first branch 3/2 + ν(2η−1) decreases in ν; the second, ¾ + v + ν(2η − ½), is non-increasing (η ≤ ¼) and at its left end
  ν = 3/2 − 2v equals 2v + 3η − 4vη ≤ R(v) (⟺ v ≤ ½).
- blocks with X^{1−4δ} < N < 2Z: Corollary L with p₀ = p₀(j), Theorem 1′, Lemma W(i): X′ ≍ N^{1+5δ}, E ≤ 23N^{(1−5δ)/2}, so the bound
  is ≤ C(δ, f)X^{O(δ)}(N^{3/2} + VN^{5/4})N^{1−2σ}, exponent max{ν/2, v + ν/4} + 2νη + O(δ) ≤ y/2 + 2yη + O(δ) for ν ≤ y + o(1)
  (it increases in ν, and v = y/4).
So Σ_DΣ_χN_2 ≪ X^{O(δ)}(X^{c₂} + X^{R(v)} + X^{y(½+2η)}). With y = 15/(8(1+η)): c₂ = y(½ + 2η) = 15(1+4η)/(16(1+η)) = F₂(σ), and
R(y/4) ≤ F₂(σ) ⟺ η ≥ 1/14 (exact computation, `thresholds.py`). □

Without the class I lift (per-D sieve above X^{1−4δ}) the same optimisation gives 1 + 7η/(2(1+2η)) and, with route F, the
threshold 21/8 below; with the lifted class I but the per-D I₂ (class II X^{2}Y^{½−σ}) route H gives (1+4η)/(1+η) and the threshold
√13 − 1 = 2.6056.

## 6. Primes with Frobenius conditions in every class

### 6.1 Proposition T¹ (m = 1, every A > 1)
**Proposition T¹.** [RS] Proposition T for m = 1 holds for every A > 1 (in place of A > 9/4): for η₀ ∈ (0,1] there are effective
δ₁ = δ₁(f, A, η₀) > 0 and, for δ ∈ (0, δ₁], κ := δ/4 and X ≥ X₀, a set E^{(1)}(X) ⊂ 𝔇(X) with #E^{(1)}(X) ≤ X^{1−δ/8} such that for
D ∈ 𝔇(X) ∖ E^{(1)}(X), every C ∈ G, x ≥ |D|^A and H ∈ [x|D|^{−κ}, x]: |hS_C^{(1)}| ≤ η₀ν_DH.

*Proof.* [RS] §5 for m = 1 with these changes (the statement for A implies it for A′ > A, so let 1 < A ≤ 1.2).
*Parameters.* ε′ := (A−1)/80; δ₁ := min{ε′/100, A/(8c_{(1)}), δ_e} with the BTZ exponent c_{(1)} = 10⁷·2⁴ and δ_e as in [RS] §5.1; κ := δ/4, ε₀ := ε₁ := δ/100,
σ₀ := ½ + 4δ, σ_j := σ₀ + jε′/A (0 ≤ j ≤ J, J minimal with σ_J ≥ 1 − ε′/A), T_X := 2X^{κ+ε₁}, Q₀(X) := (X/2)^δ.
*E₁ (Markov).* For σ_j ≤ 9/10: Theorem 2₁ (with 40δ ≤ ε′/2) on the ≤ 3T_X unit windows gives Σ_DN_D^{c,1}(σ_j,T_X) ≤ C T_X X^{2−σ_j+40δ},
and 2ε′ + (κ+ε₁) + 40δ + (2 − σ_j) − A(1−σ_j) = 1 + 2ε′ + 41δ − (A−1)(1−σ_j) ≤ 1 − 8ε′ + 2ε′ + 41δ ≤ 1 − ε′ ((A−1)(1−σ_j) ≥ (A−1)/10 = 8ε′).
For σ_j > 9/10: [RS] Proposition F (m = 1), exponent ≤ 2ε′ + 11e(σ_j) + 11ε₀ + 12(κ+ε₁) ≤ 2ε′ + 11/12 + 4δ ≤ 1 − ε′ (11e(σ) ≤ 11/12 on
[9/10, 1); ε′ ≤ 1/400). So #E₁(X) ≤ CJX^{1−ε′}.
*E₂* as [RS] §5.2 (BTZ, genus factors |d| > Q₀).
*(b) zeros with β < σ₀*: per window ≤ |D|^{½+o(1)} zeros ((R4)), each with x^{β−1} ≤ |D|^{−A(½−4δ)}: over ≤ 3T₁ windows,
E_b ≤ H|D|^{½−A(½−4δ)+κ+ε₁+o(1)} ≤ H|D|^{−(A−1)/2+6δ}.
*(a), (c), (d), (e), (f)* as [RS] §5.4; (f) and the non-degree-one ideals need only x ≥ |D|^{1+3δ} (h(H/x + x^{−½})log(|D|x) ≤ H|D|^{−2δ}
and H/√x·log²x ≤ ν_DH|D|^{−δ}/h), true as A > 1 + 3δ.
*Conclusion* as [RS] §5.5: ε′ > δ/2, (A−1)/2 − 6δ > δ/2, A/(2c_{(1)}) − κ − ε₁ > δ/2, so E_a + E_b + E_c + E_d + E_f ≤ H|D|^{−δ/2−η₁} ≤ η₀ν_DH/8,
and with (e), |hS_C^{(1)}| ≤ η₀ν_DH/4. □

### 6.2 Proposition T² (m = 2, every A > A₂*)
**Proposition T².** [RS] Proposition T for m = 2 holds for every A > A₂* := (23 − 3√5)/8 = 2.0364745… (in place of A > 3).

*Proof.* Let A₂* < A ≤ A₂* + ¼. *Parameters.* σ₀ := 1 − 1/(2A₂*) = 0.75447… (≥ ¾), s₁ := A(1−σ₀) − ½ = (A − A₂*)/(2A₂*),
ε′ := min{s₁/4, (A − A₂*)/112, 1/32}, δ₂ := min{1/2840, s₁/10, ε′/(10C), A/(8c_{(2)}), δ_e} (C the constant of Proposition F₂'s X^{Cδ},
c_{(2)} = 10⁷·3⁴, BTZ), the rest as in 6.1. *E₁.* For σ_j ∈ [σ₀, 13/14], with η_j := 1 − σ_j ∈ [1/14, 1/(2A₂*)]: Proposition F₂ on ≤ 3T_X
windows; F₂(σ_j) − Aη_j = 1 + η_j(f(η_j) − A) with f(η) := (44η − 1)/(16η(1+η)), and max_{η>0}f = f((1+3√5)/44) = A₂* (exact,
`thresholds.py`); so the exponent is ≤ 2ε′ + (κ+ε₁) + Cδ + 1 − (A − A₂*)/14 ≤ 1 − 8ε′ + 2ε′ + 2Cδ ≤ 1 − ε′. For σ_j ∈ (13/14, 1):
[RS] Proposition F (m = 2), 12e(σ_j) ≤ 12η/(1+2η) ≤ 3/4 for η ≤ 1/14, exponent ≤ 2ε′ + 3/4 + O(δ) ≤ 1 − ε′. *(b)* β < σ₀: the trivial
count as in [RS] §5.4 (m = 2), each window ≤ |D|^{½−A(1−σ₀)+o(1)} = |D|^{−s₁+o(1)}. Everything else as in [RS] §5 (m = 2); (f) and
the non-degree-one ideals need x ≥ |D|^{1+3δ}. □

At A = A₂* the family bound F₂ meets the bar exactly at σ* = 1 − (1+3√5)/44 = 0.82481…; there route F fails (12e(σ*) = 1.557 > 1.357 = 1 + A₂*(1−σ*),
with the lifted I₂ still 11.5e(σ*) = 1.492 > 1.357), and inside route H the class II/class I balance is optimal for the sieves
available (class I at the top block is ≥ y/2 + 2yη for every sieve of §7, class II decreases in y).

### 6.3 m = 3: the free choice of Y in [CS] Proposition F′ (PT)
[CS] Proposition F′ fixed Y := (X^{m+6}T^{2m+2})^{1/(2(3−2σ))}, the balance of class II against a *diagonal* class I; its proof works
for every Y ≥ V and gives, for ½ + 3δ ≤ σ ≤ 1,
  Σ_{D∈𝔇(X)}Σ_{χ∈Ĝ_c(D)}N_m(σ;t,χ) ≪ X^{O(δ)}[X^{(m+6)/4}Y^{½−σ} + max_{V/2<N<2Z}S(N)N^{1−2σ}],  S(N) := min{X(4N + X^{1/2}), N + X^{7/2}, NX^{3/4} + X^{5/2}}
((LS1)–(LS3) of [CS] Lemma F3; (LS3) is PT Corollary 1.2). For m = 3 take Y := X^{3/(1+2η)} (η = 1 − σ ≤ 5/14): class II is
X^{¾+6e(σ)}; in class I, per-D gives ≤ 1 + 2νη ≤ 1 + 3η for ν ≤ 3/2, PT gives 5/2 − ν + 2νη ≤ 1 + 3η on [3/2, 7/4] and 3/4 + 2νη ≤ 3/4 + 2yη
on [7/4, y]. So
  Σ_DΣ_{χ∈Ĝ_c}N_3(σ;t,χ) ≪ X^{max{4−3σ, ¾+6e(σ)}+O(δ)}  (9/14 ≤ σ ≤ 1).
**Proposition T₃.** [CS] Proposition T_m for m = 3 holds for every A > 3 (was 10/3). *Proof.* [CS] §2, route H′, with σ₀ := 1 − 1/6 = 5/6,
σ₂ := 1 − 1/52 (as in [CS] §2.1 for m = 3), [CS] Proposition F″ for σ_j ≥ σ₂ exactly as in [CS] §2.2, and the bound above for
σ_j ∈ [5/6, σ₂) (margins: ε′ := min{(A−3)/208, (min{A, 27/8} − 3)/24, 1/32}): 1 + (3 − A)η ≤ 1 − (A − 3)/52, and with A′ := min{A, 27/8} ∈ (3, 27/8], g(η) := ¾ + 6η/(1+2η) − A′η is increasing on (0, 1/6] (g′ = 6/(1+2η)² − A′ ≥ 27/8 − A′ ≥ 0)
with g(1/6) = 3/2 − A′/6 < 1, and ¾ + 6η/(1+2η) − Aη ≤ g(η). The per-D count
covers σ < 5/6 with s₁ = (A−3)/6. The rest is [CS] §2 verbatim (genus twists by BTZ and Thorner's (J5)). □ Without PT the bound is
X^{1+5(1−σ)} in the middle range and 10/3 (route F) stays. For m = 4, A < 11/3 forces Y > X^{11/4} at σ = 1 − 1/(2A) (§7), so 11/3 stays.

### 6.4 Theorems
**Theorem D′** ([RS] Theorem D, new thresholds). Let P(x) = c₀ + c₁x + c₂(x² − 1) with c₀ > 0 and P ≤ 1_S on [−2, 2]. Put A₀ := 1 if
c₂ = 0, A₀ := A₂* = (23 − 3√5)/8 if c₂ ≠ 0. For every A > A₀ there are effective c, κ > 0 and X₀ (depending on f, A, P) such that all
D ∈ 𝔇 with |D| ≤ X outside a set of size ≤ X^{1−c} satisfy, for every C ∈ G, x ≥ |D|^A and H ∈ [x|D|^{−κ}, x],
  Σ_{𝔭∈C, N𝔭 = p ∤ ND, x < p ≤ x+H, λ_f(p) ∈ S} log p ≥ c₀ν_DH/(32h(D)).
*Proof.* [RS] proof of Theorem D with: m = 0, [R5] Theorem 3 (every A > 1; its proof is [R] §2 with changed steps, ending with
hS_C^{(0)} ≥ ν_DH/8 outside its exceptional set, [R] §2.5); m = 1, Proposition T¹; m = 2, Proposition T²; δ := the minimum of the
admissible δ's (each is an upper bound); and x ≥ |D|^{1+3δ} instead of x ≥ |D|² for the non-degree-one ideals (their number is
≤ C(H/√x + ω(ND) + log x)log x ≤ c₀ν_DH/(64h·max|P_𝔞|) once √x ≥ |D|^{½+2δ}). □

**Corollary D′1.** Let E/Q have no CM, f = f_E. (i) For every A > A₂*, all but O(X^{1−c}) D ∈ 𝔇 with |D| ≤ X satisfy: every primitive
positive definite form of discriminant D represents primes p, p′ ∤ ND, p, p′ ≤ 2|D|^A, with a_E(p) > 0 > a_E(p′); the same for
a_E(p) ≠ 0, |a_E(p)| < β√p (β > 1), |a_E(p)| > β√p (β < 1), a_E(p) > α√p (α < ½); and for all but O(N^{1−c}) squarefree n ≤ N,
n ≡ 1, 2 (mod 4), there is a prime p = a² + nb² ≤ 2(4n)^A with a_E(p) > 0 (and one with a_E(p) < 0). (ii) For every A > 1, the same
with a_E(p) > −α√p and a_E(p) < α√p (α > 0 fixed). *Proof.* [RS] Corollary D1 (minorants of degree 2, resp. 1, table of [RS] §5). □

**Sato–Tate in every class.** [CS] Theorems 1, 2 (both signs of D), Corollaries 1.1, 2.1 hold with 𝒜 replaced by
  𝒜′(0) = 𝒜′(1) = 1, 𝒜′(2) = A₂*, 𝒜′(3) = 3 (PT; 10/3 without), 𝒜′(M) = 𝒜(M) for M ≥ 4,
which is nondecreasing. *Proof.* [CS] proof of Theorem 1 with Proposition T¹, T², T₃ for m ≤ 3, [CS] §2 for m ≥ 4, [R5] Theorem 3 for
m = 0, and x ≥ |D|^{1+3δ} for the non-degree-one ideals; Theorem 2 needs R ≤ |D|^{−½}S^{(0)}, i.e. x ≥ |D|^{2+3δ}, true as
𝒜′(2n) ≥ A₂* > 2 for n ≥ 1. □ For an interval J, A(J) = 𝒜′(M*(J)); new only when M*(J) ≤ 3. The discrepancy bound of [CS] Theorem 2
is nontrivial only for n ≥ 5 (5/(n+1) < 1), where 𝒜′(2n) = 𝒜(2n): the exponent 10/(3A−3) (A > 11/3) is unchanged.

### 6.5 D > 0
Theorem D′, Corollary D′1 (indefinite forms, a² − nb², n ≡ 2, 3 (mod 4)) and the Sato–Tate statements hold for fundamental D > 0 and
narrow classes with the same thresholds. Changes: Lemma L⁺ for Lemma L; [R5] Theorem 1⁺ for Theorem 1′ (its main term is [R3] §3.7,
so (B2) applies verbatim); Lemma W does not see D; in Theorem MS the gamma factor of L(s, Sym^m f × π_χ) depends on the archimedean type
a ∈ {0, 1} of π_χ ([RS] §6), so the two types are treated separately, and the per-D tail uses [D2] Lemma 5; the detector is [CS] §3,
"D > 0"; Propositions T¹, T², T₃ follow [RS] §6 and [CS] §6 (m = 0: [R5] Theorem 3⁺; genus factors as in [R2] Theorem 4).
Sectors ([CS] Theorem 3, D < 0, 𝒜_sec(1) = 8/3) were not redone: the lift of angular characters to orders was not checked.

## 7. Thresholds, limits, and what GRH gives

| M | 1 | 2 | 3 | 4–8 | ≥ 8 |
|---|---|---|---|---|---|
| [RS]/[CS] | 9/4 | 3 | 10/3 | 11/3 | (M+3)/3 |
| this note | **1** | **(23−3√5)/8 = 2.03647** | **3** (PT) | 11/3 | (M+3)/3 |
| what carries it | Theorem 2₁ (density hypothesis on average) + Prop F near σ = 1 | Prop F₂ (lifted I₂, lifted class I, V = Y^{1/4}) on [σ₀, 13/14] + Prop F | free Y in [CS] F′ | [CS] | [CS] |
| binding term | the explicit-formula limit (zeros near ½) | F₂ = bar at σ* = 0.8248 (class II with I₂ = X^{9/4} against class I at N = Y via the lift) | per-D at N = X^{3/2} and PT at N = Y | class I at N = X^{11/4} | class II slope |

Intermediate versions (exact, `thresholds.py`): m = 2 with the lifted I₂ but per-D class I: 21/8 (routes F with b = 23/2 and H meet); with the
lifted class I but per-D I₂: √13 − 1 = 2.6056; with both: A₂*. Numerical optimisation over all parameters (`thr2.py`, `thr4.py`, grids 0.0025–0.0125)
gives 1 (m = 1), 2.039 (m = 2; grid of V), 3 (m = 3), 11/3 (m = 4, 5, 8) and reproduces the old values (§8).
- *M ≥ 4 is the limit of this method.* At σ = 1 − 1/(2A), A < 11/3, m = 4: class II below the bar needs y > 2A/(A−1) > 11/4; at N = X^{11/4}
  every available sieve has exponent ≥ ¾ + (11/2)η (per-D 1 + 5.5η, HT ¾ + 5.5η, PT ¾ + 5.5η, lift 11/8 + 5.5η), and
  ¾ + 11/(4A) < 3/2 ⟺ A > 11/3. The length X^{11/4} is X^{7/4} beyond the conductor; Lemma P reaches the conductor only, and the lift pays
  (N/X)^{3/2}.
- *The slope ⅓ per degree (M ≥ 8)* is class II with I₂ = X^{(m+3)/2}; the lift gives X^{3(m+1)/4}, which is better only for m ≤ 2.
- *GRH* (for all L(s, Sym^m f ⊗ θ_χ), m ≤ M, and L(s, χ)): every D, every M, every A > 1 (explicit formula with absolute values:
  h·x^{1/2}log² ≤ |D|^{−κ}x/h for x ≥ |D|^{1+3κ}). So M ≤ 1 reaches the GRH value; below 1 nothing is possible by this method ([R5] §7).
- *A family sieve at the family size for every length* (numerically true for the twisted sequences, §8) would give by route H
  (class I diagonal above X^{3/2}, class II with the per-D I₂) 𝒜(M) = (M+3)/3 for M ≥ 3, and with the lifted I₂ 𝒜(2) = 3/2 (balance
  15/8 − y(½−η) = 2yη, i.e. X^{(15/2)e(σ)}; checked in `thresholds.py`). Lindelöf on average for every Sym^m family as well would give 1 for every M ([CS] §8).

## 8. Numerical checks ([checks/](checks/))
- `thresholds.py` → `thresholds.out`: 51 exact checks (sympy), 0 failures: the m = 1 exponents; y* = 15/(8(1+η)), F₂, the lift exponent
  max{ν/2, ν/4 + v} + 2νη, R(y/4) ≤ F₂ ⟺ η ≥ 1/14, f(η) = (44η−1)/(16η(1+η)), η* = (1+3√5)/44, A₂* = (23−3√5)/8, η* ∈ [1/14, 1/(2A₂*)],
  route F on (0, 1/14] and its failure at η*; 21/8 and √13 − 1 for the intermediate versions; the hypothetical 3/2; m = 3: y = 3/(1+2η), A = 3 ± 1/100; m = 4: the
  forced length 11/4 and the sieve exponents there; I₂ lift vs per-D for m ≤ 60; the discrepancy choice of n.
- `thr2.py`, `thr4.py` → `thr.out`: numerical optimisation of routes F and H over V, Y, X′ (sieves per-D, HT, PT, and the lift with the
  crude main-term budget). It reproduces the old values 9/4 (m = 1, Theorem C + route F), 3 (m = 2), 11/3 (m = 4, 5, 8), 13/3 (m = 10),
  and without PT 10/3, 4, 5 (m = 3, 5, 8), and finds 3 for m = 3 with free Y, 1 for m = 1 and 2.039 for m = 2 with the lift (grid
  optimum; the exact value is (23 − 3√5)/8 = 2.0365). `drv.sh` runs `famqE.c`; `apcheck.py` checks a_p(37a).
- `liftcheck.py` → `liftcheck.out`: Lemma L on 43 pairs (D, p₀) (§2).
- `famqE.c` (from [R5] `checks/famq.c`) → `famqE.out`: the family sums by brute force over all orders, X = 20000, for f = f_E, E = 37a
  (a_p by point counting, checked against a_5…a_43): mode 3: b = λ_E(n)n^{−½}φ₂(n/N) (unmollified, as in the AFE); mode 4: λ_{Sym²}(n)n^{−½}φ₂;
  mode 5: a_V(n)λ_E(n)n^{−5/6}φ₂ (V = √X/2); mode 6: a_V(n)λ_{Sym²}(n)n^{−5/6}φ₂ (V = X^{1/4}); mode 2: n^{−5/6}φ₂ (pole present).
  Entries 𝒬′/(X‖b‖²) / F/(X^{3/2}‖b‖²) (F = the fundamental family sum Σ_DΦh(D)Σ_C|S_C|², the left side of Lemma L with p₀ = 1):

| N | mode 3 | mode 4 | mode 5 | mode 6 | mode 2 (pole) |
|---|---|---|---|---|---|
| 2500 | 0.574 / 0.159 | 0.582 / 0.160 | 0.575 / 0.156 | 0.578 / 0.156 | 18.68 / 3.92 |
| 10000 | 0.573 / 0.159 | 0.586 / 0.160 | 0.583 / 0.157 | 0.586 / 0.156 | 72.69 / 15.03 |
| 20000 | 0.571 / 0.159 | 0.591 / 0.161 | 0.582 / 0.158 | 0.587 / 0.157 | 144.7 / 29.92 |
| 40000 | 0.573 / 0.160 | 0.604 / 0.163 | 0.588 / 0.159 | 0.596 / 0.159 | 289.5 / 59.91 |
| 80000 | 0.569 / 0.160 | 0.603 / 0.163 | 0.585 / 0.158 | 0.603 / 0.161 | 578.8 / 119.8 |

  The twisted sequences (no mollifier needed) stay at the family size from N = X/8 to N = 4X (support up to 16X), as (B2) predicts
  for N ≤ X and better than Corollary L ((N/X)^{3/2}) beyond; the pole-carrying sequence grows like N. Runtime ≤ 2 min per entry, ≤ 3 MB.

## 9. Novelty (2026-10-21)
- arXiv monthly listings math.NT 2026-09 (812 entries) and 2026-10 (625 entries): all titles scanned for Sato–Tate, Hecke eigenvalue,
  ideal class, class group, quadratic form, Frobenius, sign, Chebotarev, ring class, dihedral, theta, represent, least prime, Linnik,
  large sieve, zero density, imaginary quadratic (262 hits), abstracts read for 2609.07080 (Yang, Chebotarev–Sato–Tate bounds for
  multiplicative functions over polynomial values), 2609.10446 (Jääsaari, signs of GL_n Hecke eigenvalues in short intervals),
  2609.17806 (Atkin–Lehner equidistribution), 2610.03976, 2610.11210 (non-vanishing of twisted GL₂ L-values), 2610.11974 (prime
  quadratic twists), 2610.04254 (ideal classes in Ennola cubic fields), 2610.18491 (Frei–Sofos–Widmer). None has a Hecke-eigenvalue or
  Sato–Tate condition on primes in ideal classes of varying quadratic fields, a family mean value of L(½, f ⊗ θ_χ) over D and χ, or a
  family sieve for class group characters beyond the conductor.
- arXiv search (approach "arxiv", 6 queries: second moment Rankin–Selberg × class group characters × average over discriminants;
  large sieve × class group characters × theta series; zero density × class group twists; primes represented × Hecke eigenvalue ×
  ideal class; Lindelöf on average × L(½, f×θ_χ); ring class lift × large sieve): no relevant hit (two queries refused as irrelevant).
- The thresholds 9/4, 3, 10/3 being improved are this project's ([RS], [CS], audited NEW there); nothing in the literature was found
  with any power of |D| for these problems.

## 10. Dead ends, NOGO candidates, pointers
- *Lemma P against the N³X^{−1} term in [RS]/[CS]*: there is none; neither uses [R3] Theorem 1. NOGO as a diagnosis.
- *𝒜(M) < 11/3 for 4 ≤ M ≤ 8* with per-D, HT, PT and the lift: impossible (§7). Needs a family sieve at N ≍ X^{11/4} better than X^{7/2}.
- *Lifting the HT mollifier* (route F with a family mollifier of length < X^{7/2}): the mollifier coefficients μ_{Sym^m f}(n)λ_χ(n)
  violate (B2) (their twists are 1/L(·, Sym^m f ⊗ ψ̃), no continuation to Re < 1), so Theorem 1′ does not apply. NOGO.
- *Lift for m ≥ 3 mean squares*: X^{3(m+1)/4} ≥ X^{(m+3)/2}; no gain. NOGO.
- *Lead*: the twisted family sums stay at the family size far beyond the conductor (§8). A proof (the shifted-convolution structure of
  the y ≠ 0 lattice points when Poisson in T fails) would give 𝒜(2) = 3/2 and 𝒜(M) = (M+3)/3 for M ≥ 3 (§7).
- *Pointer*: Theorem MS (m = 1) and Theorem 2₁ give, by the [R5] argument, A > 1 for the twisted problem; [CS] Theorem 3 (sectors,
  M = 1, 8/3) should drop to 1 once Lemma P^∠ is lifted to orders (not checked).
