# RCApps3: a family large sieve for ring class characters, and the almost-all exponent (3+√69)/10

Date 2026-10-21. Brief: `.scratch/orca/rcapps3-brief.md`. Files in `.scratch/rcapps3/`: `dualcheck.py`, `dualcheck2.py`,
`dualcheck3.py` → `dualcheck*.out` (the dual identities of §1, exact on all 28 test pairs), `localcheck.py` (Lemma 7),
`poissoncheck.py` (Lemma 6 against brute force), `famls.c` + `mkb.py` → `numerics1.md`, `numerics2.txt` (the family sum 𝒬(b)
computed by brute force over all orders, §7), `thresh.py` (thresholds), `notes0.md`–`notes2.md` (working notes), `ax*.xml`,
`ax3_queries.txt`, `axparse*.py`, `gxl3_*.json`, `rutgers.txt` (novelty). References: [C] = [Combos NOTE](../Combos/NOTE.md), [D2] = [Density2 NOTE](../Density2/NOTE.md),
[R] = [RCApps NOTE](../RCApps/NOTE.md), [R2] = [RCApps2 NOTE](../RCApps2/NOTE.md). Paper proofs, not refereed, no Lean, no
2026 half-planes.

## Verdicts

| Item | Result | Status |
|---|---|---|
| tool | Theorem 1: for the class I sequences b of a new zero detector (§2), Σ_{D∈𝔇(X)}Σ_{χ∈Ĝ(D)}\|Σ_n b(n)λ_χ(n)\|² ≪ X^{O(δ)}(X^{3/2} + N³X^{−1})N^{1−2σ} for N ≤ X^{1−4δ} | new; the sum of the per-D sieves is X·N·N^{1−2σ}; X^{3/2} is the size of the family |
| density | Theorem 2: Σ_{D∈𝔇(X)}N_D(σ;t) ≪ X^{(13−8σ)/(7−2σ)+ε} for 1/2+ε ≤ σ ≤ 9/10, \|t\| ≤ X^{ε} (zeros with β ≤ 9/10 of all χ ∈ Ĝ(D)∖{1} in a unit window) | new; at σ = 5/6: 19/16 (was 5/4) |
| least primes, D < 0 | Theorem 3: [R] Theorem 2 (primes in every ring class in (x, x + x\|D\|^{−κ}], x ≥ \|D\|^A, for all but O(X^{1−c}) D) and its corollary for primes a² + nb² (without the arc condition of [R2] Cor 4′) hold for every A > A₀ := (3+√69)/10 = 1.130662… | new record; was 3/2 ([R2]); GRH gives 1 |
| A₀ → 1 | needs Theorem 1 up to N ≤ X^{1−δ} without the N³/X term; numerically the family sum stays at the family size even for N = 1.5X | open (§6) |
| sectors, D > 0 | unchanged at 3/2 ([R2] Thm 3′, Thm 4): the dual count needs angles (sectors) or unit-reduced lattice points (D > 0) | not done |

## 0. Notation

- 𝒟: negative discriminants Δ (Δ ≡ 0, 1 mod 4), not necessarily fundamental. Δ = f²D_K with D_K fundamental, f = f(Δ) the
  conductor, O_Δ = Z + fO_K, G_Δ := Cl(O_Δ) ≅ I_K(f)/P_{K,Z}(f) (classes of O_K-ideals prime to f), h(Δ) = |G_Δ|. For |Δ| > 4,
  O_Δ^× = {±1}, and h(Δ) ≤ C|Δ|^{1/2}log(3|Δ|) ([C] (F3)). For odd p: p | f ⟺ p² | Δ.
- 𝔇(X) as in [R]: negative fundamental D with X/2 < |D| ≤ X; Ĝ(D) = Ĝ_Δ for Δ = D.
- For M ≥ 1: T_M := Σ_{N𝔞=M, (𝔞,f)=1}[𝔞] ∈ Z[G_Δ] (O_K-ideals). For χ ∈ Ĝ_Δ and n squarefree, λ_χ(n) := χ(T_n) (for D
  fundamental, the n-th Hecke eigenvalue of the theta series θ_χ). [C](T_M) := #{𝔞 ∈ C : N𝔞 = M, (𝔞,f)=1}; [1](·) is the
  coefficient of the identity.
- (a/c) is the Jacobi symbol (c odd); for Δ ∈ 𝒟 and odd squarefree c, (Δ/c) = Π_{p|c}(Δ/p).
- Φ ∈ C^∞, supp Φ ⊂ [1/4, 2], 0 ≤ Φ ≤ 1, Φ = 1 on [1/2, 1]. K(v) := v^{−1}∫_0^1 Φ(s/v)(1−s)^{−1/2}ds for v > 0; K extends
  to a C^∞ function on [0, ∞) with K(0) = ∫Φ > 0 and K(v) = 0 for v ≥ 8 (substitute s = vu: K(v) = ∫Φ(u)(1−vu)^{−1/2}du for
  v < 1/2, and the formula ∫_0^1Φ(s/v)(1−s)^{−1/2}ds is smooth in v > 0).
- P := Π_{p≤17}p. "n admissible" means n odd squarefree with (n, P) = 1 (all prime factors ≥ 19).
- δ ∈ (0, 1/100] is a small parameter; C, C(δ) are absolute resp. δ-dependent effective constants; ε-losses are X^{O(δ)}.

## 1. Dual identities

**Lemma 1 (Hecke relations in the group ring).** Let Δ ∈ 𝒟, |Δ| > 4, and n, n' admissible, g := (n,n'), n = gm, n' = gm'.
Then T_nT_{n'} = Σ_{c|g}(Δ/c)T_{mm'(g/c)²} in Z[G_Δ].

*Proof.* T_{ab} = T_aT_b for coprime a, b (unique factorisation of ideals prime to f). For a prime p ∤ f: if p splits,
p = 𝔭𝔭̄, then T_p² = [𝔭]² + [𝔭̄]² + 2[(p)] and T_{p²} = [𝔭²] + [𝔭̄²] + [(p)], [(p)] = 1, so T_p² = T_{p²} + 1; if p is inert,
T_p = 0 and T_{p²} = [(p)] = 1, so T_p² = T_{p²} − 1; if p ramifies, T_p = [𝔭], T_p² = [(p)] = 1 = T_{p²}. In all cases
T_p² = T_{p²} + (Δ/p). If p | f then T_p = T_{p²} = 0 and (Δ/p) = 0. Hence T_g² = Π_{p|g}(T_{p²} + (Δ/p)) = Σ_{c|g}(Δ/c)T_{(g/c)²},
and T_nT_{n'} = T_{mm'}T_g². □

**Lemma 2 (the dual identity).** Let b be supported on admissible n, X ≥ 1, and S_C := Σ_n b(n)[C](T_n) = Σ_{𝔞∈C}b(N𝔞). Then
(i) Σ_{χ∈Ĝ_Δ}|Σ_n b(n)λ_χ(n)|² = h(Δ)Σ_{C∈G_Δ}|S_C|²;
(ii) 𝒬(b) := Σ_{Δ∈𝒟, |Δ|>4}Φ(|Δ|/X)Σ_{C∈G_Δ}|S_C|² = ½ Σ_{n,n'}b(n)b̄(n') Σ_{c,w: cw|(n,n'), (c,w)=1} μ(w) 𝒩_c(nn'/(cw)²; X/w²),
where for M ≥ 1, X' > 0 and odd squarefree c,

  𝒩_c(M; X') := Σ_{(T,y,Δ)} Φ(|Δ|/X')(Δ/c), over T, y ∈ Z and Δ ∈ 𝒟 with T² − Δy² = 4M and T ≡ yΔ (mod 2).

*Proof.* (i) is Parseval on G_Δ. (ii) Σ_C|S_C|² = [1](SS*) with S := Σ_n b(n)T_n and S* := Σ_n b̄(n)T_n, because the involution
𝔞 ↦ 𝔞̄ preserves norms and coprimality to f and [𝔞̄] = [𝔞]^{−1}. By Lemma 1, [1](T_nT_{n'}) = Σ_{c|g}(Δ/c)[1](T_{M_c}),
M_c := nn'/c². Now [1](T_M) counts the 𝔞 ∈ P_{K,Z}(f) prime to f with N𝔞 = M; such 𝔞 = (γ) with γ ∈ O_Δ unique up to ±1, and
(γ) is prime to f iff (M, f) = 1. So [1](T_M) = ½#{γ ∈ O_Δ : Nγ = M}·1[(M,f)=1], and γ = (T + y√Δ)/2 with T ≡ yΔ (2) runs over
O_Δ. Write M_c = mm'r², r = g/c. If p ∥ M_c and p | f, no γ exists (p² | Δ forces T² ≡ 4M_c mod p², impossible since
v_p(4M_c) = 1), so on the support 1[(M_c,f)=1] = 1[(r,f)=1] = Σ_{w|r}μ(w)1[w|f] = Σ_{w|r}μ(w)1[w²|Δ]. For w | r the map
(T,y,Δ) ↦ (T/w, y, Δ/w²) is a bijection from the solutions of T² − Δy² = 4M_c with w² | Δ onto the solutions of
T₂² − Δ₂y² = 4M_c/w² (w² | Δy² and w² | 4M_c give w | T; w is odd, so parities and Δ₂ ≡ Δ (mod 4) are preserved), and
(Δ/c) = (Δ₂/c) since (c, w) = 1, Φ(|Δ|/X) = Φ(|Δ₂|/(X/w²)). The terms |Δ| ≤ 4 are absent because Φ(|Δ|/X) = 0 for X ≥ 32. □

*Check.* `dualcheck.py` (5 coprime pairs; its 4 pairs with (n,n') > 1 show that the naive lattice count fails there),
`dualcheck2.py` (Lemma 1 summed over Δ, 7 pairs) and `dualcheck3.py` (Lemma 2 with a sharp range of Δ, 16 pairs incl. n = n' and
three common primes) compute the left side from reduced forms of all discriminants in [20, 1500] and agree exactly with the right side.

For coprime n, n' Lemma 2 says: the number of pairs of ideals of norms n, n' in one ring class, summed over all orders, equals
the number of lattice points (T, y) with y² | T² − 4nn' and (T² − 4nn')/y² a discriminant in the range — no class numbers and
no squarefree condition. Restricting to fundamental discriminants would bring back a squarefree sieve on (4nn' − T²)/y² with
moduli up to X^{1/2}; positivity (§3.1) lets us include all orders instead.

## 2. A zero detector with pole-free, D-independent class I coefficients

Fix the prime ℓ := 2. For V ≥ 1 put a_V(n) := Σ_{k|n, k≤V}μ(k). For D ∈ 𝔇(X) and χ ∈ Ĝ(D)∖{1} put

  𝒢(s,χ) := Σ_{n admissible} a_V(n)λ_χ(n)n^{−s},   𝒢^F(s,χ) := F(s)𝒢(s,χ),   F(s) := 1 − ℓ^{2−2s}.

F is a Dirichlet polynomial in the rational ideal (ℓ), on which every χ is 1; F(1) = 0 kills the pole of the trivial character.

**Lemma 3.** (a) For Re s > 1, 𝒢 = L(s,χ)·E·H·M̃ with E := Π_{p|2P}L_p(s,χ)^{−1}, H := Π_{p∤2P}(1 + λ_χ(p)p^{−s})L_p(s,χ)^{−1},
M̃ := Σ_{k≤V adm.}μ(k)k^{−s}Π_{p|k}ψ_p(λ_χ(p)), ψ_p(x) := x/(1 + xp^{−s}); L_p(s,χ)^{−1} = 1 − λ_χ(p)p^{−s} + (D/p)p^{−2s}.
E is entire; H and M̃ are holomorphic in Re s > 1/2; on Re s ≥ 1/2 + δ, |E| ≤ C and |H| ≤ C(δ). So 𝒢 continues to Re s > 1/2
and vanishes at every zero of L(s,χ) there.
(b) The coefficient of 𝒢 at n = 1 is 1, and at every 1 < n ≤ V it is 0.
(c) On Re s = 1/2 + δ, for V ≤ |D|^{1/2}: Σ_{χ∈Ĝ(D)}|M̃(s,χ)|² ≤ C(δ)|D|^{1/2}log^C(3|D|).
(d) |F(ρ)| ≥ ℓ^{2(1−β)} − 1 ≥ ℓ^{1/5} − 1 if Re ρ = β ≤ 9/10; |F(s)| ≤ 1 + ℓ on Re s ≥ 1/2.

*Proof.* (a) λ_χ is multiplicative on coprime squarefree arguments, so Σ_{n adm.}a_V(n)λ(n)n^{−s} = Σ_{k≤V adm.}μ(k)λ(k)k^{−s}
Π_{p∤2Pk}(1 + λ(p)p^{−s}). Further (1 + λp^{−s})L_p^{−1} = 1 + ((D/p) − λ²)p^{−2s} + λ(D/p)p^{−3s} with |λ_χ(p)| ≤ 2, so H converges
absolutely and uniformly on Re s ≥ 1/2 + δ. For p ≥ 19 and Re s ≥ 1/2, |λp^{−s}| ≤ 2/√19 < 1/2, so (1 + λp^{−s})^{−1} is
holomorphic and bounded by 2 there. (b) For admissible n ≤ V all divisors are ≤ V, so a_V(n) = Σ_{k|n}μ(k). (d) is clear.
(c) Fix s with Re s = 1/2 + δ. For p ≥ 19 and |x| ≤ 2 write ψ_p(x) = x(1 − θ_p(x)), θ_p(x) := xp^{−s}/(1 + xp^{−s}), so
|θ_p| ≤ 4p^{−1/2}. Expanding Π_{p|k}(1 − θ_p) = Σ_{t|k}(−1)^{ω(t)}Π_{p|t}θ_p and writing k = tk' gives

  M̃(s,χ) = Σ_{t≤V adm.}Θ_t(χ)·M_t(χ),  Θ_t(χ) := μ(t)(−1)^{ω(t)}λ_χ(t)t^{−s}Π_{p|t}θ_p(λ_χ(p)),  M_t(χ) := Σ_{k'≤V/t, adm., (k',t)=1}μ(k')λ_χ(k')k'^{−s},

with |Θ_t(χ)| ≤ 8^{ω(t)}t^{−1}. Cauchy–Schwarz (for each χ) gives |M̃|² ≤ (Σ_{t≤V}8^{ω(t)}t^{−1})Σ_t 8^{ω(t)}t^{−1}|M_t(χ)|², and
the first factor is ≤ C log^8(3V). M_t(χ) = Σ_𝔞 c_𝔞χ(𝔞) over ideals of squarefree norm k' ≤ V/t (all primitive), c_𝔞 = μ(N𝔞)N𝔞^{−s};
by [C] Lemma 2, Σ_χ|M_t|² ≤ C log²(3|D|)(V/t + |D|^{1/2})Σ_{k'≤V}d(k')k'^{−1} ≤ C log⁴(3|D|)(V/t + |D|^{1/2}). Summing over t:
Σ_χ|M̃|² ≤ C log^{C}(3|D|)(V + |D|^{1/2}). □

**Lemma 4 (detector).** Let 1/2 + 2δ ≤ σ ≤ 9/10, Y ≥ V ≥ 2, Z := 2Y log²Y, χ ∈ Ĝ(D)∖{1}, |D| ≥ D₀(δ), and let ρ = β + iγ be a
zero of L(s,χ) with σ ≤ β ≤ 9/10. Put c₀ := (2^{1/5} − 1)/8. Then either
(I) |Σ_j F_j(χ,ρ)| ≥ c₀, j over 2^{j+1} > V, 2^{j−1} < Z, where
  F_j(χ,ρ) := Σ_{n adm.}a_V(n)λ_χ(n)[φ(n/2^j)n^{−ρ}e^{−n/Y} − ℓ²φ(ℓ²n/2^j)(ℓ²n)^{−ρ}e^{−ℓ²n/Y}]
  (φ the dyadic partition of unity of [C] Theorem 1); or
(II) ∫_R|𝒢^F(½+δ+i(γ+u),χ)|e^{−|u|}du ≥ c₁Y^{β−1/2−δ}, c₁ > 0 absolute.

*Proof.* As [C] Lemma 5 for 𝒢^F(s,χ) = Σ_n c_nn^{−s}, c_n := 1_{adm}(n)a_V(n)λ(n) − ℓ²1_{adm}(n/ℓ²)a_V(n/ℓ²)λ(n/ℓ²):
(2πi)^{−1}∫_{(2)}𝒢^F(ρ+w)Γ(w)Y^wdw = Σ_n c_nn^{−ρ}e^{−n/Y}. Moving to Re w = 1/2 + δ − β < 0 crosses only w = 0, with residue
𝒢^F(ρ) = 0 (Lemma 3(a)). By Lemma 3(b) the left side is e^{−1/Y} − ℓ^{2−2ρ}e^{−ℓ²/Y} + Σ_j F_j(χ,ρ) + (tail n > Z, ≤ c₀ for
|D| ≥ D₀), and |e^{−1/Y} − ℓ^{2−2ρ}e^{−ℓ²/Y}| ≥ |F(ρ)| − 2ℓ²/Y ≥ 8c₀ − o(1) by Lemma 3(d). If (I) fails, the right side, which is
≤ (2π)^{−1}∫|𝒢^F(½+δ+i(γ+u))||Γ(½+δ−β+iu)|Y^{1/2+δ−β}du, is ≥ 4c₀; Stirling gives (II). □

**Removing β; the class I sequences.** As in [C] Theorem 1, with φ₂ ∈ C_c^∞((1/4,4)), 0 ≤ φ₂ ≤ 1, φ₂ = 1 on supp φ (inserted
for free) and Φ_β(v) := φ(v)v^{σ−β} (smooth, compact support, uniformly for 0 ≤ β − σ ≤ 1/2, so |Φ̂_β(τ)| ≤ C_A(1+|τ|)^{−A}):
|F_j(χ,ρ)| ≤ (2π)^{−1}∫|Φ̂_β(τ)||P_j(χ,γ+τ)|dτ, where

  P_j(χ,t') := Σ_n b_{j,t'}(n)λ_χ(n),  b_{j,t'}(n) := 1_{adm}(n)a_V(n)n^{−σ−it'}W_{j,t'}(n/2^j),
  W_{j,t'}(v) := e^{−v2^j/Y}φ₂(v) − ℓ^{2−2σ−2it'}e^{−ℓ²v2^j/Y}φ₂(ℓ²v).

(Apply the β-removal to u ↦ φ(u/2^j)u^{−ρ} at u = n and at u = ℓ²n: φ(u/2^j)u^{−ρ} = 2^{j(σ−β)}(2π)^{−1}∫Φ̂_β(τ)2^{ijτ}u^{−σ−i(γ+τ)}dτ,
with |2^{j(σ−β)}| ≤ 1; with t' := γ + τ, ℓ²(ℓ²n)^{−σ−it'} = ℓ^{2−2σ−2it'}n^{−σ−it'}, so the two pieces combine into W_{j,t'}; φ₂ = 1 on the
supports.) The key property is

  (Z)  ∫_0^∞ W_{j,t'}(v)v^{−σ−it'}dv = 0   (substitute u = ℓ²v in the second term).

For the Sobolev step we also need ∂_{t'}b_{j,t'}(n) = 1_{adm}(n)a_V(n)n^{−σ−it'}W^♯(n/2^j), W^♯(v) := −i log(2^jv)W_{j,t'}(v)
+ ∂_{t'}W_{j,t'}(v); (Z) holds for W^♯ because ∫W^♯(v)v^{−σ−it'}dv = −i(log 2^j)·0 + d/dt'[∫W_{j,t'}(v)v^{−σ−it'}dv] = 0. Both
weights are smooth, supported in [1/(4ℓ²), 4], and (v d/dv)^kW ≪_k log(2+2^j) uniformly for 2^j ≤ 4Z.

## 3. The family large sieve

**Theorem 1.** Let 0 < δ ≤ 1/100, X ≥ X₀(δ), 1/2 + 2δ ≤ σ ≤ 9/10, |t'| ≤ X^δ, 2 ≤ V ≤ X^{1/2}, and N := 2^j with
V/4 ≤ N ≤ X^{1−4δ}. Let b be b_{j,t'} or ∂_{t'}b_{j,t'} (§2). Then

  Σ_{D∈𝔇(X)}Σ_{χ∈Ĝ(D)}|Σ_n b(n)λ_χ(n)|² ≤ C(δ)X^{40δ}(X^{3/2} + N³X^{−1})N^{1−2σ}.

For comparison, the per-D sieve ([C] Lemma 2) gives X(N + X^{1/2})N^{1−2σ}, and the left side is ≥ the trivial-character
terms Σ_D|Σ_n b(n)r_D(n)|², which for a weight without (Z) are ≍ X N^{2−2σ}(log V)^{−2} (numerics: §7). Write ‖b‖² := N^{1−2σ}.
All b considered satisfy |b(n)| ≤ C log(2N)·d(n)n^{−σ}, supp b ⊂ [N/(4ℓ²), 4N], so Σ_n d(n)^A|b(n)|² ≤ C(A)X^{δ}‖b‖².

### 3.1 Positivity and the dual identity
Every D ∈ 𝔇(X) is in 𝒟 with f = 1 and Φ(|D|/X) = 1. Adding the nonnegative terms of all other Δ ∈ 𝒟 and using Lemma 2(i)
and h(Δ) ≤ CX^{1/2}log X:

  LHS ≤ Σ_{Δ∈𝒟}Φ(|Δ|/X)Σ_{χ∈Ĝ_Δ}|Σ_n b(n)λ_χ(n)|² ≤ CX^{1/2}(log X)·𝒬(b).

By Lemma 2(ii), 𝒬(b) = ½Σ_{c,w}μ(w)Q_{c,w}, over odd squarefree coprime c, w, where with β := b(cw·) (supported on admissible
n₂ ≤ 4N/(cw) prime to cw) and X' := X/w²,

  Q_{c,w} := Σ_{n₂,n₂'}β(n₂)β̄(n₂')𝒩_c(n₂n₂'; X').

So it suffices to show Σ_{c,w}|Q_{c,w}| ≤ C(δ)X^{39δ}(X + N³X^{−3/2})‖b‖².

### 3.2 Lattice point counts

**Lemma 5.** For M ≥ 1, X' ≥ 1, odd squarefree c: (a) the terms y = 0 of 𝒩_c(M;X') vanish unless M = □, and then contribute at
most 4X'; (b) the terms y ≠ 0 contribute at most C(ε)M^ε·min{√M + 1, √X' + √(M/X')}.

*Proof.* (a) y = 0 forces T² = 4M, and Σ_{Δ}Φ(|Δ|/X') ≤ 2X'. (b) Φ(|Δ|/X') ≠ 0 forces X'/4 ≤ |Δ| ≤ 2X', so y² ≤ 16M/X'. Each
(T, y) determines Δ, and for fixed T, y² | 4M − T² ≠ 0 leaves ≤ d(4M − T²) values of y: the count is ≤ Σ_{|T|≤2√M}d(4M − T²)
≪ M^ε(√M + 1). For fixed y, T lies in ≤ 2^{ω(y)+2}(y², 4M)^{1/2} classes mod y² and in two intervals of total length
≤ y²X'/√M·2 (as 4M − T² ∈ y²[X'/4, 2X']), so the count is ≤ Σ_{y≤4√(M/X')}2^{ω(y)+2}(y²,4M)^{1/2}(2X'/√M + 1)
≪ M^ε(√X' + √(M/X')). □

### 3.3 Large w
Let G₂ := X^{1−2δ}/(64N) ≥ X^{2δ}/64. For w > G₂ use Lemma 5 with √M ≤ 4N/(cw), √(M/X') ≤ 4N/(c√X) and (Σ|β|)² ≤
(4N/(cw))Σ|β|²:

  |Q_{c,w}| ≪ X^{δ}[(N/(cw))·min{N/(cw) + 1, √X/w + N/(c√X)} + X/w²]·Σ_{n₂}|b(cwn₂)|².

Sum over c and w > G₂ by grouping n = cwn₂: Σ_{c, w>G₂}|Q_{c,w}| ≪ X^{δ}Σ_n|b(n)|²Σ_{cw|n, w>G₂}[…]. For w ≤ X^{1/2} the bracket
is ≤ N√X/w² + N²/(w√X) + X/w² ≤ N√X/G₂² + N²/(G₂√X) + X/G₂² ≪ X^{4δ}(N³X^{−3/2} + N²X^{−1}); for w > X^{1/2} it is
≤ 2N²/X + N/√X + 1. Since N ≤ X, N²X^{−1} ≤ X. Hence Σ_{c, w>G₂}|Q_{c,w}| ≪ X^{6δ}(N³X^{−3/2} + X)‖b‖².

### 3.4 Small w: Poisson in T

**Lemma 6.** Let c be odd squarefree, (c, M) = 1, X' ≥ 64c√M X^{δ}, M ≤ X². Then the y ≠ 0 part of 𝒩_c(M; X') equals

  (μ(c)/c)(X'/√M)Σ_{y≥1, (y,c)=1}A_M(y)K(y²X'/(4M)) + O_δ(X^{−10}),
  A_M(y) := ¼#{T₀ mod 4y² : y² | T₀² − 4M, (T₀² − 4M)/y² ≡ 0, 1 (mod 4)}.

*Proof.* For fixed y ≥ 1 (and −y), Δ = (T² − 4M)/y² ∈ 𝒟 iff T lies in the classes counted by A_M(y) mod 4y²; then
T ≡ yΔ (2) is automatic (T² ≡ Δy² (4)), and (Δ/c) depends on T mod 4y²c. Poisson summation over each class mod q := 4y²c of
f(t) := Φ((4M − t²)/(y²X')) gives main terms q^{−1}∫f and errors bounded by Σ_{k≠0}|f̂(k/q)|. Put u(t) := (4M − t²)/(y²X');
on supp f, u ∈ [1/4, 2], so y²X' ≤ 16M, |u'| ≤ 4√M/(y²X'), |u''| = 2/(y²X') ≤ (4√M/(y²X'))², u''' = 0; hence
‖f^{(i)}‖₁ ≪_i √M(4√M/(y²X'))^i and |f̂(k/q)| ≪_i √M(q·4√M/(|k|y²X'))^i = √M(16c√M/(|k|X'))^i ≤ √M(X^{−δ}/(4|k|))^i.
With i := ⌈12/δ⌉ the errors are ≪ X^{−11}, also after summing over the ≤ q classes and the ≤ 4√M ≤ 4X values of y. The main
term: ∫f = 2√M∫_{−1}^{1}Φ((1−τ²)·4M/(y²X'))dτ = (y²X'/(2√M))K(y²X'/(4M)) (substitute t = 2√Mτ, s = 1 − τ²). The local sum:
if (y, c) = 1, CRT and Jacobsthal (Σ_{T mod p}((T² − 4M)/p) = −1 for p ∤ M) give Σ_{T₀ mod q}[…]((Δ/c)) = 4A_M(y)μ(c); if a
prime p divides (y, c), with k := v_p(y) the p-part is Σ over the roots T₁ mod p^{2k} of T₁² ≡ 4M and over j mod p of
(((T₁ + jp^{2k})² − 4M)/p^{2k}·u/p) (u a fixed unit), = Σ_j((A₁ + 2T₁j)u/p) = 0 because p ∤ 2T₁. Collecting the factor 2 for ±y
gives the formula. □

For c, w with w ≤ G₂ and admissible n₂, n₂' ≤ 4N/(cw) prime to c: M = n₂n₂' ≤ 16N²/(cw)², so 64c√M X^δ ≤ 256NX^{δ}/w ≤
X/w² = X' (w ≤ G₂ = X^{1−2δ}/(64N)). Lemma 6 applies.

**Lemma 7 (local densities).** For odd M, A_M(y) = A^{(2)}_M(2^{v₂(y)})·Π_{p^k∥y, p odd}A_M(p^k), where A^{(2)}_M(1) = A^{(2)}_M(2) = 1,
A^{(2)}_M(4) = 2·1[M ≡ 1 (4)], A^{(2)}_M(2^a) = 4·1[M ≡ 1 (8)] (a ≥ 3); and for odd p, k ≥ 1: A_M(p^k) = 1 + (M/p) if p ∤ M,
0 if p ∥ M, p if p² ∥ M and k = 1, p(1 + ((M/p²)/p)) if p² ∥ M and k ≥ 2.

*Proof.* CRT; for odd p ∤ M, T² ≡ 4M (p^{2k}) has 1 + (M/p) roots (Hensel); if p ∥ M there are none; if p² ∥ M then
T = pT₁, T₁² ≡ 4M/p² mod p^{2k−2}, with p roots mod p² for k = 1 and p(1 + ((M/p²)/p)) roots mod p^{2k} for k ≥ 2. The 2-adic
values are a finite computation (T = 2T', T'² ≡ M mod 4^{a−1}, and the condition mod 4 on the quotient). `localcheck.py` checks
all cases against the definition. □

Consequently, on odd y, A_M = 1 ∗ h_M with h_M multiplicative, h_M(p^k) = 0 for k ≥ 3, and for odd p:
h_M(p) = (M/p) if p ∤ M, −1 if p ∥ M, p − 1 if p² ∥ M; h_M(p²) = p((M/p²)/p) if p² ∥ M, 0 otherwise.
For M = n₂n₂' with admissible n₂, n₂' one checks case by case

  h_M(p) = (n₂/p)(n₂'/p) − 1[p|n₂] − 1[p|n₂'] + (p+1)1[p|n₂]1[p|n₂'],  h_M(p²) = p·1[p|n₂]1[p|n₂']((n₂/p)/p)((n₂'/p)/p).

Expanding the product over p | e: for odd cube-free e prime to c,

  (S)  h_{n₂n₂'}(e) = Σ_{e = T₁T₂T₃T₄T₅²} ϑ(T₂,T₃,T₄,T₅)·F_{T₂T₄T₅, T₁, T₅}(n₂)·F_{T₃T₄T₅, T₁, T₅}(n₂'),

over factorisations into pairwise coprime squarefree T_i, with ϑ := (−1)^{ω(T₂T₃)}Π_{p|T₄}(p+1)·T₅ (so |ϑ| ≤ 2^{ω(T₄)}T₄T₅) and
F_{κ,q,s}(n) := 1[κ|n]·(n/q)·((n/s)/s) (for s | κ). Each factorisation is determined by the triple (κ, κ', q) := (T₂T₄T₅,
T₃T₄T₅, T₁T₅): T₅ = (q, κ), T₄ = (κ, κ')/T₅, T₂ = κ/(T₄T₅), T₃ = κ'/(T₄T₅), T₁ = q/T₅. For n = κm with (m, κ) = 1 one has
F_{κ,q,s}(κm) = (κ/T₁)((κ/s)/s)·(m/q): a fixed sign times the quadratic character (·/q) of the cofactor. The 2-adic factor is
A^{(2)}_{n₂n₂'}(2^a) = Σ_{ψ mod 8}â_a(ψ)ψ(n₂)ψ(n₂') with |â_a(ψ)| ≤ 4 (it depends on n₂n₂' mod 8 only).

### 3.5 Small w: Mellin in y
Let g(x) := K(x²) (even, smooth, supported in |x| ≤ √8, g(0) > 0) and S(u) := Σ_{k odd ≥ 1}g(ku). Its Mellin transform is
Š(z) = (1 − 2^{−z})ζ(z)ǧ(z), ǧ(z) := ∫_0^∞g(x)x^{z−1}dx; ǧ is meromorphic with simple poles at 0, −2, −4, …, so Š is
meromorphic in Re z > −2 with a single simple pole, at z = 1, residue c₁ := ½ǧ(1) > 0 (the pole of ǧ at 0 is cancelled by
1 − 2^{−z}), and Š decays faster than any power of |Im z| on −δ ≤ Re z ≤ 2. Also S(u) = 0 for u ≥ √8.

With v := X'/(4M), Σ_{y≥1,(y,c)=1}A_M(y)K(y²v) = Σ_{a≥0}A^{(2)}_M(2^a)Σ_{e odd,(e,c)=1}h_M(e)Σ_{t|c}μ(t)S(2^ate√v) (write
y = 2^aek, k odd prime to c, and remove (k, c) = 1 by Möbius). On our ranges √v ≥ c√X/(8N), so S(2^ate√v) = 0 unless
(a,t,e) ∈ 𝔈 := {2^ate ≤ E₀}, E₀ := 23N/(c√X); restrict to this finite set. Insert S(u) = (2πi)^{−1}∫_{(2)}Š(z)u^{−z}dz, swap
with the finite sums, and move to Re z = −δ, crossing only z = 1. With v^{−1/2} = 2√M/√X' and Lemma 6:

  Σ_{n₂,n₂'}β(n₂)β̄(n₂')𝒩_c^{y≠0}(n₂n₂';X') = (μ(c)/c)[2c₁√X'·𝔐_{c,w} + X'(2πi)^{−1}∫_{(−δ)}Š(z)(4/X')^{z/2}𝔖_{c,w}(z)dz] + O(X^{−8}),
  𝔐_{c,w} := Σ_{(a,t,e)∈𝔈}μ(t)(2^ate)^{−1}Σ_{n₂,n₂'}β(n₂)β̄(n₂')A^{(2)}_{n₂n₂'}(2^a)h_{n₂n₂'}(e),
  𝔖_{c,w}(z) := Σ_{(a,t,e)∈𝔈}μ(t)(2^ate)^{−z}Σ_{n₂,n₂'}β(n₂)β̄(n₂')(n₂n₂')^{(z−1)/2}A^{(2)}_{n₂n₂'}(2^a)h_{n₂n₂'}(e).

By (S), each inner double sum equals Σ_ψ â_a(ψ)Σ ϑ·B_{κ,q}(z)·B'_{κ',q}(z), summed over the triples (κ, κ', q) with
T₁⋯T₅² = e, where (L := 4N/(cw) bounds the support of β)

  B_{κ,q}(z) := Σ_{n: κ|n} β(n)n^{(z−1)/2}ψ(n)F_{κ,T₁,T₅}(n),  B'_{κ',q}(z) := Σ_{n': κ'|n'} β̄(n')n'^{(z−1)/2}ψ(n')F_{κ',T₁,T₅}(n').

Put A_κ := Σ_{κ|n}d(n)|β(n)|².

### 3.6 The secondary term (Re z = −δ)
On Re z = −δ, |(2^ate)^{−z}| ≤ E₀^{δ}, |(4/X')^{z/2}| ≤ X'^{δ/2}, and |n^{(z−1)/2}| ≤ 64L^{−(1+δ)/2} on supp β ⊂ [L/64, L]. Fix (a, t, ψ, T₂, T₃, T₄, T₅) and sum over T₁ (it only enters through q): by Cauchy–Schwarz and Heath-Brown's
quadratic large sieve [HB95, Thm 1],

  Σ_{q≤Q odd sqfree}|Σ_{m≤L', m odd sqfree}γ_m(m/q)|² ≤ C(ε)(QL')^ε(Q + L')Σ_m|γ_m|²,

applied with Q = E₀, L' = L/κ (the coefficients γ_m of m = n/κ are fixed: β(κm)(κm)^{(z−1)/2}ψ(κm)·sign), we get
Σ_{T₁}|B_{κ,q}B'_{κ',q}| ≪ X^{δ}(E₀ + L/κ)^{1/2}(E₀ + L/κ')^{1/2}L^{−1−δ}A_κ^{1/2}A_{κ'}^{1/2}. Since κ ≤ e ≤ E₀ and
E₀²/L ≤ 23²Nw/(cX) ≤ 1 (w ≤ G₂), E₀ + L/κ ≤ 2L/κ, so this is ≪ X^{δ}L^{−δ}(κκ')^{−1/2}A_κ^{1/2}A_{κ'}^{1/2}. With
|ϑ| ≤ 2^{ω(T₄)}T₄T₅ and T₄T₅(κκ')^{−1/2} = (T₂T₃)^{−1/2}, and dropping the constraint (a,t,e) ∈ 𝔈 in the remaining sums of
nonnegative terms (only 2^a, t ≤ E₀ matter, ≤ X^δ choices):

  |𝔖_{c,w}(z)| ≪ X^{3δ}Σ_{G=T₄T₅}2^{ω(G)}(Σ_{T₂}T₂^{−1/2}A_{T₂G}^{1/2})² ≪ X^{4δ}Σ_G2^{ω(G)}Σ_{T₂}A_{T₂G} ≪ X^{5δ}Σ_n d(n)^4|β(n)|²,

using (Σ_{T₂≤E₀}T₂^{−1/2}A_{T₂G}^{1/2})² ≤ (Σ_{T₂}T₂^{−1})(Σ_{T₂}A_{T₂G}). Hence the secondary part of Q_{c,w} is
≪ X^{7δ}c^{−1}(X/w²)Σ_{cw|n}d(n)^4|b(n)|², and summing over c, w gives ≪ X^{8δ}X‖b‖².

### 3.7 The main term (z = 1)
At z = 1, B_{κ,q}(1) = ±Σ_{m: (m,κ)=1}b(cwκm)ψ(κm)(m/q) (the factor ((n/s)/s) is ((κ/s)/s)(m/s), absorbed: the twist on m
is (·/q), q = T₁T₅). By Lemma 9 below (with κ replaced by cwκ),

  |B_{κ,q}(1)| ≤ C(δ)X^{3δ}N^{1/2−σ}(cwκ)^{−1/2}d(cwκ)²V^{1/2}𝓛(q,ψ),  𝓛(q,ψ) := ∫|L(½+δ+i(t'+u), ψ(·)(·/q))|(1+|u|)^{−4}du.

The weight (2^ate)^{−1} = (2^at)^{−1}(T₁T₂T₃T₄T₅²)^{−1}, |ϑ| ≤ 2^{ω(T₄)}T₄T₅, and (κκ')^{−1/2} = (T₂T₃)^{−1/2}(T₄T₅)^{−1}. So

  |𝔐_{c,w}| ≪ X^{8δ}N^{1−2σ}(cw)^{−1}V·Σ_{a,t,ψ}(2^at)^{−1}Σ_{T₂,T₃,T₄}(T₂T₃)^{−3/2}2^{ω(T₄)}T₄^{−1}·Σ_{q≤E₀}q^{−1}σ_{−1}(q)𝓛(q,ψ)²
          (using Σ_{T₁T₅=q}T₁^{−1}T₅^{−2} = q^{−1}Σ_{T₅|q}T₅^{−1})
          ≪ X^{9δ}N^{1−2σ}V(cw)^{−1}·max_ψ∫Σ_{q≤E₀}q^{−1}|L(½+δ+i(t'+u),ψ(·)(·/q))|²(1+|u|)^{−4}du ≪ X^{11δ}N^{1−2σ}V/(cw),

using 𝓛² ≤ C∫|L|²(1+|u|)^{−4}du and the classical mean value Σ_{q≤Q}^♭|L(½+δ+iu,ψ(·)(·/q))|² ≤ C(δ)(Q(1+|u|))^{1+δ}
(Jutila 1981; or [HB95] with the approximate functional equation). Therefore
Σ_{c,w}c^{−1}2c₁√X'|𝔐_{c,w}| ≪ X^{11δ}√X V N^{1−2σ}Σ_{c,w}c^{−2}w^{−2} ≪ X^{11δ}X‖b‖² (V ≤ X^{1/2}).

**Lemma 9 (twisted sums; the pole is killed).** Let κ ≤ 4N be admissible, q odd squarefree (not necessarily prime to κ), ψ a character mod 8,
ψ̃ := ψ(·)(·/q). Then

  |Σ_{m adm., (m,κ)=1} b(κm)ψ̃(m)| ≤ C(δ)X^{3δ}N^{1/2−σ}κ^{−1/2}d(κ)²V^{1/2}∫|L(½+δ+i(t'+u),ψ̃)|(1+|u|)^{−4}du.

*Proof.* Take b = b_{j,t'} (the ∂_{t'} case is identical with W^♯). Write b(κm) = a_V(κm)(κm)^{−s}W(κm/N), s := σ + it'.
Mellin: the sum is (2πi)^{−1}∫_{(2)}𝒟(s+u)κ^{−s}(N/κ)^uŴ(u)du with Ŵ(u) := ∫_0^∞W(v)v^{u−1}dv and
𝒟(s') := Σ_{m adm., (m,κ)=1}a_V(κm)ψ̃(m)m^{−s'}. Since a_V(κm) = Σ_{k₁|κ}Σ_{k₂|m, k₁k₂≤V}μ(k₁)μ(k₂),

  𝒟(s') = [L(s',ψ̃)/L(2s',ψ̃²)]·Π_{p|2Pκ}(1 + ψ̃(p)p^{−s'})^{−1}·Σ_{k₁|κ}μ(k₁)Σ_{k₂≤V/k₁ adm.,(k₂,κ)=1}μ(k₂)ψ̃(k₂)k₂^{−s'}Π_{p|k₂}(1 + ψ̃(p)p^{−s'})^{−1}.

On Re s' = 1/2 + δ: |L(2s',ψ̃²)^{−1}| ≤ ζ(1+2δ), the product over p | 2Pκ is ≤ C·d(κ), and the double sum is ≤ d(κ)Σ_{k₂≤V}2^{ω(k₂)}
k₂^{−1/2} ≤ Cd(κ)V^{1/2}log V (each (1 + ψ̃(p)p^{−s'})^{−1} has modulus ≤ (1 − 19^{−1/2})^{−1} < 2 for p ≥ 19). Move the u-line to
Re(s + u) = 1/2 + δ. A pole is crossed only if ψ̃ is principal (q = 1, ψ = ψ₀), at s + u = 1, with residue
Res_{s'=1}𝒟·κ^{−s}(N/κ)^{1−s}Ŵ(1 − s), and Ŵ(1 − s) = ∫W(v)v^{−σ−it'}dv = 0 by (Z). On the new line
|κ^{−s}(N/κ)^u| = N^{1/2+δ−σ}κ^{−1/2−δ} and |Ŵ(u)| ≤ C(1 + |Im u|)^{−4}log⁴(2N) (four integrations by parts in log v). □

### 3.8 Conclusion
The y = 0 terms (n₂ = n₂', Lemma 5(a)) give Σ_{c,w}4(X/w²)Σ|b(cwn₂)|² ≪ X^{2δ}X‖b‖². With §3.3, §3.6, §3.7:
Σ_{c,w}|Q_{c,w}| ≪ X^{12δ}(X + N³X^{−3/2})‖b‖², and §3.1 gives Theorem 1 (with room in the exponent of X^δ). □

*Remarks.* (1) Without (Z) the residue in Lemma 9 is ≍ N^{1−σ}κ^{−1}(log V)^{−1} for q = 1, and 𝔐_{c,w} contains
√X·|Σ_n b(n)|²-type terms of size √X N^{2−2σ}(log V)^{−2}: the trivial character's contribution; this is the X·N of the per-D sieve.
(2) The only place where the size of N enters beyond X^{1−4δ} is G₂ (Poisson needs X' ≫ c√M): pairs n, n' sharing a factor w > G₂ are
bounded trivially in §3.3. This costs N³X^{−3/2}, and is the only obstruction to the full family size (§6).

## 4. Density over the family

For D ∈ 𝔇(X), t ∈ R, 1/2 < σ ≤ 9/10 let N^♮_D(σ;t) := Σ_{χ∈Ĝ(D)∖{1}}#{ρ : L(ρ,χ) = 0, σ ≤ β ≤ 9/10, |γ − t| ≤ 1}
(this includes the complex characters Ĝ_c(D) of [R]).

**Theorem 2.** For 0 < δ ≤ 1/100, 1/2 + 3δ ≤ σ ≤ 9/10, |t| + 3 =: T ≤ X^{δ}/4:

  Σ_{D∈𝔇(X)} N^♮_D(σ;t) ≤ C(δ)X^{(13−8σ)/(7−2σ) + 60δ}.

At σ = 5/6 the exponent is 19/16 = 1.1875; the sum of the per-D bounds ([R2] Thm 1′) and the family bound of [R] Prop 1 both give
5/4 there. For comparison, 2 − σ (the density hypothesis on average) would be 7/6.

*Proof.* Let D ∈ 𝔇(X), D ≥ D₀(δ), V := X^{1/2}/2 (≤ |D|^{1/2}), Y := X^{5/(7−2σ)−5δ}, Z := 2Y log²Y (so 4Z ≤ X^{1−4δ}).
Each χ ∈ Ĝ(D)∖{1} has ≤ C₀log(X T²) zeros in the window ([C] (F4)); for each χ with a zero ρ_χ, σ ≤ β_χ ≤ 9/10, |γ_χ − t| ≤ 1,
Lemma 4 puts χ in 𝒳_I or 𝒳_II. Then N^♮_D(σ;t) ≤ C₀log(XT²)(|𝒳_I| + |𝒳_II|).

*Class II.* By Lemma 4(II), e^{−|u|} ≤ e^{1−|v−t|} for v = γ_χ + u, Cauchy–Schwarz at each v over all χ ∈ Ĝ(D)∖{1}, the mean
value Σ_{χ≠1}|L(½+δ+iv,χ)|² ≤ C(δ)|D|^{1/2+δ}(|v|+3)^{1+δ} ([D2] (2.1), mollifier length 1), Lemma 3(a),(c),(d):
|𝒳_II| ≤ C(δ)Y^{1/2+δ−σ}X^{1/2+2δ}T. Summed over D: ≤ C(δ)X^{3/2+2δ}Y^{1/2+δ−σ}T ≤ C(δ)X^{(13−8σ)/(7−2σ)+10δ}
(since 3/2 + (5/(7−2σ))(1/2 − σ) = (13−8σ)/(7−2σ)).

*Class I.* As in [C] Theorem 1 (with k = 1): each χ ∈ 𝒳_I has |F_j(χ,ρ_χ)| ≥ c₀/(2J) for some j, and, by Cauchy–Schwarz in τ
and |Φ̂_β(τ)| ≤ C_A(1+|τ|)^{−A}, |𝒳_{I,j}| ≤ C_AJ²∫(1+|τ|)^{−A}Σ_{χ∈Ĝ(D)}sup_{t'∈I_τ}|P_j(χ,t')|²dτ, I_τ := [t+τ−1, t+τ+1].
Sobolev: sup_I|P|² ≤ ∫_I((3/2)|P|² + |∂_{t'}P|²)dt'. Sum over D ∈ 𝔇(X). For |τ| ≤ X^{δ}/2 (so |t'| ≤ X^δ) apply Theorem 1 to
b_{j,t'} and ∂_{t'}b_{j,t'}; for |τ| > X^{δ}/2 use the per-D bound Σ_D Σ_χ|P_j|² ≤ CX^{1+δ}(N + X^{1/2})N^{1−2σ} ([C] Lemmas 2, 3) and
A := 1 + 3/δ: that part is ≤ CX^{1+δ}·X·N^{1−2σ}·X^{−3} ≤ N^{1−2σ}. With N = 2^j, V/4 ≤ N ≤ 2Z:

  Σ_D|𝒳_I(D)| ≤ C(δ)X^{41δ}Σ_{V/4≤2^j≤2Z}(X^{3/2} + 8^jX^{−1})2^{j(1−2σ)} ≤ C(δ)X^{42δ}(X^{3/2}V^{1−2σ} + Z^{4−2σ}X^{−1}).

Here X^{3/2}V^{1−2σ} ≤ 4X^{2−σ} ≤ 4X^{(13−8σ)/(7−2σ)} (the difference of the exponents is (2σ−1)(1−σ)/(7−2σ) ≥ 0), and
Z^{4−2σ}X^{−1} ≤ C X^{(5/(7−2σ))(4−2σ)−1} = CX^{(13−8σ)/(7−2σ)}. Adding the classes and the factor log(XT²) proves the theorem. □

## 5. Least primes in ring classes for almost all D < 0

**Theorem 3.** Let A₀ := (3 + √69)/10 = 1.1306623… . [R] Theorem 2 holds for every A > A₀: there are effective c, κ, δ > 0 and X₀
such that all D ∈ 𝔇 with |D| ≤ X, outside a set of size ≪_A X^{1−c}, satisfy (2.1) of [R] for every class C, every x ≥ |D|^A
and every H ∈ [x|D|^{−κ}, x]. In particular every primitive positive definite binary quadratic form of discriminant D represents
a prime ≤ 2|D|^A, and (Corollary) for all but O(N^{1−c}) squarefree n ≤ N with n ≡ 1, 2 (mod 4) the least prime a² + nb² is
≤ 2(4n)^A, and (x, x + x(4n)^{−κ}] contains such a prime for every x ≥ (4n)^A.

A₀ is where the two family bounds meet the Markov bar: with F₁(σ) := (13−8σ)/(7−2σ) (Theorem 2) and F₂(σ) := 10(1−σ)/(3−2σ)
([R] Prop 1), F_i(σ) < 1 + A(1−σ) ⟺ A > 6/(7−2σ) resp. A > 10/(3−2σ) − 1/(1−σ); these are equal at σ_c = 5 − √69/2 = 0.84669,
where both are (3+√69)/10 (σ_c is the smaller root of 4σ² − 40σ + 31 = 0). `thresh.py` checks this and the sup over a grid.

*Proof.* [R] §2 (with [R2]'s bookkeeping), with three changes; it suffices to treat A₀ < A ≤ 1.2.
1. *Parameters.* Put ε' := (1 − σ_c)(A − A₀)/8, δ := min{ε'/100, (A − 1)/80, A/4400}, κ := δ/4, ε₀ := ε₁ := δ/100,
   σ₀ := 1/2 + 4δ, σ_j := σ₀ + jε'/A (0 ≤ j ≤ J, J minimal with σ_J ≥ 1 − ε'/A); σ_g, T_X, Q₀ as in [R] §2.1.
2. *E₁ (Markov).* E₁(X) := {D ∈ 𝔇(X) : some j has N^c_D(σ_j, T_X) > X^{A(1−σ_j)−2ε'}}. For every σ ∈ [σ₀, 1),
   Σ_{D∈𝔇(X)}N^c_D(σ, T_X) ≤ C(δ)X^{F*(σ) + 61δ}, F*(σ) := min{F₁(σ), F₂(σ)} for σ ≤ 9/10 and F₂(σ) for σ > 9/10: indeed zeros with
   β ≤ 9/10 are counted by Theorem 2 (summed over the ≤ 2T_X unit windows), and those with β > 9/10 by [R] Prop 1 with
   F₂(9/10) = 0.834 < F₁(σ). The function φ(σ) := F*(σ) − A(1−σ) satisfies φ ≤ 1 − 8ε' on [σ₀, 1): on [σ₀, σ_c],
   F₁(σ) − A(1−σ) = 1 − (1−σ)(A − 6/(7−2σ)) ≤ 1 − (1−σ_c)(A − A₀); on [σ_c, 1), F₂ − A(1−σ) is decreasing (its derivative is
   A − 10/(3−2σ)² < 0) and equals 1 − (1−σ_c)(A − A₀) at σ_c. Markov then gives #E₁(X) ≤ C(δ)J X^{1−ε'} as in [R] §2.2 (the
   61δ and the T_X-powers are ≤ ε').
3. *Zeros with β < σ₀ ([R] §2.4(b)).* Each unit window holds ≤ C|D|^{1/2}log³(3|D|)log(|k|+3) zeros of all χ ∈ Ĝ(D) ([C] Thm 1(b)),
   and x^{β−1} ≤ x^{σ₀−1} ≤ |D|^{−A(1/2−4δ)}. Summing over |k| ≤ T₁: E_b ≤ H|D|^{1/2−A/2+5δ+κ+ε₁+o(1)} ≤ H|D|^{−(A−1)/2+6δ}.
Everything else ([R] §2.2 E₂, §2.3, §2.4(a),(c)–(f), §2.5) is unchanged; the conclusion needs ε' > δ/2, (A−1)/2 − 6δ > δ/2 and
0.005A − κ − ε₁ > δ/2, which hold by the choice of δ. The corollaries follow as in [R] §3 and [R2] Corollary 4′ (D = −4n,
principal class). □

Sectors ([R2] Theorem 3′, angular Hecke characters) and D > 0 ([R2] Theorem 4) are not improved: Theorem 1 is proved only for
ring class characters of imaginary quadratic orders (§6).

## 6. What is left, and the road to A₀ = 1

- *The only loss* is §3.3: pairs n, n' with a common factor w > G₂ ≈ X/N, where the T-intervals (length ≈ y²X/(w²√M)) are
  shorter than the modulus 4y²c; there Lemma 5 is used without cancellation. If Theorem 1 held with X^{3/2} alone for all
  N ≤ X^{1−δ}, Theorem 2 would give X^{2−σ+O(δ)} (the density hypothesis on average over D), and A₀ would drop to 1, the GRH
  exponent (F₁ replaced by 2 − σ meets F₂ at A = 1, σ = 0.8508; `thresh.py`). The pieces with w > G₂ are themselves family sums at
  scale X/w² with N/(cw) > X/w², i.e. the family large sieve beyond the conductor. Numerically (§7) the pole-killed family sum stays
  at the family size up to N = 1.5X. A proof would need an asymptotic for Σ_{n,n'}b(n)b(n')𝒩_c(nn'/w²; X/w²) in short
  T-intervals, i.e. cancellation in Salié-type sums Σ_{T² ≡ 4M (y²)}e(kT/y²) averaged over M = nn' (half-integral weight Kuznetsov).
- *Sensitivity.* With the cruder error term N² in Theorem 1 (instead of N³X^{−1}) the same argument gives A₀ = (√420 − 6)/12 = 1.2078;
  the bound min{N/(cw), √X/w + N/(c√X)} of Lemma 5(b) is what gives N³X^{−1}.
- *Sectors.* For angular characters ψ = χ·(β/|β|)^ℓ the dual count acquires the weight cos(ℓ arg γ); Poisson in T then sees a phase
  e(ℓθ(T)), which is harmless for ℓ ≪ X^{δ} but was not written out.
- *D > 0.* The dual count is over γ ∈ O_Δ modulo the unit group; (T, y) with T² − Δy² = 4M, Δ > 0, is an infinite set and needs a
  fundamental domain for the units, i.e. a count in a hyperbolic region; not attempted.

## 7. Numerical check of Theorem 1 (`famls.c`, `mkb.py`, `numerics1.md`)

𝒬(b) = Σ_{Δ∈𝒟, X/2≤|Δ|≤X}Φ(|Δ|/X)Σ_{C∈Cl(O_Δ)}|S_C|² (Φ a sin²-bump, all orders, classes from reduced forms, S_C from all
representations by the reduced form, ideals prime to f), σ = 5/6, t' = 0, ℓ = 2, V = √X/2, and n odd squarefree prime to 3 (instead of
prime to P). Rows X = 20000, 80000 were rerun on 2026-10-21 (`numerics2.txt`, also N = 2000, 8000, 16000 at X = 20000: molkill 0.122,
0.116, 0.118; rand 0.135, 0.135, 0.138).
Sequences: molkill = the b of §2 (weight with (Z)); molnokill = same without the ℓ-term (pole present); sqfkill / sqf = μ²(n)n^{−σ}
with / without the ℓ-term; prime = primes, no ℓ-term; rand = random signs. Ratio 𝒬(b)/(X‖b‖²):

| X | N | molkill | molnokill | sqfkill | sqf | prime | rand |
|---|---|---|---|---|---|---|---|
| 40000 | 500 | 0.125 | 0.136 | 0.118 | 0.438 | 0.294 | 0.125 |
| 40000 | 2000 | 0.123 | 0.246 | 0.121 | 1.399 | 0.712 | 0.130 |
| 40000 | 8000 | 0.121 | – | 0.132 | 5.263 | – | 0.133 |
| 40000 | 16000 | 0.117 | – | 0.112 | 10.20 | – | 0.133 |
| 20000 | 30000 | 0.115 | – | – | 27.07 | – | 0.148 |
| 80000 | 16000 | 0.120 | – | – | 7.377 | – | 0.132 |

Pole-carrying sequences grow like N/√X (the trivial-character term √X(Σ b)², ratio 𝒬/(√X(Σ|b|)²) ≈ 0.42 for sqf throughout);
pole-killed ones stay at the level of random signs, i.e. at the family size, even for N up to 1.5X (beyond what Theorem 1 proves).
Since max h(Δ) ≍ √X log X, 𝒬 ≈ 0.12X‖b‖² corresponds to Σ_{(D,χ)}|P|² ≈ X^{3/2}‖b‖².

## 8. Novelty audit (2026-10-21)
- arXiv API (`ax*.xml`, `axparse*.py`): 20 queries answered (class group L-functions × moments; large sieve × class group,
  ring class, dihedral, Hurwitz class numbers; zero density × class group / Hecke L-functions; least prime × quadratic forms × almost
  all); 22 more (`ax8`–`ax15`, `ax3_*`, `ax3_queries.txt`) were refused with "Rate exceeded". Relevant hits: Templier 0808.1476
  (twisted second moment, one D), Liyang Yang 2607.04476 (spectral reciprocity, uniform subconvexity), Baier 2610.15706
  (x₁² + x₂² − x₃² = D, spectral), Ditchen 1312.1502 (primes represented by forms on average over D; |D|^{20/3+ε} outside a
  density-zero set), Zaman 1710.08914, Gaudet 2408.01610, Thorner–Zaman 1909.01338, Chen An 2004.14410.
- Web search: Sardari 1803.03218 (one D: a positive proportion of classes represent a prime in [X, 2X]); Gaudet's thesis (Rutgers
  2023, `rutgers.txt`: a per-D large sieve for class group characters over almost-primes, Linnik-type D^L).
- GXL evidence API, 3 claims (`gxl3_*.json`): family large sieve at size X^{3/2} (one weak support, Ditchen; two unrelated); least
  prime in every class with A < 3/2 for almost all D (no evidence either way); zero density for class group L-functions on average
  over D (only one-level-density papers and Nordentoft, Wide moments of L-functions I, ANT 18 (2024)).
- None of these has a family large sieve at the family size, a density estimate on average over D below [R], or an almost-all
  exponent below 3/2 ([R2], ours).
- The ingredients are classical (Hecke relations, positivity, Poisson summation, Heath-Brown's quadratic large sieve, Jutila's mean
  value). New here: (i) the dual identity of Lemma 2 (all orders at once turn the dual of the family sieve into a plain lattice-point
  count with no class numbers and no squarefree condition); (ii) the pole-killing weight w(v) − ℓ^{2−2σ−2it'}w(ℓ²v), which removes
  the trivial-character main term from the class I sums of every D simultaneously; (iii) the squarefree-supported detector that keeps
  the class I coefficients independent of D.

## 9. Dead ends (NOGO candidates)
- Matching the main terms of the family sieve with the trivial and genus characters (Siegel–Weil on average, trace formula,
  [R2] §4 lead): not needed; with the weight (Z) all main terms are small and are bounded, not matched.
- Restricting the dual count to fundamental discriminants: the squarefree sieve on (4nn' − T²)/y² needs moduli up to X^{1/2}
  with O(1) error each; positivity over all orders removes the condition (Lemma 2).
- Standard ideal-based mollifier coefficients a_𝔞 = Σ_{𝔡|𝔞}μ_K(𝔡): they depend on D through non-primitive ideals, so the dual is
  no longer a function of N𝔞 alone; the squarefree series 𝒢 (Lemma 3) avoids this.
- Pole killing by smoothing the mollifier (Selberg weights) or by dyadic differences: only log savings, or destroys the head of the
  detector; the rational ideal (ℓ) is the only universal choice (F(s) = 1 − ℓ^{2−2s}).
- Per-D improvements below e(σ) of [R2]: unchanged NOGO ([R2] §7).
