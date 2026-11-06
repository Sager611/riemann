# Explicit polylogarithmic least prime in every Chebotarev class: derivation

Everything below is classical analytic number theory except one input: a zero-free half-plane for the Dedekind
zeta function of the top field. Every numerical constant comes from `cert.py` (Arb ball arithmetic), whose log is
`logs/cert_output.txt`. Section numbers (§) are referred to in the comments of `cert.py`.

## 0. Statements

**Hypothesis H(δ)** (for a number field L): ζ_L(s) ≠ 0 for Re s > 1 − δ.
- H(1/32) holds for every L: Cureton–Thorner, Hexagon 2610.00200v2 (l. 172/179), Lean development
  github.com/AxiomMath/DedekindZetaZeroFree (confirmed by Zenodo 23373237).
- H(δ_2), δ_2 := 1 − B_2, B_2 = (777 − √1081)/852 = 0.87338196661320752, holds for every L if Liu, Hexagon
  2611.00082 (unrefereed; our audit Riemann/Deep/Hexagon2/REPORT.md) is right (his theorem for the trivial character
  of L). 1/δ_2 = 7.8977691664599….

Only ζ_L of the top field is used.

**Theorem A (from H(1/32)).** Let L/K be a Galois extension of number fields, L ≠ ℚ, G = Gal(L/K), and C ⊂ G a
conjugacy class. There is a prime 𝔭 of K, unramified in L, of degree one over ℚ (N𝔭 = p prime), with Frob_𝔭 = C and

  N𝔭 ≤ (0.484 · log(3Δ_L))^32,  Δ_L := |disc L|.

Moreover N𝔭 ≤ (c·log 3Δ_L)^32 with c = 0.29 if log Δ_L ≥ 10, c = 0.151 if log Δ_L ≥ 100, c = 0.108 if
log Δ_L ≥ 1000, and c = 0.094 if log Δ_L ≥ 10⁴.

**Theorem B (from H(δ_2), i.e. conditional on Liu).** Same setting:

  N𝔭 ≤ (0.62 · log(3Δ_L))^{1/δ_2} ≤ (0.62 · log(3Δ_L))^7.90,

and c' = 0.595, 0.426, 0.366, 0.344 for log Δ_L ≥ 10, 100, 1000, 10⁴. The last inequality holds because
0.62·log 9 > 1.

**Optimality of the shape at small fields.** For L = ℚ(√−3), K = ℚ, C = {1}, the least split prime is 7. So no
constant below 7^{1/32}/log 9 = 0.48365 (resp. 7^{δ_2}/log 9 = 0.58228) can work in a bound of this shape.
Theorem A's 0.484 is within 0.07 % of that, and Theorem B's 0.62 within 6.5 %.

## 1. Set-up and two classical lemmas

Fix g ∈ C. Put H := ⟨g⟩ (cyclic), E := L^H, n_L = [L:ℚ] = r_1 + 2r_2, Λ := log Δ_L, ℓ := Λ + log 3.
By Artin reciprocity each character χ of H = Gal(L/E) is a primitive finite-order Hecke character of E, and the
Artin L-function of χ equals the primitive Hecke L-function L(s, χ), including the local factors at ramified
primes (the Artin and Hecke conductors agree). The Artin formalism for the regular representation of H gives

  ζ_L(s) = ∏_{χ ∈ Ĥ} L(s, χ)  exactly.   (1.1)

Here L(s, 1) = ζ_E(s), and L(s, χ) is entire for χ ≠ 1.

**Lemma 1.1 (Minkowski).** For L ≠ ℚ: Λ ≥ n_L log(π/2), i.e. n_L ≤ Λ/0.451582….
*Proof.* Minkowski's bound gives √Δ_L ≥ (π/4)^{r_2} n^n/n!, with n = n_L. Since 2r_2 ≤ n, Δ_L ≥ (π/4)^n (n^n/n!)².
It suffices that n^n/n! ≥ 2^{n/2}. Equality holds at n = 2. Going from n to n + 1, the left side is multiplied by
(1 + 1/n)^n ≥ 2 and the right side by √2. ∎

**Lemma 1.2 (small fields).** If Δ_L ≤ 12 then n_L = 2, and L is one of ℚ(√d), d ∈ {−3, −4, 5, −7, 8, −8, −11, 12}.
*Proof.* Minkowski: for n = 3, r_2 = 1 gives √Δ ≥ (π/4)(27/6) = 3.534, so Δ ≥ 12.49; r_2 = 0 gives Δ ≥ 20.25. For
n = 4, √Δ ≥ (π/4)²(256/24) = 6.58, so Δ ≥ 43.3. For n ≥ 5 the bound (π/4)^{n/2}n^n/n! keeps increasing, because the
step ratio (π/4)^{1/2}(1 + 1/n)^n ≥ 0.886·2.44 > 1. The list is the fundamental discriminants with |d| ≤ 12. ∎

**Lemma 1.3 (direct range for Theorem B).** If Δ_L ≤ 257 then n_L ≤ 4.
*Proof.* For n = 5, Minkowski gives √Δ ≥ (π/4)^{r_2}·5^5/5! ≥ (π/4)²·26.04 = 16.064, so Δ ≥ 258.05.
For n = 6, √Δ ≥ (π/4)³·6^6/6! = 31.39. Larger n gives larger bounds, as in Lemma 1.2. ∎

The transitive groups of degree ≤ 4 are C2; C3, S3; C4, V4, D4, A4, S4. PARI/GP 2.19 `nflist(G, [1, 257])` lists
all fields with these groups and |disc| ≤ 257. It finds 157 + 3 + 32 + 1 + 3 + 2 + 0 + 2 = 200 fields. The quartic
ones are 117, 189 (D4), 125 (C4), 144, 225, 256 (V4), 229, 257 (S4), as in the classical tables.

**Zeros.** Assume H(δ) with 0 < δ ≤ 1/2. The nontrivial zeros of ζ_L are those of the completed
ξ_L(s) = s(s − 1)Δ_L^{s/2}γ_L(s)ζ_L(s), where γ_L = Γ_ℝ(s)^{r_1}Γ_ℂ(s)^{r_2}, Γ_ℝ(s) = π^{−s/2}Γ(s/2) and
Γ_ℂ(s) = 2(2π)^{−s}Γ(s). Since ξ_L(s) = ξ_L(1 − s), every nontrivial zero ρ = β + iγ satisfies δ ≤ β ≤ 1 − δ.
By (1.1), the multiset union over χ ∈ Ĥ of the nontrivial zeros of L(s, χ) is the multiset of nontrivial zeros of ζ_L.
The same holds for the trivial zeros. Write m_j(χ) for the order of the zero of L(s, χ) at s = −j. Then

  Σ_χ m_j(χ) = m_j(ζ_L) = r_1 + r_2 − 1 (j = 0), r_1 + r_2 (j ≥ 2 even), r_2 (j odd),

which is ≤ n_L − 1 for j = 0 and ≤ n_L for every j. These are the orders of the poles of γ_L, minus the simple pole of
ξ_L/(s(s−1)) at s = 0.

## 2. The kernel

Parameters: h > 0, s* ∈ ℝ, σ_0 > 1. Put U := 3h and B_3 := 1_{[0,h]} ∗ 1_{[0,h]} ∗ 1_{[0,h]}, so that
B_3(u) = h²M_3(u/h), where M_3(v) = v²/2 on [0,1], (−2v² + 6v − 3)/2 on [1,2] and (3 − v)²/2 on [2,3]. Define

  F(u) := e^{s* u} B_3(u),  k(t) := F(log(1/t)) for 0 < t ≤ 1, k(t) := 0 for t > 1,
  Φ(s) := ∫_0^∞ k(t) t^{s−1} dt = ∫_0^U F(u) e^{−su} du = ((1 − e^{−h(s−s*)})/(s − s*))³.   (2.1)

The closed form holds because the Laplace transform of a convolution is the product of the transforms. Properties:
- (K1) F ≥ 0, F ∈ C¹ (B_3 is C¹), supp F = [0, U]; k is continuous on (0, ∞) with support [e^{−U}, 1].
  F is log-concave (B_3 is), hence unimodal.
- (K2) Φ is entire, and Φ(s̄) is the conjugate of Φ(s). For real σ, Φ(σ) = ∫F e^{−σu} > 0 is decreasing in σ, and
  Φ(−j) ≤ e^{jU}Φ(0).
- (K3) For γ ≠ 0: |Φ(σ + iγ)| ≤ (1 + e^{h(s*−σ)})³/|γ|³, since |1 − e^{−hz}| ≤ 1 + e^{−h Re z} and |z| ≥ |γ|.
- ‖k‖ := max F; S_j := max_u F(u)e^{−u/j}. The certificate encloses these by Arb balls on the three
  polynomial pieces, with no cancellation.

In `cert.py`, Φ is evaluated as (h e^{−w} sinc(iw))³ with w = h(s − s*)/2, which is the same function.

## 3. Analytic side

### 3.1 Explicit formula
**Lemma 3.1.** Let χ ∈ Ĥ and x > e^U. Then

  S_χ(x) := Σ_𝔞 Λ_E(𝔞)χ(𝔞)k(N𝔞/x) = 1_{χ=1}·xΦ(1) − Σ_ρ x^ρ Φ(ρ) − Σ_{j≥0} m_j(χ) x^{−j} Φ(−j).

Here 𝔞 runs over the integral ideals of E, and χ(𝔞) = 0 when 𝔞 is not coprime to the conductor. ρ runs over the
nontrivial zeros of L(s, χ) with multiplicity. Both sums converge absolutely.

*Proof.* By (K1) and (K3), Mellin inversion k(t) = (2πi)^{−1}∫_{(c)}Φ(s)t^{−s}ds holds for every real c and t > 0,
with absolute convergence. With c = 2, summing against Λ_E(𝔞)χ(𝔞) is justified by absolute convergence and gives

  S_χ(x) = (2πi)^{−1} ∫_{(2)} (−L'/L)(s, χ) x^s Φ(s) ds.

Move the line to Re s = −J − 1/2, J ∈ ℕ, through rectangles with horizontal sides at ±T_k → ∞. Choose T_k in
[k, k + 1] so that both +T_k and −T_k are at distance ≫ 1/log(A(χ)T_k^{n_E}) from every ordinate. (For complex χ
the ordinates need not be symmetric.) This is possible because the number of zeros with |γ ∓ T| ≤ 1 is
O(log A(χ) + n_E log(|T| + 2)), Iwaniec–Kowalski Prop. 5.7. Then |L'/L(σ ± iT_k)| ≪_χ log²T_k for −1 ≤ σ ≤ 2, by the partial-fraction formula IK Prop. 5.7(3); Lagarias–Odlyzko
1977 §5 do the same for Dedekind and Hecke L-functions. For σ ≤ −1, the functional equation
L'/L(s, χ) = −log A(χ) − γ_χ'/γ_χ(s) − γ_χ̄'/γ_χ̄(1 − s) − L'/L(1 − s, χ̄) gives |L'/L(s, χ)| ≪_χ log(|s| + 2) at
distance ≥ 1/4 from the nonpositive integers. By (K3) the horizontal integrals are O_{χ,J}(T_k^{−3}log²T_k) → 0.

On Re s = −J − 1/2 we have |x^sΦ(s)| ≤ x^{−J−1/2} min(Φ(−J−1/2), (1 + e^{h(s*+J+1/2)})³/|t|³)
≪ e^{|s*|U}(e^U/x)^{J+1/2}min(1, |t|^{−3}). The left integral is therefore O_χ((e^U/x)^{J}·log(J + 2)) → 0 as J → ∞,
because x > e^U.

Φ is entire, so the residues are: xΦ(1) at the pole s = 1 (χ = 1 only); −m x^ρΦ(ρ) at each nontrivial zero of
multiplicity m; −m_j(χ)x^{−j}Φ(−j) at s = −j. The zero sum converges absolutely, because the number of zeros with
|γ| ≤ T is O(T log T) and (K3) holds. The trivial series is ≤ n_EΦ(0)Σ_j(e^U/x)^j by (K2). ∎

### 3.2 The Deuring sum
Put θ(𝔞) := |H|^{−1}Σ_{χ∈Ĥ} χ̄(g)χ(𝔞) and S(x) := Σ_𝔞 Λ_E(𝔞)θ(𝔞)k(N𝔞/x) = |H|^{−1}Σ_χ χ̄(g)S_χ(x).

**Lemma 3.2.** For a prime 𝔓 of E with inertia group I ⊂ H and Frobenius Frob_𝔓 ∈ H/I in L/E:

  θ(𝔓^m) = |I|^{−1}·1[g ∈ Frob_𝔓^m I] ∈ [0, 1].

In particular θ(𝔓^m) = 1[Frob_𝔓^m = g] when 𝔓 is unramified in L/E. (H is abelian, so I and Frob do not depend on
the prime of L above 𝔓.)
*Proof.* χ(𝔓) = χ(Frob_𝔓) if χ|_I = 1 (𝔓 ∤ 𝔣_χ), and χ(𝔓) = 0 otherwise. So Σ_χ χ̄(g)χ(𝔓)^m is the sum over the
characters of H/I, which equals |H/I|·1[g ≡ Frob^m mod I] by orthogonality. ∎

Combining Lemma 3.1 with (1.1) and §1 (the zeros of all L(s, χ) together are those of ζ_L), for every x > e^U
(assumed from here on, and asserted by cert.py):

  |H|·S(x) ≥ xΦ(1) − Σ_{ρ: ζ_L(ρ)=0} x^β|Φ(ρ)| − T(x),
  T(x) := Σ_j m_j(ζ_L)x^{−j}Φ(−j) ≤ (n_L − 1)Φ(0) + n_LΦ(0)·τ_x,
  τ_x := Σ_{j≥1} x^{−j}Φ(−j)/Φ(0) ≤ ε_x := e^U/(x − e^U).   (3.1)

cert.py bounds τ_x with the closed form Φ(−j) = ((e^{h(j+s*)} − 1)/(j + s*))³ for 1 ≤ j ≤ 30. For j > 30 it uses
Φ(−j) ≤ e^{jU}Φ(0) and sums the geometric tail. τ_x decreases in x.

### 3.3 The zero sum: Hadamard identity at a real point
For real σ > 1 (Stark 1974, Lemma 1; Lagarias–Odlyzko 1977, (5.1)–(5.2)):

  Σ_ρ Re 1/(σ − ρ) = ½Λ + 1/σ + 1/(σ − 1) + γ_L'/γ_L(σ) + ζ_L'/ζ_L(σ).   (3.2)

The sum runs over the nontrivial zeros of ζ_L, and every term (σ − β)/|σ − ρ|² is > 0.
*Derivation.* ξ_L has order 1, so ξ_L(s) = e^{A+Bs}∏_ρ(1 − s/ρ)e^{s/ρ}. The functional equation ξ_L(s) = ξ_L(1 − s)
gives Re B = −Σ_ρ Re(1/ρ), hence Re ξ_L'/ξ_L(s) = Σ_ρ Re 1/(s − ρ). Finally
ξ_L'/ξ_L = 1/s + 1/(s − 1) + ½Λ + γ_L'/γ_L + ζ_L'/ζ_L. (Sanity check 4 confirms the normalisation numerically for ζ.)

Now drop ζ_L'/ζ_L(σ) = −Σ Λ_L(𝔞)N𝔞^{−σ} ≤ 0. Then
γ_L'/γ_L(σ) = −r_1 a(σ) − r_2 b(σ), with a = (log π − ψ(σ/2))/2 and b = log 2π − ψ(σ). Since
r_1a + r_2b = r_1a + (2r_2)(b/2) ≥ n_L·g(σ), where g(σ) := min(a(σ), b(σ)/2),

  Σ_ρ w(ρ) ≤ B(σ_0) := ½Λ + 1/σ_0 + 1/(σ_0 − 1) − n_L g(σ_0),  w(ρ) := Re 1/(σ_0 − ρ).   (3.3)

(For the 405 kernels used, 1.0376 ≤ σ_0 ≤ 1.7184, so 0.8 < g(σ_0) ≤ 1.18. The certificate uses a rigorous lower bound for g.)

### 3.4 One zero against its Hadamard weight: maximum principle
Let 𝒮 := {δ ≤ Re s ≤ 1 − δ} and R(s) := |x^{s−(1−δ)}Φ(s)| / Re(1/(σ_0 − s)) on 𝒮.
The denominator is a positive harmonic function on 𝒮, because σ_0 > 1. So log R = Re((s − 1 + δ)log x) +
log|Φ(s)| − log Re(1/(σ_0 − s)) is subharmonic: log|Φ| is subharmonic, and −log of a positive harmonic function is
subharmonic. R is continuous on the closed strip, and R(s) ≤ C|Im s|^{−1} → 0 as |Im s| → ∞, by (K3) and
1/Re(1/(σ_0 − s)) = |σ_0 − s|²/(σ_0 − β). By the maximum principle on rectangles, sup_𝒮 R equals the sup on the two
boundary lines. With (K2) (symmetry γ ↦ −γ):

  x^{β−(1−δ)}|Φ(ρ)| ≤ M_x·w(ρ) for every ρ ∈ 𝒮,  M_x := max(A_1, x^{2δ−1}A_2),
  A_1 := sup_{γ≥0} |Φ(1−δ+iγ)|·((σ_0−1+δ)² + γ²)/(σ_0−1+δ),  A_2 := sup_{γ≥0} |Φ(δ+iγ)|·((σ_0−δ)² + γ²)/(σ_0−δ).

With (3.3): Σ_ρ x^β|Φ(ρ)| ≤ x^{1−δ}M_x B(σ_0).   (3.4)

**Certified sups.** cert.py bounds A_1 and A_2 by branch and bound. It starts with 256 γ-intervals on [0, G] and
evaluates R on each complex ball β + i[γ_a, γ_b]. It bisects every interval whose upper bound exceeds the best point
value by more than a factor 1 + tol (tol = 10⁻³ for A_1, 0.5 for A_2) and keeps the maximum of all final upper
bounds. For γ ≥ G it uses the tail bound (1 + e^{h(s*−β)})³((σ_0−β)² + G²)/((σ_0−β)G³), which follows from (K3) and is
decreasing in γ. G is doubled until the tail is below a quarter of R(β + 0i).

## 4. Arithmetic side

Split the terms (𝔓, m) of S(x) (𝔞 = 𝔓^m) into three classes. Let 𝔭 := 𝔓 ∩ O_K and N𝔓 = p^f.
- **Good:** m = 1, f = 1, 𝔭 unramified in L. Here θ(𝔓) = 1[Frob_𝔓 = g]. Write Good(x) := Σ_{good, Frob_𝔓 = g} log p·k(p/x).
- **Ram:** 𝔭 ramified in L. This includes every 𝔓 ramified in L/E.
- **Pow:** 𝔭 unramified in L, and mf ≥ 2, i.e. N𝔓^m = p^j with j ≥ 2.

So S(x) = Good(x) + Ram(x) + Pow(x), and Ram, Pow ∈ [0, ∞) by Lemma 3.2.

### 4.1 Deuring's step: a good term gives the prime we want
Let 𝔓 be good with Frob_𝔓(L/E) = g and N𝔓 = p ≤ x. Then f(𝔭/p) divides f(𝔓/p) = 1, so N𝔭 = p and 𝔭 has degree one.
Let 𝔔 be a prime of L above 𝔓. Frob_𝔔(L/E) = g means g ∈ D(𝔔/𝔓) ⊂ D(𝔔/𝔭) and g(a) ≡ a^{N𝔓} mod 𝔔 for
a ∈ O_L. Since 𝔭 is unramified in L, Frob_𝔔(L/K) is the unique element of D(𝔔/𝔭) acting as a ↦ a^{N𝔭}. Since
N𝔭 = N𝔓, it is g. Hence Frob_𝔭 = C, 𝔭 is unramified, and N𝔭 ≤ x.
So **Good(x) > 0 implies the theorem with the bound x**: k(p/x) > 0 forces p ≤ x.

### 4.2 Ramified terms
Let 𝔭 be ramified in L, with e ≥ 2 and residue degree f in L/K. Every prime above 𝔭 has different exponent ≥ e − 1, so
v_𝔭(d_{L/K}) ≥ (|G|/(ef))·f·(e − 1) ≥ |G|/2. The number of primes of E above 𝔭 is ≤ [E:K] = |G|/|H|.
Fix 𝔓 | 𝔭 and put Δ := log N𝔓. Three bounds hold for Σ_{m≥1} Δ·k(N𝔓^m/x) = Σ_m Δ·F(log x − mΔ):
- ≤ ‖k‖ log x: at most log x/Δ values of m have N𝔓^m ≤ x;
- ≤ ‖k‖(U + Δ): at most U/Δ + 1 values of m have N𝔓^m in the support [xe^{−U}, x];
- ≤ Φ(0) + ‖k‖Δ: a Riemann sum of the unimodal F with step Δ exceeds ∫F = Φ(0) by at most Δ·max F.

Sum over 𝔓 | 𝔭, using #{𝔓 | 𝔭} ≤ [E:K], Σ_{𝔓|𝔭} f(𝔓/𝔭) ≤ [E:K], and Σ min ≤ min Σ. This gives
≤ [E:K]·m(N𝔭) with m(y) := min(‖k‖ log x, ‖k‖(U + log y), Φ(0) + ‖k‖ log y). Each of the three
functions m_i(y)/log y decreases in y ≥ 2, so m(y)/log y ≤ m(2)/log 2 =: R_f. Then

  |H|·Ram(x) ≤ |H|·Σ_{𝔭 ram} (2/|G|)v_𝔭(d_{L/K})·(|G|/|H|)·m(N𝔭)
           ≤ 2·Σ_𝔭 v_𝔭(d_{L/K}) log N𝔭 · R_f ≤ 2Λ·R_f,
  R_f = min(‖k‖ log x, ‖k‖(U + log 2), Φ(0) + ‖k‖ log 2)/log 2.   (4.1)

Here we used Σ_𝔭 v_𝔭(d_{L/K}) log N𝔭 = log N(d_{L/K}) = Λ − [L:K] log Δ_K ≤ Λ.

### 4.3 Prime powers and higher-degree primes
For a rational prime p and j ≥ 2, the pairs (𝔓 | p, m) with mf_𝔓 = j contribute
Σ_{𝔓|p, f_𝔓|j} f_𝔓 log p·k(p^j/x) ≤ n_E log p·k(p^j/x), because Σ_{𝔓|p} f_𝔓 ≤ n_E and θ ≤ 1. With n_E|H| = n_L:

  |H|·Pow(x) ≤ n_L P_k(x),  P_k(x) := Σ_{j≥2} Σ_p log p·k(p^j/x).   (4.2)

- **x ≤ 2·10¹²:** P_k(x) is summed exactly over the primes p ≤ √x ≤ 1.42·10⁶ (sieve). Each k-value is enclosed in
  Arb, and only p^j in the support [xe^{−U}, x] count.
- **x > 2·10¹²:** fix j, put φ(y) := k(y^j/x). φ is C¹, compactly supported in (1, ∞) (because x > e^U), and
  unimodal by (K1).
  Partial summation gives Σ_p log p·φ(p) = −∫θ(y)φ'(y)dy ≤ ∫θ(y)(−φ'(y))_+ dy ≤ 1.01624∫y(−φ'(y))_+ dy, using
  θ(y) < 1.01624y for y > 0 (Rosser–Schoenfeld 1962, Thm 9). For unimodal φ with maximum at y_m this equals
  1.01624(y_mφ(y_m) + ∫_{y_m}^∞ φ) ≤ 1.01624(sup_y yφ(y) + ∫φ). Substituting y = (xe^{−u})^{1/j}:
  sup_y yφ(y) = x^{1/j}S_j and ∫φ dy = (x^{1/j}/j)Φ(1/j). For j > J := ⌊log x/log 2⌋ there is no prime
  p ≤ x^{1/j}. For j ≥ 4 use S_j ≤ ‖k‖ and Φ(1/j) ≤ Φ(0):

  P_k(x) ≤ 1.01624·[x^{1/2}(S_2 + Φ(½)/2) + x^{1/3}(S_3 + Φ(⅓)/3) + (J − 3)x^{1/4}(‖k‖ + Φ(0)/4)].   (4.3)

## 5. The inequality, and how the certificate covers all fields

### 5.1 The master inequality
From (3.1), (3.4), (4.1) and (4.2), for every L/K, C, every kernel and every x > e^U:

  |H|·Good(x)/x ≥ LHS(x; Λ, n_L) := Φ(1) − x^{−δ}M_x B(σ_0)
                 − x^{−1}[(n_L−1)Φ(0) + n_LΦ(0)τ_x + 2Λ·R_f + n_L P_k(x)].   (5.1)

If LHS > 0 then, by §4.1, a prime as in the theorem exists with N𝔭 ≤ x. For fixed x and kernel, LHS is affine in
n_L and decreasing in Λ. Λ enters only through ½Λ in B(σ_0) and the ramified term.

### 5.2 Interval checks
A chain is a constant c, a start Λ_0 and intervals [Λ_a, Λ_b] (contiguous: each Λ_a ≤ the previous Λ_b), each with
its own kernel (h, s*, σ_0). For an interval, put x := (c(Λ_a + log 3))^{1/δ} and require

  LHS(x; Λ_b, 2) > 0 and LHS(x; Λ_b, Λ_b/log(π/2)) > 0.

Let L have Λ ∈ [Λ_a, Λ_b]. By Lemma 1.1, n_L ∈ [2, Λ_b/log(π/2)], so by affinity and monotonicity
LHS(x; Λ, n_L) > 0. Hence N𝔭 ≤ x ≤ (c·log 3Δ_L)^{1/δ}.

*Rounding and the sign of B* (raised by the reviewer). The certificate replaces Φ(1) by a lower bound P and
x^{−δ}M_x by an upper bound q, as exact numbers. It then evaluates the affine function
f(n) := P − q·B(σ_0; Λ_b, n) − x^{−1}[(n−1)Φ(0)⁺ + nΦ(0)⁺τ_x⁺ + 2Λ_bR_f⁺ + nP_k⁺] (⁺ = upper bound) at the two
endpoints, and min f over [2, n_max] is attained at an endpoint. At the artificial endpoint n_max, B can be negative;
that does not matter. For an actual field, B(σ_0; Λ_b, n_L) ≥ B(σ_0; Λ, n_L) ≥ Σ_ρ w(ρ) ≥ 0 by (3.3), so
|H|·Good(x)/x ≥ Φ(1) − x^{−δ}M_x·Σ_ρ w − (…) ≥ f(n_L) ≥ min(f(2), f(n_max)) > 0. Upper-rounding q is safe exactly
because it multiplies the nonnegative quantity B(σ_0; Λ_b, n_L) of the actual field.

### 5.3 Tail check
The last element of a chain is a kernel valid for all Λ ≥ Λ_N with x := (cℓ)^{1/δ}, ℓ = Λ + log 3, which grows with
Λ. Use n_L ≥ 2 in B (needs g(σ_0) > 0, asserted) and n_L ≤ Λ/log(π/2) in the remaining terms:

  x^{−δ}B ≤ (½ℓ − ½log 3 + 1/σ_0 + 1/(σ_0−1) − 2g)/(cℓ) ≤ 1/(2c) + max(0, 1/σ_0 + 1/(σ_0−1) − ½log 3 − 2g)/(cℓ_N).

Also M_x ≤ M_{x_N}, τ_x ≤ τ_{x_N}, and R_f ≤ ‖k‖ log x/log 2. Each remaining term is a multiple of Λ/x, Λ log x/x, Λx^{−1/2}, Λx^{−2/3} or
Λ log x·x^{−3/4} (by (4.3), with J − 3 ≤ 2 log x). Along x = (cℓ)^{1/δ}, Λx^{−a} decreases in ℓ iff
δℓ/(ℓ − log 3) < a, and Λ log x·x^{−a} decreases if δ(ℓ/(ℓ − log 3) + 1/log(cℓ)) < a. The certificate asserts
(1 − δ)ℓ_N > log 3, 2δℓ_N/(ℓ_N − log 3) < 1 and δ(ℓ_N/(ℓ_N − log 3) + 1/log(cℓ_N)) < 3/4. These conditions are
monotone in ℓ, so they hold for all ℓ ≥ ℓ_N. So the lower bound for LHS evaluated at Λ_N is valid for every
Λ ≥ Λ_N.   (5.3)

### 5.4 Direct ranges
**Theorem A: Δ_L ≤ 12.** By Lemma 1.2, L is quadratic and K ∈ {ℚ, L}.
- K = ℚ: the classes are "split" and "inert".
- K = L: every prime is unramified in L/L, so the least degree-one prime of L is min(least split, least ramified).

cert.py computes these with Kronecker symbols:

| d | −3 | −4 | 5 | −7 | 8 | −8 | −11 | 12 |
|---|---|---|---|---|---|---|---|---|
| split, inert, degree-one | 7, 2, 3 | 5, 3, 2 | 11, 2, 5 | 2, 3, 2 | 7, 3, 2 | 3, 5, 2 | 3, 2, 3 | 11, 5, 2 |

The worst ratio N^{1/32}/log(3Δ_L) is 0.48365, at d = −3 and N = 7. It is below 0.484 (asserted in cert.py).

**Theorem B: Δ_L ≤ 257.** By Lemma 1.3 there are 200 fields; `direct.gp` (PARI) handles every subfield K with L/K
Galois:
- K = L: the least degree-one prime of L.
- K = ℚ when L/ℚ is Galois (C2, C3, C4, V4): the least unramified p for each Frobenius element, via
  `idealfrobenius`.
- K = k, each quadratic subfield of a quartic L (C4, V4, D4): the least degree-one prime of k that splits, resp.
  stays inert, in L/k, via `rnfidealprimedec`. A prime counts as inert only if its absolute ramification index is 1.
  This may skip some inert primes, so the reported p is a valid witness and an upper bound for the least one.

That gives 563 (L, K, C) cases, in logs/direct_range.txt. cert.py checks the field counts per group, the number of
classes per field, and N ≤ (0.62 log 3Δ_L)^{1/δ_2} for every case. The worst ratio N^{δ_2}/log(3Δ_L) is again
ℚ(√−3) with N = 7: 0.58228. The next worst are 0.523 (ℚ(√−3), K = L, N = 3) and 0.500 (ℚ(√5), split, N = 11).
The completeness of nflist for these groups is the one computational input taken from PARI.

### 5.5 What the certificate checks, and rounding directions
- Coverage: the uniform chain for A starts at Λ_0 = 2.5649 < log 13, and the one for B at 5.5529 < log 258. Both
  are asserted in Arb. (The float log 13 is 2.6·10⁻¹⁸ too large, which is why 2.5649 is used.)
- Lower bound for Φ(1). Upper bounds for Φ(0), Φ(½), Φ(⅓), ‖k‖, S_2, S_3, A_1, A_2, P_k(x), 1/σ_0 + 1/(σ_0 − 1),
  x^{−δ}, τ_x and R_f. Lower bound for g(σ_0) (Arb digamma).
- B_2 = (777 − √1081)/852 and 1/32 enter as Arb balls. Kernel parameters are the exact binary values of the floats in
  `chains.json`. `design.py` (floating point, not part of the proof) produced the chains.
- The margins LHS/Φ(1) reported per interval in logs/cert_output.txt are all ≥ 0.0029.

| chain | δ | c | range of log Δ_L | intervals + tail | min LHS/Φ(1) |
|---|---|---|---|---|---|
| CTu | 1/32 | 0.484 | ≥ 2.5649 (and Δ_L ≤ 12 direct) | 2 + tail | 0.00300 |
| CT10 / CT100 / CT1000 / CT10000 | 1/32 | 0.29 / 0.151 / 0.108 / 0.094 | ≥ 10 / 100 / 1000 / 10⁴ | 1 / 5 / 12 / 19 + tail | 0.0164 / 0.0105 / 0.0066 / 0.0034 |
| Lu | δ_2 | 0.62 | ≥ 5.5529 (and Δ_L ≤ 257 direct) | 20 + tail | 0.00300 |
| L10 / L100 / L1000 / L10000 | δ_2 | 0.595 / 0.426 / 0.366 / 0.344 | ≥ 10 / 100 / 1000 / 10⁴ | 10 / 19 / 58 / 249 + tail | 0.0032 / 0.0139 / 0.0039 / 0.0030 |

The command `uv run --no-project --with python-flint python3 cert.py chains.json` runs in 105 s on one core, with
peak RSS 38 MB (logs/cert_time.txt). It must be run from this directory, because it reads logs/direct_range.txt.

An independent floating-point re-evaluation (`xcheck.py`) uses mpmath quadrature for Φ, dense grids for the sups and
a direct prime-power sum. For the first interval of each uniform chain it gives LHS/Φ(1) = 0.00415 (A) and
0.00438 (B). The rigorous values are 0.00301 and 0.00303; they are lower because the branch and bound
overestimates A_1 by up to 10⁻³. On a dense grid xcheck.py also confirms that the interior maximum of R does not
exceed the boundary maximum (§3.4).

## 6. Remarks
- **Limit of the method.** With this kernel family, inf A_1/(2Φ(1)) over (h, s*, σ_0 → 1⁺) is 2.5958·δ for both δ
  (scale invariance; `.scratch` floating-point optimisation, not certified): 0.0811 for δ = 1/32 and 0.3287 for δ_2.
  So the method gives N𝔭 ≤ ((2.596δ + o(1))·log Δ_L)^{1/δ} as Δ_L → ∞. The certified table reaches 0.094 (3.0δ) and
  0.344 (2.7δ).
- **Shape check under GRH.** The same argument with δ = ½ would give (κ·log Δ_L)² with κ → 1.30. Bach–Sorenson's GRH
  bound (4 log Δ_L + 2.5n_L + 5)² has the same shape.
- **n_L uniformity.** n_L enters only through the Γ-factor saving −n_L g(σ_0) (which helps), the trivial zeros
  (n_L − 1)Φ(0) + n_LΦ(0)τ_x, and the prime powers n_L P_k(x). All are covered by checking both ends n_L = 2 and
  n_L = Λ/log(π/2) (Lemma 1.1).
- **No GRH and no zero counting.** The zeros enter only through (3.2) at one real point σ_0 > 1, together with the
  zero-free strip. There are no N(T) bounds and no log-free density estimates.

## 7. Self-review (adversarial)
1. *Deuring step:* §4.1. A degree-one prime 𝔓 of E over an unramified 𝔭 has f(𝔓/𝔭) = 1, so the two Frobenius
   elements coincide. The case K = L (H trivial, E = L) is included, and there every prime of L is unramified in L/K.
2. *Ramified primes of L/E:* Lemma 3.2 shows θ ∈ [0, ½] there, and they lie in "Ram" (§4.2) because
   e(𝔔/𝔭) = e(𝔔/𝔓)e(𝔓/𝔭) ≥ 2.
3. *Primes of E above 𝔭 ramified in L/K but unramified in L/E:* also in "Ram". The bound (4.1) counts all primes of E
   above each ramified 𝔭.
4. *Truncations:* there are none. The explicit formula is an identity (Lemma 3.1). The trivial-zero series is summed
   with τ_x. The zero sum is bounded globally by (3.4). The B&B sups include an explicit tail.
5. *Uniformity in n_L:* §6; the zero sum has no n_L term at all.
6. *Downward roundings:* none. Every inequality in cert.py uses Arb upper or lower endpoints in the safe direction.
   Floats appear only as exact parameters and in loop limits. Loop limits are padded: prime ranges ±1, and the
   support test keeps any p^j whose ball meets [xe^{−U}, x].
7. *Contiguity and coverage:* asserted (Λ_a ≤ previous Λ_b; every chain ends with a tail; uniform chains start
   below log 13).
8. *Is a uniform c ≥ 0.4837 really forced?* Yes. ℚ(√−3) has K = ℚ, C = {1}, least split prime 7, Δ = 3.
9. *Possible weak point:* Lemma 3.1 uses standard qualitative facts (zero counting in unit windows, L'/L at
   well-chosen heights). They only make the contour integrals tend to 0 and contribute no constant.
10. *Riemann-sum bound in (4.1):* it needs F unimodal, which holds by (K1): log-concave B-spline times exponential.
    The three bounds per prime are combined with Σ min ≤ min Σ, never the other way round.
11. *Trivial zeros:* Φ(−j) > 0 is used in closed form, and the j > 30 tail is geometric with ratio e^U/x < 1.
    cert.py asserts x > e^U.
12. *Direct range B:* the reported least primes are valid witnesses (split: two primes above 𝔭; inert: one
    unramified prime of relative degree 2). The field counts per group are asserted against the nflist output
    (200 fields: C2 157, C3 3, S3 32, C4 1, V4 3, D4 2, A4 0, S4 2). So are the class counts per field: |G| classes
    over ℚ for abelian L/ℚ, and two per quadratic subfield.

**Reviewer passes (2026-11-06, two `reviewer` subagents).**
- §3.2–§4.3: every piece was judged OK. The one gap: x > e^U was not stated in (3.1) and §4.3. It was always assumed
  (Lemma 3.1) and asserted by cert.py; it is now stated explicitly.
- §3.1, §5.1–5.3 and the cert code: two points, both repaired in the text and neither changing any number.
  - Choose T_k so that ±T_k both avoid the ordinates (now in Lemma 3.1).
  - Upper-rounding x^{−δ}M_x multiplies B, which can be negative at the artificial endpoint n_max. This is safe;
    see the paragraph after (5.2).
  - The tail inequalities and n-endpoint logic were judged correct in direction.

## 8. References
- Bach, Sorenson, Explicit bounds for primes in residue classes, Math. Comp. 65 (1996).
- Cureton, Thorner, Hexagon 2610.00200v2 (31/32 for every ζ_K and Hecke L; Lean: AxiomMath/DedekindZetaZeroFree).
- Frei, Sofos, Widmer, arXiv 2610.18491v2, Thm 1.1 and Cor. 1.2 (constant c(δ) "effectively computable").
- Hu, Hexagon 2610.00216v1, Cor. chebotarev / (9.10) (non-explicit).
- Iwaniec, Kowalski, Analytic Number Theory, Prop. 5.7, Thm 5.11.
- Kadiri, Ng, Wong, arXiv 1902.08640 (≪ Δ_L^{16}); Kadiri, Wong, arXiv 2105.14181 (Δ_L^{310} for all L ≠ ℚ);
  Zaman, arXiv 1508.00287 (≪ Δ_L^{40}); Thorner, Zaman, arXiv 1604.01750; Ahn, Kwon (Δ_L^{12577}).
- Lagarias, Odlyzko, Effective versions of the Chebotarev density theorem (1977).
- Liu, Hexagon 2611.00082 (B_2; unrefereed).
- Rosser, Schoenfeld, Approximate formulas for some functions of prime numbers, Illinois J. Math. 6 (1962), Thm 9.
- Stark, Some effective cases of the Brauer–Siegel theorem, Invent. Math. 23 (1974), Lemma 1.
