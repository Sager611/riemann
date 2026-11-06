# Least prime primitive roots from the zero-free half-plane: derivation

Every step and every numeric constant used by cert.py (section numbers match its comments). Notation: p an odd prime,
n = p − 1, G = (ℤ/pℤ)^×, Υ = rad(n) = l_1 l_2 ⋯ l_ω with primes 2 = l_1 < l_2 < ⋯ < l_ω, ω = ω(p − 1); p_j is the j-th prime,
P_ω = p_1 ⋯ p_ω, V_ω = Π_{j≤ω}(1 − 1/p_j); g*(p) is the least prime that generates G; Λ von Mangoldt; L = log p.

## 0. Inputs

- **(H)** Bebel (QuasiRiemannTracker, Lean-kernel-checked): ζ(s) and every Dirichlet L(s, χ) are zero-free on
  Re s > B_2 + 2.2·10⁻²⁵⁴, B_2 = (777 − √1081)/852. We use **B = 0.8733819666132076 > B_2 + 10⁻¹⁷** (certified),
  b = 1 − B, 1/b = 7.89776916646…, 3/b = 23.6933075…
- **(EP)** Riemann/Deep/ExplProg/DERIVATION.md, used verbatim (cert.py imports ExplProg/cert.py by path): Lemma EF (explicit
  formula for weights supported in [t₀, t₁] ⊂ (0, ∞)), (1.1) (trivial zeros), Lemma P (pairing ρ ↦ 1 − ρ̄), Lemma Z (zero sums
  from the Bennett–Martin–O'Bryant–Rechnitzer count, primitive χ of conductor ≤ e^L, L ≥ log 9), Lemma Zζ (Σ_ζ 1/γ² ≤ 0.0463712),
  Lemma K1 and kernel 3.1: w(t) = g(log t) t^{−B}, g the triangle on [−7, 0] with peak at −21/10, support of w = [e^{−7}, 1],
  W(1) = 2.424245725, W(0) = 115.9708945, ‖w‖_∞ = 38.854381, and for Y = log x − 7 > 0 the envelope constants
  (n₀, k₁(Y), k₂(Y)), (m₀, k₁^L(Y), k₂^L(Y)) of Lemma K1. Z(L), Z_L(L) are the zero-sum bounds of Lemma Z for these envelopes;
  z₁ = 0.0231856·k₂, z₁^L = 0.0231856·k₂^L (Lemma Zζ). The majorant sup_{L≥L_m} Z(L)/L is ExplProg sec. 4 ("tail lemma").
- **(RS62)** Rosser–Schoenfeld, Illinois J. Math. 6 (1962), used only in sec. 5 and only for x ≥ 286 (θ: x ≥ 41):
  R1/R2 (Thm 5) |Σ_{p≤x} 1/p − log log x − M| < 1/(2 log² x); R3/R4 (Thm 7) e^{−γ}(1 − 1/(2log²x))/log x < Π_{p≤x}(1 − 1/p)
  < e^{−γ}(1 + 1/(2log²x))/log x; R5 (Thm 4) θ(x) > x(1 − 1/log x) (x ≥ 41). M cancels in every use.
- Sieve facts proved below (Lemmas S1, S2, T1–T3); the shifted sieve is Iwaniec's (Demonstratio Math. 11 (1978), Lemma 1),
  in the form quoted by G. Martin (Acta Arith. 80 (1997), arXiv math/9807104, §3, the shifted-sieve lemma); we re-derive it.

No GRH, no verification of GRH to a height, no zero density estimate and no unrefereed zero count is used.

## 1. Detection

**Lemma 1.1.** Let x > e⁸ and S(x) := Σ_{m ≥ 1, p ∤ m, m mod p generates G} Λ(m) w(m/x). If S(x) > 0 then g*(p) ≤ x.
*Proof.* Some m = r^k (r prime, r ≠ p) with m ∈ supp w(·/x) ⊂ [x e^{−7}, x] generates G; then r generates G (⟨r^k⟩ ⊆ ⟨r⟩),
and r ≤ r^k ≤ x. ∎ (Prime powers are therefore not removed; no ψ − θ correction is needed.)

For squarefree d | Υ let G^d be the subgroup of d-th powers (index d; G cyclic, so m ∈ G^d iff m ∈ G^l for every prime l | d) and
A_d := Σ_{p∤m, m ∈ G^d} Λ(m) w(m/x). Orthogonality on G/G^d: 1_{G^d}(m) = d^{−1} Σ_{χ mod p, χ^d = 1} χ(m), so

  A_d = A_1/d + r_d,  r_d = d^{−1} Σ_{χ^d = 1, χ ≠ χ₀} S_χ,  |r_d| ≤ (1 − 1/d) R ≤ R,  r_1 = 0,                    (1.1)

S_χ := Σ_m Λ(m) χ(m) w(m/x), R := max_{χ ≠ χ₀ mod p} |S_χ|, and A_1 = Σ_{p∤m} Λ(m) w(m/x).

**Lemma 1.2.** Let Y = log x − 7 > 0, L₉ = max(log p, log 9). Then
  R ≤ (x^B Z(L₉) + x^{1/2} Z_L(L₉))/2 + W(0)/(1 − e^{−2Y}),                                                      (1.2)
  A_1 ≥ x W(1) − x^B z₁ − x^{1/2} z₁^L − W(0)/(1 − e^{−2Y}) − ‖w‖_∞ (7 + log p).                                    (1.3)
*Proof.* p is prime, so every χ ≠ χ₀ mod p is primitive of conductor p ≥ 3. Lemma EF (x t₀ = e^Y > 1) gives
S_χ = −Σ_ρ x^ρ W(ρ) − Σ_{m∈T(χ)} x^{−m}W(−m); the trivial-zero sum is ≤ W(0)/(1 − e^{−2Y}) by (EP 1.1); Lemma P, Lemma K1 and
Lemma Z (conductor p ≤ e^{L₉}) bound the ρ-sum by (x^B Z + x^{1/2} Z_L)/2. For A_1: Lemma EF for ζ gives
Σ_m Λ(m)w(m/x) = xW(1) − Σ_ρ x^ρW(ρ) − (trivial ≥ 0) ≥ xW(1) − x^B z₁ − x^{1/2} z₁^L − W(0)/(1 − e^{−2Y}) (Lemmas P, Zζ, K1, as
EP's Ω₁ bound); the removed terms m = p^k have Λ = log p, w ≤ ‖w‖_∞ and p^k ∈ [x e^{−7}, x] for at most 1 + 7/log p values of k. ∎

## 2. Sieve lemmas

Lower-bound weights: real λ_d on d | P_ω with σ_n := Σ_{d | n} λ_d ≤ [n = 1] for every n | P_ω. Put
T := Σ_{d|P_ω} λ_d/d and C := Σ_{d|P_ω} |λ_d|.

**Lemma S1 (DFS sets).** Let D be a set of divisors of P_ω such that, writing d = q_1 q_2 ⋯ q_r with q_1 > ⋯ > q_r, (i) 1 ∈ D and
q_1⋯q_r ∈ D ⇒ q_1⋯q_{r−1} ∈ D, (ii) every e ∉ D with q_1⋯q_{r−1} ∈ D ("boundary") has r even. Then λ_d = μ(d)1_D(d) are
lower-bound weights.
*Proof.* For n | P_ω every divisor d of n either lies in D or has a unique shortest prefix e ∉ D, and then e is a boundary element
and d = e f with f | n_e := Π{primes of n/e below the least prime of e}. Hence
[n = 1] = Σ_{d|n} μ(d) = Σ_{d|n, d∈D} μ(d) + Σ_{e boundary, e|n} μ(e)[n_e = 1], and μ(e) = +1 by (ii). ∎
cert's weight sets: Rosser's linear-sieve rule (a prime q at an even position m is admitted iff q_1⋯q_{m−1} q³ < y, odd positions
are free) satisfies (i) and (ii) by construction. sanity.py verifies σ_n ≤ [n = 1] by brute force for every n | P_ω, ω ≤ 13,
11 levels y (0 violations).

**Lemma S2 (shifted sieve).** Let λ be lower-bound weights on P_ω, Υ = l_1⋯l_ω as above and d ↦ d' the multiplicative bijection
p_j ↦ l_j. For any finite weighted set with w ≥ 0, a map a ↦ ν(a) | Υ, A_{d'} := Σ_{d' | ν(a)} w(a), A_1 = X and
|A_{d'} − X/d'| ≤ R (d ≠ 1):
  S := Σ_{ν(a) = 1} w(a) ≥ Σ_d λ_d A_{d'} ≥ X (V_Υ/V_ω) T − R(C − |λ_1|),  V_Υ := Π_j (1 − 1/l_j).                        (2.1)
*Proof.* Pointwise, with m_a the preimage of ν(a): Σ_{d | m_a} λ_d = σ_{m_a} ≤ [m_a = 1] = [ν(a) = 1]; multiply by w(a) and sum:
S ≥ Σ_d λ_d A_{d'} = X Σ_d λ_d/d' + Σ_d λ_d r_{d'}. Let g(p_j) = 1/l_j and h(p) = g(p)/(1 − g(p)). For multiplicative g,
Σ_{d|P_ω} λ_d g(d) = Σ_{m|P_ω} σ_m Π_{p|m} g(p) Π_{p|P_ω/m}(1 − g(p)) = V_Υ Σ_m σ_m Π_{p|m} h(p) (expand g(d) as the probability
that a random subset of primes, each chosen independently with probability g(p), contains d). As l_j ≥ p_j,
0 < h(p_j) ≤ 1/(p_j − 1), and σ_m ≤ 0 for m ≥ 2, so Σ_m σ_m Π_{p|m} h(p) ≥ Σ_m σ_m Π_{p|m} 1/(p − 1) = T/V_ω (the same identity
with g(p) = 1/p). ∎
**Master inequality.** In our setting a = m (p ∤ m), w(a) = Λ(m)w(m/x), ν(m) = Π{l | Υ : m ∈ G^l}, so A_{d'} is the A_{d'} of
sec. 1 and S(Υ) = S(x); λ_1 = 1 for all weights used. If T ≥ 0 and A_1 ≥ 0, then since V_Υ ≥ V_ω,
  S(x) ≥ A_1 T − R (C − 1).                                                                                            (2.2)
Hence g*(p) ≤ x as soon as A_1^{lo} T > R^{up}(C − 1) with (1.2), (1.3). sanity.py checks the main-term inequality of (2.1)
(Σ λ_d/d' ≥ V_Υ T/V_ω) exactly in rationals for 200 random prime sets {l_j} per (ω, y), 40 pairs, ω ≤ 10.

## 3. Weights for ω ≤ W₁ = 10⁴ (sieve.py, grid.py)

For a level y the set D(y) (Lemma S1, Rosser's rule) restricted to divisors of P_ω gives T_ω(y), C_ω(y). sieve.py enumerates D(y)
depth first; leaves at odd positions under an even node with prod·q·8 ≥ y (no admissible grandchild since q' ≥ 2) are summed in
bulk from prefix sums of ⌊2⁶⁴/q⌋ and ⌈2⁶⁴/q⌉. T is accumulated as an exact integer lower bound in units 2⁻⁶⁴ (floor for μ = +1
terms, ceiling subtracted for μ = −1), C exactly; cumulative in the index of q_1 = P⁺(d), which gives every ω ≤ W₁ at once.
grid.py tries y = ⌊10^{k/8}⌉, k = 4, …, 82, and keeps for each ω the y minimising (C − 1)/T (any choice is valid; 62 distinct levels,
largest 4.22·10⁹). Examples: ω = 10: y = 750, T = 0.0942, C = 34; ω = 100: y = 133352, T = 0.02948, C = 1704; ω = 10⁴:
y = 4.22·10⁹, T = 0.007882, C = 11887948, (C − 1)/T = 0.1375 p_ω². Optimal levels are y ≈ p_ω^{1.9–2.0}; at finite ω, T/V_ω is
still positive below s = log y/log p_ω = 2 (the limit f(s) = 0 at s = 2 is not reached at these sizes).

## 4. Theorems for ω(p − 1) ≤ W₁ (cert.py omega_section)

Fix ω ≤ W₁ and the level y_ω of sec. 3; T = T_ω > 0 (certified), C = C_ω. Let p be a prime with ω(p − 1) = ω; then
L = log p ≥ L_ω := log(P_ω + 1) (p − 1 ≥ Π l_j ≥ Π p_j). Take x = (κL)^{1/b}, Y = log x − 7. Dividing (1.3) and (1.2) by x^B L:
  A_1/(x^B L) ≥ W(1)κ − (z₁ + x^{1/2−B} z₁^L)/L − W(0)/((1 − e^{−2Y}) x^B L) − ‖w‖_∞(7/L + 1)/x^B =: Xo(L),          (4.1)
  R/(x^B L) ≤ Z(L₉)/(2L) + x^{1/2−B} Z_L(L₉)/(2L) + W(0)/((1 − e^{−2Y}) x^B L) =: Ro(L).                          (4.2)
**Monotonicity (4.3).** Admissibility: κ > 0 with Y_m := log x(L_m) − 7 > 1 (cert rejects κ otherwise); Lemma 1.1 holds for every
x > 0, Lemma 1.2 needs Y > 0. For fixed κ, x and Y increase with L, so k₁, k₂, k₁^L, k₂^L (decreasing in Y) may be frozen at Y_m;
then z₁, z₁^L, x^{1/2−B}, x^{−B} and 1/L are non-increasing on [L_m, ∞) and Xo(L) ≥ Xo(L_m). For Ro we do not claim monotonicity:
with M_m := sup_{L≥L_m} Z(L₉)/L and M_m^L := sup_{L≥L_m} Z_L(L₉)/L (frozen k's; the ExplProg majorant at max(L_m, log 9), and for
L_m < log 9 the maximum of that and Z(log 9)/L_m, Z_L(log 9)/L_m), every L ≥ L_m has
Ro(L) ≤ R̂_m := M_m/2 + x_m^{1/2−B} M_m^L/2 + W(0)/((1 − e^{−2Y_m}) x_m^B L_m), x_m = x(L_m). Hence
  T·Xo(L_m) − (C − 1)·R̂_m > 0                                                                                         (4.4)
implies A_1 T > R(C − 1), i.e. S(x) > 0 by (2.2), for every L ≥ L_m (cert's sieve_margin computes exactly Xo(L_m) and R̂_m). cert finds the least κ = m·10^e (3 significant digits,
rounded up; fixed-point estimate, then the Arb sign of (4.4)).

**Theorem A.** For every ω ≤ 10⁴ and every prime p with ω(p − 1) = ω: g*(p) ≤ (κ_ω log p)^{1/(1−B)}, κ_ω in logs/kappa_table.txt
(L_m = L_ω). Examples: κ₁ = 3.53 (p − 1 a power of 2), κ₂ = 11.1, κ₅ = 54.1, κ₁₀ = 216, κ₂₀ = 868, κ₅₀ = 6130,
κ₁₀₀ = 2.72·10⁴, κ₁₀₀₀ = 4.32·10⁶, κ₁₀⁴ = 6.85·10⁸.
**Theorem B.** With c_ω := (κ_ω/L_ω²)^{1/3} (rounded up, Arb-checked c_ω³L_ω² ≥ κ_ω): κ_ω L ≤ c_ω³L³ for L ≥ L_ω, so
g*(p) ≤ (c_ω log p)^{3/(1−B)}. max_{ω≤10⁴} c_ω = 1.4312 (ω = 2); c_ω ≤ 0.474 for ω > 100, 0.398 at ω = 10⁴.
With L_m = max(L_ω, log 10⁷) (primes p > 10⁷) the maximum is 0.8139 → c = 407/500 = 0.814 (ω = 8); primes p ≤ 10⁷ are covered by sec. 6.

## 5. The tail ω(p − 1) > W₁ (tail.py; cert.py tail_section)

Here z := p_ω ≥ Z₁ := p_{10001} = 104743, u := log z ≥ u₁ = log 104743 = 11.55926…, L = log p ≥ L₁ := log(P_{10001} + 1) = 104403.76, and all primes
≤ z are (shifted) sieving primes. Parameters: s = 23/10, D := z^s, w₀ = 50, Q₂ = 10⁵, δ = 1/10. P(t) := Π_{p<t} p.

**Weights.** For each prime q ∈ [w₀, z], ρ^{(q)} are Selberg's optimal Λ² weights for the primes < q at level ξ_q := (D/q)^{1/2}:
ρ_d = μ(d)(d/φ(d)) G_{q,d}(ξ_q/d)/G_q(ξ_q) (d | P(q), d < ξ_q), G_q(ξ) := Σ_{d|P(q), d<ξ} μ²(d)/φ(d), G_{q,d} the same sum over
(m, d) = 1; ρ_1 = 1. Put λ_e := μ(e)[e | P(w₀)] − Σ_{w₀≤q≤z} Σ_{d₁,d₂ | P(q), q[d₁,d₂] = e} ρ^{(q)}_{d₁} ρ^{(q)}_{d₂}.

**Lemma T1.** λ are lower-bound weights on P_ω with λ_1 = 1.
*Proof.* σ_n = [(n, P(w₀)) = 1] − Σ_{q | n, q ≥ w₀} (Σ_{d | (n, P(q))} ρ^{(q)}_d)² (q ∤ [d₁, d₂]). If n has a prime factor < w₀,
σ_n ≤ 0; if n > 1 has none, its least prime q₀ ≥ w₀ has (n, P(q₀)) = 1, its square is ρ₁² = 1 and σ_n ≤ 0; σ_1 = 1. ∎

**Lemma T2 (Selberg).** Σ_{d₁,d₂} ρ_{d₁}ρ_{d₂}/[d₁, d₂] = 1/G_q(ξ_q) and Σ_d |ρ_d| ≤ c_σ ξ_q/G_q(ξ_q), c_σ := Π_p(1 + 2/(p(p − 1))) ≤ 3.2795975.
*Proof.* The first is Selberg's minimum for density 1/p (h(p) = 1/(p − 1), h(d) = 1/φ(d)). Next |ρ_d| = (d/φ(d))G_{q,d}/G_q, and
Σ_d μ²(d)(d/φ(d)) Σ_{m<ξ/d, (m,d)=1, m|P(q)} μ²(m)/φ(m) = Σ_{k<ξ, k|P(q)} μ²(k) Σ_{d|k} d/φ(k) = Σ_k μ²(k)σ(k)/φ(k).
As μ²(k)σ(k)/φ(k) = μ²(k)Σ_{e|k} μ²(e)Π_{p|e} 2/(p − 1), Σ_{k<ξ} ≤ ξ Σ_e μ²(e)Π_{p|e} 2/(p(p − 1)) = c_σ ξ. cert: exact product over
p < 3·10⁵ times exp(2/(N − 1)) (Σ_{n≥N} 2/(n(n − 1)) = 2/(N − 1)). ∎

**Main term and cost.** With V(q) := Π_{p<q}(1 − 1/p), V⁺(z) := Π_{p≤z}(1 − 1/p) and Σ_{w₀≤q≤z} V(q)/q = V(w₀) − V⁺(z) (telescoping):
  T = V(w₀) − Σ_q 1/(q G_q(ξ_q)) = V⁺(z) − Σ_{w₀≤q≤z} (V(q)/q)(Φ_q − 1),  Φ_q := 1/(V(q) G_q(ξ_q)) ≥ 1,                (5.1)
  C − 1 ≤ 2^{π(w₀)} + Σ_q (Σ_d |ρ^{(q)}_d|)² ≤ 2^{15} + c_σ² D Σ_{w₀≤q≤z} 1/(q G_q(ξ_q)²).                                (5.2)

**Lemma T3.** For ξ > 1: (a) G_q(ξ) ≥ log min(ξ, q); (b) G_q(ξ) ≥ log ξ − Σ_{q≤l<ξ, l prime}(1 + log(ξ/l))/l;
(c) V(q) G_q(ξ) ≥ 1 − ξ^{−α} Π_{p<q}(1 + (p^α − 1)/p) for every α > 0.
*Proof.* μ²(d)/φ(d) = Σ_{rad m = d} 1/m, so G_q(ξ) ≥ H_q(ξ) := Σ_{m<ξ, P⁺(m)<q} 1/m. (a) G increases with ξ; at ξ' = min(ξ, q)
every m < ξ' has P⁺(m) < q and Σ_{m<ξ'} 1/m ≥ log ξ'. (b) Σ_{m<ξ} 1/m ≥ log ξ; an m < ξ with P⁺(m) ≥ q is l·k (l ≥ q prime,
k < ξ/l), counted at least once in Σ_l (1/l)Σ_{k<ξ/l} 1/k, and Σ_{k<y} 1/k ≤ 1 + log y (y ≥ 1). (c) Rankin:
1/V(q) − G_q(ξ) = Σ_{d|P(q), d≥ξ} μ²(d)/φ(d) ≤ ξ^{−α} Π_{p<q}(1 + p^α/(p − 1)) = ξ^{−α} V(q)^{−1} Π_{p<q}(1 + (p^α − 1)/p). ∎
*(b) made explicit* (used only when ξ > q, and the bound (5.3) only when ξ > x₀; for ξ ≤ q use (a), for q < ξ ≤ 293 the exact sum,
as there are no primes in [288, 293)). Let x₀ = max(q, 288) ≥ 286, S(x) = Σ_{p≤x} 1/p. By R1/R2, S(x) − S(x₀⁻) < log(log x/log x₀) + η₀ for x ≥ x₀,
η₀ = 1/log²x₀ + 1/x₀. Stieltjes integration by parts with f(x) = 1 + log(ξ/x) (decreasing, f(ξ) = 1):
  Σ_{x₀≤l<ξ} f(l)/l = f(ξ)(S(ξ⁻) − S(x₀⁻)) + ∫_{x₀}^{ξ}(S(x) − S(x₀⁻)) dx/x ≤ log a₀ + η₀(1 + log(ξ/x₀)) + log x₀(a₀ log a₀ − a₀ + 1),   (5.3)
a₀ = log ξ/log x₀; primes l ∈ [q, 288) are summed exactly. For q > Q₂, dividing by log q:
  G_q(ξ) ≥ log q · h_Q(a),  h_Q(a) = 2a − a log a − 1 − η_B(a − 1) − (log a + η_B)/log Q₂,  a = log ξ/log q ≥ 1,  η_B = 1/log²Q₂ + 1/Q₂,  (5.4)
and h_Q may be evaluated at min(a, a_cap), a_cap = 2.609 (below the maximum of h_Q; G increases with ξ).

**Part A (50 ≤ q ≤ Q₂, numeric).** For each prime q, V(q) exactly; G_lo(q, ξ) := max of (a), (b) [(5.3) with the exact sum below
288] and (c) for α ∈ {0.03, 0.05, 0.08, 0.12, 0.18, 0.25, 0.35, 0.5} with the exact product; term (V(q)/q)(Φ_q − 1) ≤ (1/q)(1/G_lo − V(q)).
**Part B (q > Q₂, uniform).** Φ_q − 1 ≤ ψ(a_q), a_q := log ξ_q/log q = (su/log q − 1)/2, with ψ non-increasing:
  a ≤ 1: ψ(a) = e^γ/((1 − ε_M) a) − 1   [(a) and R3: 1/(V(q) log q) < e^γ/(1 − 1/(2log²q)); ε_M = 1/(2 log²Q₂)];
  a ≥ 1: ψ(a) = min{ e^γ/(1 − ε_M) − 1,  e^γ/((1 − ε_M) h_Q(min(a, a_cap))) − 1,  min_{c: r_c<1} r_c/(1 − r_c) },  r_c = exp(−ca + R_B(c)) (candidates with r_c ≥ 1 are dropped; certified r_c < 1),
  c ∈ {1/2, 1, 3/2, 2, 5/2, 3, 4, 5, 6}, R_B(c) := Σ_{p≤283}(p^{c/log Q₂} − 1)/p + Σ_{k≥1} c^k/(k·k!) + (e^c − 1)(1/(2log²Q₂) + 1/(2log²288)).
R_B(c) ≥ Σ_{p<q}(p^{c/log q} − 1)/p for q > Q₂: primes below 288 by monotonicity in the exponent; for 288 ≤ p < q, Abel summation of
f(x) = x^{c/log q} − 1 against S = log log x + M + E: the main part is ∫_{log288/log q}^{1}(e^{cv} − 1)dv/v ≤ Σ c^k/(k·k!), the
E-part is ≤ f(q)/(2log²q) + f(288)/(2log²288) + (f(q) − f(288))/(2log²288). Then (c) with α = c/log q (ξ^{−α} = e^{−ca}) and
Π(1 + y) ≤ e^{Σy}. With V(q) ≤ e^{−γ}(1 + ε_V)/log q, 1 + ε_V = (1 + 1/(2log²Q₂))/(1 − 1/Q₂) (R4):
  Σ_{Q₂<q≤z} (V(q)/q)(Φ_q − 1) ≤ e^{−γ}(1 + ε_V) Σ_{Q₂<q≤z} F(q)/q,  F(x) = ψ(a(x))/log x.
Abel summation against S on (Q₂, z] (R1/R2), t := log x/u, Ψ(t) := ψ((s/t − 1)/2) (non-decreasing):
  Σ F(q)/q ≤ (1/u) ∫_{t₂}^{1} Ψ(t) dt/t² + (2Ψ(1) + 4J + K₃)/(2u³),  t₂ = log Q₂/u,                                     (5.5)
J := ∫_0^1 Ψ/t⁴, K₃ := sup Ψ(t)/t³: the boundary terms are F(z)E(z) ≤ Ψ(1)/(2u³) and F(Q₂)E(Q₂) ≤ K₃/(2u³); |dF| ≤ dΨ/(tu) + Ψ dt/(t²u),
|E| < 1/(2t²u²) and ∫dΨ/t³ ≤ Ψ(1) + 3J give ∫|E||dF| ≤ (Ψ(1) + 4J)/(2u³). Integrals: upper Riemann sums on 4000 cells of [1/50, 1];
on (0, 1/50], Ψ ≤ r_c/(1 − r_c(1/50)), r_c(t) = e^{c/2+R_B(c)} e^{−cs/(2t)}, with ∫_0^{t_m} e^{−k/t}t^{−2}dt = e^{−k/t_m}/k,
∫_0^{t_m} e^{−k/t}t^{−4}dt = e^{−kW}(W²/k + 2W/k² + 2/k³) (W = 1/t_m), sup_{w≥W} w³e^{−kw} = W³e^{−kW} (W ≥ 3/k).
**τ.** V⁺(z) > e^{−γ}(1 − 1/(2u²))/u (R3), so (5.1) gives T ≥ τ V⁺(z) for every z ≥ Z₁ with τ := 1 − sup_u (Part A + Part B)/V⁺.
205 cells u ∈ [u₁, 60] (width 0.1, 0.25, 0.5): Part A at the left end (G increases with ξ), 1/V⁺ at the right end, Part B with
the integrals over [t_i, 1], t_i ≤ log Q₂/u_right (grid point below), for I, J and K₃ alike (the Abel-summation error in (5.5)
only involves t ∈ [t₂, 1]: ∫_{t₂}^1 dΨ/t³ ≤ Ψ(1) + 3∫_{t₂}^1 Ψ/t⁴, and F(Q₂)E(Q₂) ≤ sup_{t≥t₂} Ψ/t³/(2u³)), error terms at the left end. For u ≥ 60 (normalised bounds): Part A/V⁺ ≤ Σ_q (V(q)/q) r_q/(1 − r_q) · e^γ u e^{−κ(u−60)}/(1 − 1/(2u²)) (one α = 0.12, r_q at
u = 60, κ = αs/2; the last factor decreases on u ≥ 60 because 60κ ≥ 1), so it is ≤ its value at u = 60; Part B/V⁺ ≤
(1 + ε_V)(I + (2Ψ(1) + 4J + K₃)/(2u²))/(1 − 1/(2u²)) with the full integrals, decreasing in u, so ≤ its value at u = 60. Results: min over the 205 cells 0.111196 (at u ≈ 59.56); for u ≥ 60 Part A ≤ 0.009356, Part B ≤ 0.889211; **τ = 0.101433** (lower endpoint used).
(1 + ε_V) I = 0.8889 < 1 is also checked; I = 0.882775, J = 3.27487, K₃ = 4.72882, Ψ(1) = 1.75049 (upper bounds).
**Cost.** Σ_{w₀≤q≤z} 1/(q G_q(ξ_q)²) ≤ K_A + K_B, K_A := Σ_{50≤q≤Q₂} 1/(q G_lo(q, ξ_q(u₁))²) (G increases with u),
K_B := g_B^{−2}(1/(2log²Q₂))(1 + 2/log²Q₂), g_B = min(1, (s − 1)/2) [(a) with a_q ≥ (s − 1)/2; Σ_{q>Q₂} 1/(q log²q) by Abel
summation]. Values: K_A ≤ 0.0160555, K_B ≤ 0.00906309.

**Final inequality.** For p with ω(p − 1) = ω > W₁: z = p_ω, L ≥ L₁ and θ(z) < L, so z < L/(1 − 1/u₁) (R5),
u ≤ ℓ(L) := log(L/(1 − 1/u₁)). Take x = (c_T L)^{β/b}, β = 1 + s + δ = 17/5. By (2.2), (1.2), (1.3), (5.1), (5.2), divided by
x^B L^{1+s}, g*(p) ≤ x as soon as
  τ'(W(1) c_T^β L^δ/ℓ(L) − sub(L)/(L^s ℓ(L))) > Ro(L)(2^{15}/L^s + c_σ²(K_A + K_B)/(1 − 1/u₁)^s),  τ' = τ e^{−γ}(1 − 1/(2u₁²)),  (5.6)
with sub = W(1)κ − Xo of (4.1) (the subtracted terms) and Ro of (4.2). Extension from L₁ to all L ≥ L₁, as in (4.3)–(4.4):
freeze k₁, k₂, k₁^L, k₂^L at Y₁ := Y(L₁) (they decrease in Y, and Y increases with L). Then (i) L^δ/ℓ(L) ≥ m_δ on L ≥ L₁ (increasing
there, δ ℓ(L₁) ≥ 1); (ii) sub(L) is bounded by its frozen form, each term of which (z₁/L, x^{1/2−B} z₁^L/L, W(0)/((1 − e^{−2Y}) x^B L),
‖w‖_∞(7/L + 1)/x^B) is non-increasing, so sub(L)/(L^s ℓ(L)) ≤ sub(L₁)/(L₁^s ℓ(L₁)); (iii) for Ro no monotonicity is claimed: with
M₁ := sup_{L≥L₁} Z(L)/L and M₁^L := sup_{L≥L₁} Z_L(L)/L (frozen k's; the ExplProg majorant at L₁ ≥ log 9),
  Ro(L) ≤ R̂₁ := M₁/2 + x₁^{1/2−B} M₁^L/2 + W(0)/((1 − e^{−2Y₁}) x₁^B L₁)  for every L ≥ L₁,  x₁ = x(L₁),
since x^{1/2−B} and the trivial-zero term are non-increasing in L; (iv) L^{−s} is decreasing. So (5.6) with its left side at L₁
(m_δ, sub(L₁)) and Ro replaced by R̂₁ proves (5.6) for every L ≥ L₁; cert's tail_section margin() evaluates exactly this
(ZoverL(·, L₁) = M₁, M₁^L; k's from KER.ks(Y₁)).
m_δ = 0.272694 (δ ℓ(L₁) = 1.16 ≥ 1); the least c_T on a 10⁻³ grid with a certified positive margin of (5.6) is 391/200.
**Theorem C.** Every prime p with ω(p − 1) > 10⁴ has g*(p) ≤ (1.955 log p)^{(1+s+δ)/(1−B)} ≤ (1.955 log p)^{26.853}.

## 6. Small primes (cert.py small_section, small.py)

Exact g*(p) for all 664578 primes 3 ≤ p ≤ 10⁷: smallest-prime-factor table (numpy), the distinct primes l | p − 1, and the first
prime r with r^{(p−1)/l} ≢ 1 (mod p) for every l (Python integers, exact). Maximum 211 (p = 4022911); records 3:2, 7:3, 23:5, 41:7,
109:11, 191:19, 271:43, 2791:53, 11971:79, 31771:107, 190321:149, 2080597:151, 3545281:163, 4022911:211. Proof of the bound for 5 ≤ p ≤ 10⁷: max g*(p) = 211 < (0.814 log 5)^{3/(1−B)} = 601.44 (Arb), and (0.814 log p)^{3/(1−B)} increases with p.
(A float scan of every p, worst log g*/log bound = 0.108 at p = 5, is printed as information only.) p = 3 (g* = 2) is excluded from the 0.814 form
((0.814 log 3)^{3/(1−B)} < 1) and covered by Theorem B with 1.4312.

## 7. Final statements (cert.py uniform_c, corollary_thresholds)

**Theorem D (every prime).** For every prime p ≥ 3: g*(p) ≤ (1.956 log p)^{17/(5(1−B))} ≤ (1.956 log p)^{26.853}
(17/(5(1−B)) = 26.85241…; the base exceeds 1, so rounding the exponent up is safe).
*Proof.* ω(p − 1) > 10⁴: Theorem C (1.955 < 1.956). ω(p − 1) = ω ≤ 10⁴: Theorem A and κ_ω L ≤ c^{17/5}L^{17/5} for L ≥ L_ω whenever
c ≥ (κ_ω/L_ω^{12/5})^{5/17}; cert checks c^{17/5}L_ω^{12/5} ≥ κ_ω in Arb for all 10⁴ values (the binding case is the tail). ∎
**Theorem B'.** For every prime p ≥ 5 with ω(p − 1) ≤ 10⁴ (in particular every p ≥ 5 with log p < 104403.76):
g*(p) ≤ (0.814 log p)^{3/(1−B)} ≤ (0.814 log p)^{23.694} (3/(1−B) = 23.69330…). (Sec. 4 for p > 10⁷, sec. 6 for p ≤ 10⁷.)
**Corollary E (cross-checks).** (i) g*(p) < p for every prime p with log p ≥ 94.4783 (sufficient: p ≥ 10^{41.032}); (ii) g*(p) < √p − 2
for every prime p with log p ≥ 224.610 (sufficient: p ≥ 10^{97.547}). *Proof.* For each ω ≤ 120 whose table κ_ω does not already give (1/b) log(κ_ω L) < L
(resp. < L/2 − 10⁻³) at L_ω, cert raises L_m in steps of 1/4 from L_ω, recertifies κ at L_m (sec. 4) and stops at the first L_m
with the inequality; the difference L − (1/b) log(κL) (resp. L/2 − …) increases for L > 1/b (resp. 2/b), checked. Binding:
ω = 27 (i), ω = 52 (ii). For ω > 120 the table values already satisfy it at L_ω (checked for every ω ≤ 10⁴), and for ω > 10⁴
Theorem C gives log x = (17/(5(1−B))) log(1.955 L) < L/2 − 10⁻³ at L ≥ 10⁵, increasing in L. For (ii), log(x + 2) ≤ log x + 2/x and x > e⁸. ∎
(i) matches the threshold 10⁴¹ of Zenodo 23623384 (QRH-conditional, unrefereed, plain sieve with x = p); (ii) is the first
statement of Grosswald type for g*(p) not conditional on GRH (GRH: p > 2791, McGown–Treviño–Trudgian, arXiv 1508.05182).

**Table** (log₁₀ of the bound; "worst ω" is max over the ω allowed at that size of (κ_ω log p)^{1/(1−B)}):

| p ≈ | max ω(p−1) | worst ω (Thm A) | (0.814 log p)^{3/(1−B)} | (1.956 log p)^{3.4/(1−B)} | ω = 2 | ω = 5 |
|---|---|---|---|---|---|---|
| 10¹⁰ | 10 | 10^29.2 | 10^30.2 | 10^44.4 | 10^19.0 | 10^24.4 |
| 10²⁰ | 16 | 10^34.8 | 10^37.3 | 10^52.5 | 10^21.4 | 10^26.8 |
| 10¹⁰⁰ | 53 | 10^49.0 | 10^53.9 | 10^71.3 | 10^26.9 | 10^32.3 |
| 10¹⁰⁰⁰ | 350 | 10^71.0 | 10^77.5 | 10^98.1 | 10^34.8 | 10^40.2 |
| 10¹⁰⁰⁰⁰ | 2584 | 10^94.1 | 10^101.2 | 10^125.0 | 10^42.7 | 10^48.1 |

## 8. Certificate conventions

python-flint Arb balls (200 bits in cert.py/ExplProg, 120 bits in tail.py); every inequality used as a proof step is a certified sign
(pos(), else the run aborts). Floats only choose parameters (κ, c, cells), and every chosen value is rounded up and then certified.
Sieve main terms: exact integers (sec. 3). Inputs are used as stated (RS62: 1/2, 286, 41) or rounded in the safe direction
(B up; c_σ, K_A, K_B, I, J, K₃ as upper endpoints; τ as the lower endpoint). One command:
`cd Riemann/Deep/PrimRootNum && uv run --no-project --with python-flint --with mpmath --with numpy python3 cert.py all`
(≈ 4 min, one thread; output logs/cert_output.txt, logs/kappa_table.txt, logs/tail_output.txt). logs/sieve_best.json only records the
level y_ω chosen per ω (grid.py, 219 s); cert.py recomputes T_ω, C_ω from it, and any choice of y would be valid.

## 9. Adversarial checklist (self-review)

- **p − 1 highly composite.** The shifted sieve's worst case is p − 1 = P_ω·k (first ω primes); Theorems A/B take the worst case
  for every ω, and the tail covers ω > 10⁴. adv.py: exact g*(p) for the least prime ≡ 1 (mod P_m) (m ≤ 45, p up to 10⁸¹),
  Sophie Germain p = 2l + 1, Fermat primes, 2^a3^b + 1: all ≤ (κ_ω log p)^{1/b}, worst log g*/log bound 0.083.
- **Smallest p.** p = 3, 5, 7: Theorem A holds with L_m = log(P_ω + 1) ≥ log 3 and Lemma Z at log 9 (conductor 3, 5, 7 < 9);
  p = 3 is excluded only from the 0.814 form. p = 2 is excluded (no odd prime; g*(2) = 1 is not prime).
- **Sophie Germain / ω = 2, Fermat / ω = 1.** Table rows κ₂ = 11.1, κ₁ = 3.53 (the sieve is then the plain Vinogradov sieve).
- **Imprimitive characters.** None: p is prime, every χ ≠ χ₀ mod p is primitive of conductor p.
- **Prime powers m = r^k.** Kept in the sum; r^k generating G forces r to generate G (Lemma 1.1). Powers of p are removed in (1.3).
- **Real χ (quadratic character).** Covered by (H); no exceptional zero is assumed or excluded separately.
- **No GRH, no zero density, no unrefereed count.** Only (H), BMOR (refereed, Lemma Z), Davenport's Σ Re 1/ρ identity (Lemma Zζ), RS62.
- **Downward roundings.** None: κ, c rounded up; τ lower endpoint; all bounds upper endpoints. Every max/min over Arb balls
  (τ over cells, K₃, thresholds Λ, G_lo candidates) compares exact endpoints, and ψ's branches are decided on exact rationals.
- **Monotonicity claims** (sec. 4.3, 5.6, Part A in u, Ψ in t, h_Q on [1, a_cap]) are argued in the text; Ro is never claimed
  monotone: (4.4) and (5.6) use its supremum majorants R̂_m, R̂₁ over L ≥ L_m, L₁ with frozen kernel constants; a_cap = 2.609 lies below the
  maximiser of h_Q (h_Q′(2.609) > 0 since 1 − log 2.609 − η_B − 1/(2.609 log Q₂) > 0).
- **RS62 ranges.** adv.py checks R1–R5 numerically at every prime x ≤ 10⁷ (max |E| log²x = 0.464 < 1/2, max |V e^γ log x − 1| log²x
  = 0.453 < 1/2, min θ(x)/(x(1 − 1/log x)) = 1.066 > 1); the proof uses the published statements.
- **Lemma S1 and the main-term inequality of Lemma S2** are checked by brute force / exact rationals in sanity.py (0 violations);
  **sieve.py's bulk-leaf enumeration** is checked against a plain DFS listing every element (exact rational T, count C) for
  (ω, y) = (20, 10⁴), (20, 10⁶), (50, 10⁴), (50, 10⁵), (100, 10⁵), (100, 3·10⁵), (200, 10⁵): identical C, T within the 2⁻⁶⁴ rounding.

## 10. Constants (where they live)

| constant | value | file (rg) |
|---|---|---|
| B, b, 1/b | 0.8733819666132076, 0.1266180…, 7.8977691664… | ExplProg/cert.py (Fr(8733819666132076, 10 ** 16)) |
| kernel ℓ, rise | 7, 49/10 (peak at −21/10) | ExplProg/cert.py (Tri(7, Fr(49, 10))) |
| W(1), W(0), ‖w‖_∞ | 2.424245725, 115.9708945, 38.854381 | computed in ExplProg/cert.py (KER.W1, KER.W0, KER.winf) |
| Σ_ζ 1/γ²/2 | 0.0231856 | ExplProg/cert.py (SUM_INV_G2) |
| sieve unit | 2⁻⁶⁴ (K = 64) | sieve.py (K = 64) |
| levels y | ⌊10^{k/8}⌉, k = 4…82 | grid.py (10 ** (k / 8)) |
| W₁ | 10⁴ | logs/sieve_best.json (Wmax), grid.py run |
| small-p range | 10⁷ | cert.py (10 ** 7) |
| tail s, δ, w₀, Q₂, u_A | 23/10, 1/10, 50, 10⁵, 60 | cert.py (Fr(23, 10), Fr(1, 10), w0=50, Q2=100000, uA=60) |
| α grid (Part A) | 0.03 … 0.5 | tail.py (Fr(3, 100) …) |
| c grid (Part B) | 1/2 … 6 | tail.py (Fr(1, 2) …) |
| split points | 288, 293, 283 | tail.py |
| a_cap | 2609/1000 | tail.py (computed: argmax − 0.001) |
| quadrature | 4000 cells on [1/50, 1] | tail.py (n=4000, tmin=Fr(1, 50)) |
| decay α for u ≥ 60 | 12/100 | cert.py (Fr(12, 100)) |
| c_σ | ≤ 3.2795975 (primes < 3·10⁵) | tail.py (300000) |
| RS62 | 1/2, 286, 41 | tail.py (docstring), used via 288, Q2 |
| results | κ_ω table; 1.4312; 407/500; 391/200; 489/250; 94.4783; 224.610 (upper ends) | logs/cert_output.txt, logs/kappa_table.txt |
