# Explicit Littlewood bound for L(1, χ) from the zero-free half-plane — proof

Notation. χ is a primitive Dirichlet character mod q ≥ 3, a ∈ {0, 1} with χ(−1) = (−1)^a, ρ = β + iγ runs over the zeros
of L(s, χ) with multiplicity, γ_E = 0.5772… is Euler's constant, ψ = Γ'/Γ, M = 0.2614972128… is the Meissel–Mertens
constant, E = −1.3325822757… = −γ_E − Σ_{n≥2} Σ_p (log p)/p^n.

Hypothesis H(B): L(s, χ) ≠ 0 for Re s > B, for every Dirichlet character χ. It is a theorem for
B = 0.87338196661320752 (Lean-checked, Bebel, 31 Oct 2026; B = 7/8 by OpenAI's quasi-Riemann-hypothesis theorem, Lean
port in `Riemann/Deep/QRH`). Throughout B is this number and b := 1 − B = 0.12661803338679248.

## Results

**Theorem 1.** Assume H(B). For every primitive character χ mod q ≥ 3 (real or complex):

1. |L(1, χ)| ≥ 0.11 / max(3, log log q);
2. |L(1, χ)| ≥ 0.115 / log log q if log log q ≥ 3.13, and |L(1, χ)| ≥ 0.116 / log log q if log log q ≥ 4;
3. |L(1, χ)| ≥ c_N / log log q if log log q ≥ Λ_N, for (Λ_N, c_N) = (40, 0.115659), (100, 0.116424), (1000, 0.116873),
   (10⁴, 0.116933); and liminf_{q→∞} min_χ |L(1, χ)| log log q ≥ (π²/6) e^{−γ_E} (1 − B) = 0.1169398….

**Corollary 2.** For every fundamental discriminant −d < −4,
h(−d) ≥ 0.0350 √d / max(3, log log d), h(−d) ≥ 0.0366 √d / log log d for log log d ≥ 3.13, and
h(−d) > 100, 1000, 2000, 10⁴, 10⁵ as soon as d ≥ 6.973·10⁷, 7.292·10⁹, 2.957·10¹⁰, 7.627·10¹¹, 7.976·10¹³.

**Theorem 3 (side result, upper bound).** Assume H(B). For every primitive χ mod q with log log q ≥ 3,
|L(1, χ)| ≤ 15 log log q; the constant tends to e^{γ_E}/(1 − B) = 14.0665 (14.2223 for log log q ≥ 40, 14.1288 for ≥ 100,
14.0746 for ≥ 1000). This beats the classical |L(1, χ)| ≤ ½ log q + O(1) from about log q = 110 on
(certified 59.35 < 74.2 = ½ log q at log q = e⁵).

With B = 7/8 instead (OpenAI's theorem) the certificate (`cert_output_7_8.txt`) gives 0.105/max(3, log log q),
0.11/log log q for log log q ≥ 3, 0.114/log log q for log log q ≥ 4 (limit (π²/6)e^{−γ_E}/8 = 0.11545), and
h(−d) > 2000, 10⁵ for d ≥ 3.047·10¹⁰, 8.137·10¹³.

Sanity check of the shape: with B = 1/2 (GRH) the same argument gives the constant (π²/6)e^{−γ_E}/2 = π²/(12 e^{γ_E}),
which is Littlewood's 1928 constant under GRH. With the half-plane the price is the factor 2(1 − B) = 0.2532.

## 1. The weight and its transform

Fix L ≥ 8 and 0 < h ≤ 2, put x = e^L and

v(u) = 1 (u ≤ 1), (1 + h − u)/h (1 ≤ u ≤ 1 + h), 0 (u ≥ 1 + h);  V(w) = (e^{(1+h)Lw} − e^{Lw}) / (h L w²).

V is entire except for a simple pole at w = 0 with residue 1 (the Taylor expansion gives V(w) = 1/w + (1 + h/2)L + O(w)).

**Lemma 1.** For n ≥ 1 and c > 0: (1/2πi) ∫_{(c)} V(w) n^{−w} dw = v(log n / L).

Proof. (1/2πi) ∫_{(c)} e^{wy} w^{−2} dw = max(y, 0) for real y. Apply it with y = (1 + h)L − log n and y = L − log n and
divide by hL. ∎

**Lemma 2 (explicit formula).** For every real σ ≥ 1,

−L'/L(σ, χ) = Σ_n Λ(n) χ(n) v(log n / L) n^{−σ} + Σ_ρ V(ρ − σ),

the sum over all zeros of L(s, χ) (trivial ones included), absolutely convergent.

Proof. By Lemma 1 the left-hand sum equals (1/2πi) ∫_{(c)} −L'/L(σ + w, χ) V(w) dw for c > 0 (absolute convergence:
|V(w)| ≪ |w|^{−2} on Re w = c). Move the line to Re w = −U, U → ∞ avoiding trivial zeros, along horizontal segments at
heights ±T_j where |L'/L(s, χ)| ≪ log²(qT_j) for −1 ≤ Re s ≤ 2, and use |L'/L(s, χ)| ≪ log(q|s|) for Re s ≤ −1 at
distance ≥ 1/2 from trivial zeros (Davenport, Multiplicative Number Theory, Ch. 16–17). Since |V(w)| ≪ e^{L Re w}|w|^{−2}
for Re w ≤ 0, the new integrals tend to 0. The residues are −L'/L(σ, χ) at w = 0 and −V(ρ − σ) at w = ρ − σ. ∎

These bounds for L'/L are only used qualitatively, to justify the contour shift; the identity itself is exact.

**Lemma 3.** Define log L(1, χ) := −∫_1^∞ L'/L(σ, χ) dσ (H(B) gives L(σ, χ) ≠ 0 for σ ≥ 1; Re log L(1, χ) = log|L(1, χ)|
for every branch). Then

log L(1, χ) = P_0 + Σ_ρ G(ρ),  P_0 := Σ_n Λ(n) χ(n) v(log n / L) / (n log n),  G(ρ) := ∫_1^∞ V(ρ − σ) dσ.

Proof. Integrate Lemma 2 over σ ∈ [1, ∞). The interchange is justified by Lemma 5 below (Σ_ρ ∫ |V(ρ − σ)| dσ < ∞) and
∫_1^∞ n^{−σ} dσ = 1/(n log n). ∎

## 2. The zeros

**Lemma 4 (Hadamard).** Z := Σ_{ρ nontrivial} Re 1/(1 − ρ) = Re L'/L(1, χ) + ½ log(q/π) + ½ ψ((1 + a)/2).

Proof. Davenport Ch. 12, (17) and (18): L'/L(s, χ) = −½ log(q/π) − ½ ψ((s + a)/2) + B(χ) + Σ_ρ (1/(s − ρ) + 1/ρ) with
Re B(χ) = −Σ_ρ Re(1/ρ). Take real parts at s = 1. ∎

Under H(B) every nontrivial zero has 0 < β ≤ B, so each term of Z is (1 − β)/|1 − ρ|² > 0. The functional equation
makes ρ ↦ 1 − ρ̄ an involution of the multiset of nontrivial zeros (Λ(s, χ) = ε Λ(1 − s, χ̄) and the zeros of L(s, χ̄) are
the conjugates of those of L(s, χ)), so also β ≥ 1 − B.

**Lemma 5 (one zero).** Let σ ≥ 1 and ρ = β + iγ nontrivial. Then |ρ − σ| ≥ |1 − ρ| and

|V(ρ − σ)| ≤ (e^{−(1+h)L(σ−β)} + e^{−L(σ−β)}) / (hL|1 − ρ|²).

Hence, for the weight (L, h), |G(ρ)| ≤ g_0(β)/|1 − ρ|² with g_0(β) = (e^{−(1+h)L(1−β)}/(1 + h) + e^{−L(1−β)})/(hL²),
and for a second weight (L_1, h_1) with transform V_1, |V_1(ρ − 1)| ≤ g_1(β)/|1 − ρ|² with
g_1(β) = (e^{−(1+h_1)L_1(1−β)} + e^{−L_1(1−β)})/(h_1 L_1). Both g_0 and g_1 are increasing in β.

Proof. |e^{(1+h)Lw} − e^{Lw}| ≤ e^{(1+h)L Re w} + e^{L Re w}; |ρ − σ|² = (σ − β)² + γ² increases in σ ≥ β; integrate
e^{−cL(σ−β)} over σ ≥ 1. ∎

**Lemma 6 (pairs).** If g ≥ 0 is increasing on (0, 1) and |F(ρ)| ≤ g(β)/|1 − ρ|² for every nontrivial zero, then
Σ_ρ |F(ρ)| ≤ K Z with K := g(B)·B/b + g(1/2).

Proof. Since ρ ↦ ρ' := 1 − ρ̄ permutes the zeros, Σ_ρ |F(ρ)| = ½ Σ_ρ (|F(ρ)| + |F(ρ')|) and Z = ½ Σ_ρ (w(ρ) + w(ρ')),
w(ρ) = (1 − β)/|1 − ρ|². It suffices to show |F(ρ)| + |F(ρ')| ≤ K(w(ρ) + w(ρ')) for each ρ, and by symmetry we may take
β ≥ 1/2 (ρ' = 1 − β + iγ). Write α = 1 − β ≤ c = β, u = γ²; then |1 − ρ|² = α² + u, |1 − ρ'|² = c² + u,
w(ρ) = α/(α² + u), w(ρ') = c/(c² + u), |F(ρ)| ≤ g(β)/(α² + u), |F(ρ')| ≤ g(1 − β)/(c² + u).
Multiplying by α² + u and putting r = (α² + u)/(c² + u) ∈ [α²/c², 1), the claim reads g(β) + g(1 − β) r ≤ K(α + c r), linear
in r, so it suffices to check r = 1 and r = α²/c². At r = α²/c² it is K ≥ K(β) := g(β) c/α + g(1 − β) α/c; at r = 1 it is
K ≥ g(β) + g(1 − β), which follows from K(β) ≥ g(β) + g(1 − β) ⇔ (c − α)(g(β)/α − g(1 − β)/c) ≥ 0, true because
g(β) ≥ g(1 − β) and α ≤ c. Finally g(β)β/(1 − β) ≤ g(B)B/b (increasing in β ≤ B) and g(1 − β)(1 − β)/β ≤ g(1/2). ∎

The pairing improves the single-zero bound g(B)/b by the factor B + b·g(1/2)/g(B), essentially B ≈ 0.873. For each β
the worst height is γ = 0 (r = α²/c²), so the low zeros are the extremal case; the real pair (B, 1 − B) gives the ratio
g(B)B/b + g(b)b/B, slightly below K.

**Lemma 7 (trivial zeros).** The trivial zeros are −m, m ≥ 0, m ≡ a (mod 2), simple. For them V(−m − σ) < 0 and
|V(−m − σ)| ≤ e^{−L(m+σ)}/(m + σ). Hence −Σ_triv V_1(−m − 1) ∈ [0, T_1] and Σ_triv G(−m) ∈ [−T_0, 0] with
T_1 = e^{−L_1}/(1 − e^{−2L_1}) and T_0 = e^{−L}/(L(1 − e^{−2L})).

Proof. e^{−(1+h)Ly} − e^{−Ly} ∈ [−e^{−Ly}·hLy, 0) for y > 0; sum the geometric series. ∎

## 3. The bound for Z and for log|L(1, χ)|

Put κ_0 := g_0(B)B/b + g_0(1/2), κ_1 := g_1(B)B/b + g_1(1/2) (Lemma 6) and assume κ_1 < 1. Let P_1 be the sum of
Lemma 2 at σ = 1 for the weight (L_1, h_1): P_1 = Σ_n Λ(n)χ(n) v_1(log n/L_1)/n.

From Lemma 2 at σ = 1 (weight (L_1, h_1)), Lemmas 5–7:
Re L'/L(1, χ) = −Re P_1 − Σ_ρ Re V_1(ρ − 1) − Σ_triv V_1 ≤ −Re P_1 + κ_1 Z + T_1. With Lemma 4 and ψ((1 + a)/2) ≤ ψ(1) = −γ_E:

(1 − κ_1) Z ≤ −Re P_1 + ½ log(q/π) − γ_E/2 + T_1.   (Z-bound)

From Lemma 3 (weight (L, h)), Lemmas 5–7: log|L(1, χ)| ≥ Re P_0 − κ_0 Z − T_0. With κ' := κ_0/(1 − κ_1):

log|L(1, χ)| ≥ Re P_0 + κ' Re P_1 − κ'(½ log(q/π) − γ_E/2 + T_1) − T_0.   (★)

## 4. The prime sums

Write z = χ(p), |z| ≤ 1, x = e^L, x_1 = e^{L_1}, and let τ(y) := Σ_{p<y} Σ_{k≥2, p^k>y} p^{−k}.

**Lemma 8.** (i) Re P_0 ≥ −Σ_p v(log p/L)/p + Σ_{p<x} (1/p − log(1 + 1/p)) − τ(x)/2 − 1/(2(x − 1)).
(ii) Re P_1 ≥ −Σ_p v_1(log p/L_1)(log p)/p − L_1 τ(x_1) − (1 + h_1)L_1/(2(x_1 − 1)).

Proof. (i) For p < x, v(k log p/L) = 1 when p^k < x, so Σ_k z^k v_k/(k p^k) = −log(1 − z/p) − Σ_{p^k>x} z^k(1 − v_k)/(k p^k)
with real part ≥ −log(1 + 1/p) − ½ Σ_{k≥2, p^k>x} p^{−k}. For p ≥ x the k = 1 term has real part ≥ −v(log p/L)/p and the
k ≥ 2 terms are at most ½ Σ_{k≥2} p^{−k} = 1/(2p(p − 1)) in absolute value; Σ_{p≥x} 1/(p(p − 1)) ≤ 1/(x − 1).
(ii) For p < x_1: Re Σ_k z^k (log p)/p^k = (log p) Re(w/(1 − w)) with |w| ≤ 1/p, which is ≥ −(log p)/(p + 1) ≥ −(log p)/p;
the corrections have p^k > x_1 and log p < L_1. For p ≥ x_1 the k ≥ 2 terms vanish unless k log p < (1 + h_1)L_1, so they
are at most ((1 + h_1)L_1/2) Σ_{k≥2} p^{−k}. ∎

**Lemma 9 (explicit Mertens).** With ε_D(L) = 1/(10L²) + 4/(15L³), E_h := (1/h)∫_1^{1+h} log t dt = ((1 + h)log(1 + h) − h)/h:

(a) Σ_p v(log p/L)/p ≤ log L + M + E_h + ε_D(L) for L ≥ 10;
(b) Σ_p v_1(log p/L_1)(log p)/p ≤ L_1(1 + h_1/2) + E + 0.2/L_1 + 0.2/L_1² for L_1 ≥ 8;
(c) Σ_{p<x} (1/p − log(1 + 1/p)) ≥ M − γ_E + log(π²/6) − 1/(2(x − 1));
(d) τ(y) ≤ τ̄(y) := y^{−1/2} + ((log y/log 2)(y^{1/3} + 3)/2 + 8)/y. In ℓ = log y, τ̄ = e^{−ℓ/2} + (ℓ/(2 log 2))(e^{−2ℓ/3} + 3e^{−ℓ})
+ 8e^{−ℓ}; it is decreasing for ℓ ≥ 3/2, and ℓτ̄ is decreasing for ℓ ≥ 3 (each of ℓe^{−ℓ/2}, ℓ²e^{−2ℓ/3}, ℓ²e^{−ℓ}, ℓe^{−ℓ}
is). τ itself is not monotone: it jumps up by Σ_{k≥2} p^{−k} = 1/(p(p − 1)) when y passes a prime p. Only τ̄ is used below.

Proof. (a) v(log p/L) = (1/h) meas{t ∈ [1, 1 + h] : p < e^{Lt}}, so Σ_p v(log p/L)/p = (1/h)∫_1^{1+h} S(e^{Lt}−) dt with
S(y) = Σ_{p≤y} 1/p ≤ log log y + M + ε_D(log y) for y ≥ 10372 (Dusart, arXiv:1002.0442, Thm 6.11); ε_D decreases.
(b) The same with R(y) = Σ_{p≤y}(log p)/p ≤ log y + E + 0.2/log y + 0.2/log² y for y ≥ 2974 (Dusart, Thm 6.12);
e^{L_1} ≥ e^8 > 2974. (c) 1/p − log(1 + 1/p) = [1/p + log(1 − 1/p)] − log(1 − 1/p²); summed over all p this is
(M − γ_E) + log ζ(2); the tail over p ≥ x is ≤ Σ_{p≥x} 1/(2p²) ≤ 1/(2(x − 1)).
(d) For 2 ≤ k ≤ K := ⌊log y/log 2⌋ let M_k = ⌊y^{1/k}⌋ + 1; then Σ_{m≥M_k} m^{−k} ≤ M_k^{−k}(1 + M_k/(k − 1))
≤ (1/y)(1 + (y^{1/k} + 1)/(k − 1)); for k > K, Σ_{m≥2} m^{−k} ≤ 3·2^{−k} and Σ_{k>K} 3·2^{−k} ≤ 6/y. Sum. ∎

## 5. Master inequality

**Proposition 10.** Assume H(B). For every primitive χ mod q ≥ 3, every L ≥ 10, L_1 ≥ 8, 0 < h, h_1 ≤ 2 with κ_1 < 1:

log|L(1, χ)| ≥ F := C_0 − log L − E_h − [ε_D(L) + 1/(x − 1) + τ̄(x)/2 + T_0] − κ' W,

C_0 = log(π²/6) − γ_E = −0.0795153…, W = L_1(1 + h_1/2) + E + 0.2/L_1 + 0.2/L_1² + L_1 τ̄(x_1) + (1 + h_1)L_1/(2(x_1 − 1))
+ ½ log(q/π) − γ_E/2 + T_1.

Proof. Insert Lemmas 8 and 9 into (★), with τ ≤ τ̄ (Lemma 9(d)); C_0 = (M − γ_E + log(π²/6)) − M. ∎

F depends on q only through −κ'·½ log q, so for fixed parameters F is decreasing in log q. The certificate `cert.py`
evaluates exactly this F (with τ̄) in ball arithmetic (python-flint arb, 256 bits) at exact rational parameters, and
returns no value (cell fails) unless the ball for κ_1 lies below 1.

**Proposition 11 (upper bound).** Under the same hypotheses,

log|L(1, χ)| ≤ F_up := γ_E + log L + E_h + ε_D(L) + 1/(2(x − 1)) + κ' W_up,
W_up = L_1(1 + h_1/2) − γ_E + 0.2/L_1 + 0.2/L_1² + ½ log(q/π) − γ_E/2 + T_1.

Proof. Lemma 3 with Σ_triv G ≤ 0 and Lemma 6 give log|L(1, χ)| ≤ Re P_0 + κ_0 Z ≤ Re P_0 − κ' Re P_1 + κ'(½ log(q/π) − γ_E/2 + T_1).
For p < x, Re Σ_k z^k v_k/(k p^k) ≤ −log(1 − 1/p); for p ≥ x it is ≤ v(log p/L)/p + 1/(2p(p − 1)). Since
Σ_p (−log(1 − 1/p) − 1/p) = γ_E − M has positive terms, Re P_0 ≤ Σ_p v/p + γ_E − M + 1/(2(x − 1)), and Lemma 9(a) applies.
Next −Re P_1 ≤ Σ_p v_1 (log p)/p + Σ_p Σ_{k≥2} (log p)/p^k = Σ_p v_1 (log p)/p − γ_E − E, and Lemma 9(b) applies. ∎

The proof of Theorem 3 is the same as that of Theorem 1: 7400 cells on [3, 40] (worst margin 0.0188), and for Λ ≥ Λ_N the
tail Lemma T' (upper form): under the hypotheses of Lemma T (Λ_N/b ≥ 10, 0 < h ≤ 2, κ_1(Λ_N) < 1), with L = L_1 = Λ/b,
h_1 = h, |L(1, χ)| ≤ (e^{γ_E}/b) e^{E_h + Err_up(Λ_N)} Λ, Err_up = ε_D(L) + 1/(2(x − 1)) + κ'(L(1 + h/2) + e^Λ/2). Proof as for
Lemma T: κ_1 < 1 on the whole tail by Lemma T(i); W_up ≤ L(1 + h/2) + ½ e^Λ because −γ_E + 0.2/L + 0.2/L² − ½ log π − γ_E/2
+ T_1 is non-increasing in L and ≤ −1.4097 < 0 at L = 8 (certified, section (0) of cert_output.txt); Err_up decreases by
Lemma T(iii). Certified: 14.2223 for Λ_N = 40, h = 0.015 (κ_1(Λ_N) ≤ 9.6·10⁻¹⁸).

## 6. Proof of Theorem 1

Let Λ = log log q.

(1) Λ ≤ 3: with P = (12.854, 1.016, 8.000, 1.832), F(e³) = −3.287948 > log(0.11/3); since F decreases in log q, this covers
every q ≥ 3 with Λ ≤ 3, and |L(1, χ)| ≥ 0.03733.
3 ≤ Λ ≤ 40: 7400 cells [Λ_i, Λ_i + 0.005]; on each cell F(exp(Λ_{i+1}); P_i) > log(0.11/Λ_i) ≥ log(0.11/Λ) (worst margin
0.0174). Λ ≥ 40: Lemma T below with Λ_N = 40, h = 0.015 (κ_1(Λ_N) ≤ 9.6·10⁻¹⁸ < 1) gives 0.115659 > 0.11.
(2) Same with 0.115 on [3.13, 40] (7374 cells, worst margin 0.00047) plus Lemma T (40); with 0.116 on [4, 100] (9600 cells
of 0.01, worst margin 0.0228) plus Lemma T (100, h = 0.006; κ_1(Λ_N) ≤ 8.4·10⁻⁴⁴): 0.116424.
(3) Lemma T.

**Lemma T (tail).** Fix 0 < h ≤ 2 and Λ_N with Λ_N/b ≥ 10, and assume κ_1 < 1 at L = L_1 = Λ_N/b, h_1 = h (write
κ_1(Λ_N)). For Λ = log log q ≥ Λ_N apply Proposition 10 with L = L_1 = Λ/b, h_1 = h, x = e^L. Then κ_1 < 1 for these
parameters, so Proposition 10 applies, and

log|L(1, χ)| ≥ C_0 − log(Λ/b) − E_h − Err(Λ_N),  Err(Λ) := ε_D(L) + 1/(x − 1) + τ̄(x)/2 + T_0 + κ'(L(1 + h/2) + e^Λ/2).

Proof. (i) κ_1 on the whole tail. For fixed β < 1, g_1(β) = (e^{−(1+h)L(1−β)} + e^{−L(1−β)})/(hL) is a sum of terms
e^{−cL}/(hL) with c > 0, hence decreasing in L; so κ_1 = g_1(B)B/b + g_1(1/2) is decreasing in L = Λ/b, and
κ_1(Λ) ≤ κ_1(Λ_N) < 1 for every Λ ≥ Λ_N.
(ii) The term W. With L_1 = L, W = L(1 + h/2) + ½ e^Λ + ω(L), ω(L) := E + 0.2/L + 0.2/L² + Lτ̄(e^L) + (1 + h)L/(2(e^L − 1))
− ½ log π − γ_E/2 + T_1, T_1 = 1/(2 sinh L). Every L-dependent term of ω is non-increasing for L ≥ 8 (Lemma 9(d) for Lτ̄;
L/(e^L − 1) and 1/sinh L decrease), and ω increases with h, so ω(L) ≤ ω(8)|_{h=2} ≤ −1.7237 < 0 (certified, section (0) of
cert_output.txt, with E replaced by the upper bound −1.33258227573). Hence W ≤ L(1 + h/2) + ½ e^Λ and κ'W ≤ κ'(L(1 + h/2) + e^Λ/2).
(iii) Monotonicity of Err in Λ. ε_D(L), 1/(x − 1) and T_0 = 1/(2L sinh L) decrease in L; τ̄(x) decreases (Lemma 9(d); this
is why τ̄ and not τ appears). Each term of κ_0 = g_0(B)B/b + g_0(1/2) has the form e^{−cL}/(hL²), c > 0, so κ_0 and κ_0 L
decrease. With e^{−bL} = e^{−Λ}:
κ_0 e^Λ = (B/b)(e^{−hΛ}/(1 + h) + 1)/(hL²) + (e^{−((1+h)/2 − b)L}/(1 + h) + e^{−(1/2 − b)L})/(hL²),
which decreases because b < 1/2. By (i), 1/(1 − κ_1) decreases too, so κ' = κ_0/(1 − κ_1), κ'L and κ'e^Λ decrease.
Therefore Err(Λ) ≤ Err(Λ_N), and Proposition 10 with (ii) gives the claim (log L = log(Λ/b)). ∎

Certified (cert_output.txt, section (C)): κ_1(Λ_N) ≤ 9.6·10⁻¹⁸, 8.4·10⁻⁴⁴, 6.1·10⁻⁴³⁵, 1.4·10⁻⁴³⁴³ for (Λ_N, h) = (40, 0.015),
(100, 0.006), (1000, 0.001), (10⁴, 10⁻⁴), giving the constants c_N = e^{C_0} b e^{−E_h − Err(Λ_N)} of Theorem 1(3).
A numerical spot check (section (0), not part of the proof) confirms that κ_1 and Err decrease on Λ = 40·{1, …, 10}.

Hence |L(1, χ)| ≥ e^{C_0} b e^{−E_h − Err(Λ_N)} / Λ for Λ ≥ Λ_N.

Limit. Take h = 1/Λ_N. Then hL = 1/b at L = Λ_N/b, so κ_1(Λ_N) ≤ 2B e^{−Λ_N} + 2b e^{−Λ_N/(2b)} → 0 (the hypothesis of
Lemma T holds for large Λ_N); κ_0 e^{Λ_N} ≤ 2(B/b + 1)/(hL²) = 2b/Λ_N → 0; κ_0 L ≤ 2e^{−Λ_N} → 0; E_h ≤ h/2 → 0
(since (1 + h)log(1 + h) − h ≤ h²/2); and ε_D, 1/(x − 1), τ̄, T_0 → 0. So E_h + Err(Λ_N) → 0 and the constants tend to
e^{C_0} b = (π²/6)e^{−γ_E}(1 − B) = 0.1169398…. ∎

## 7. Proof of Corollary 2

For a fundamental discriminant −d < −4 the character χ_{−d} is primitive, real and odd, of conductor d, and
h(−d) = √d L(1, χ_{−d})/π (w = 2). Theorem 1 gives the first two bounds (0.11/π = 0.03501…, 0.115/π = 0.03661…).
Thresholds: on cells [ℓ_i, ℓ_i + 0.01] of ℓ = log d from the stated threshold up to ℓ_end (20.09, 22.80, 24.23, 27.53,
32.22) the certificate checks e^{ℓ_i/2} e^{F(ℓ_{i+1}; P_i)}/π > H; beyond ℓ_end, √d·0.11/(π log log d) = e^{ℓ/2}·0.11/(π log ℓ) is
increasing (ℓ log ℓ > 2 for ℓ ≥ 20) and exceeds H at ℓ_end. ∎

## 8. Remarks

* Complex χ cost nothing: the per-prime bounds use only |χ(p)| ≤ 1 and are attained at χ(p) = −1, which real characters
  can realise. So Theorem 1 is uniform over all primitive χ.
* The branch of log L(1, χ) never matters: only log|L(1, χ)| = Re log L(1, χ) is used.
* Region |t| small: there is no separate treatment; in Lemma 6 the worst height for every β is γ = 0, and
  Z = Σ(1 − β)/|1 − ρ|² weighs every height correctly.
* What limits the constant: x must be about (log q)^{1/(1−B)}, so the Euler product runs to log x ≈ log log q/(1 − B); the
  worst case χ(p) = −1 for all p ≤ x then gives (π²/6)e^{−γ_E}/log x. A better B improves the constant linearly in 1 − B.
* Inputs: H(B); Dusart's Theorems 6.11, 6.12 (arXiv:1002.0442, Ramanujan J. 45 (2018)); Davenport Ch. 12, 16, 17
  (qualitative). No zero counts, no convexity bounds, no Borel–Carathéodory.
