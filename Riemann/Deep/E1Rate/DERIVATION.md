# E1Rate — derivation (e1rate, 14 Nov 2026)
Sources: SL = Santibáñez-Leal, Zenodo 23134787 v0.03, TeX main.tex from the CAOS zip (Zenodo 23135354), line numbers "l." as in E1Audit/FULLPROOF.md;
E1 = Riemann/Deep/SantiAudit/E1.md; Baluyot = thesis (Rochester 2017) ch. 3, mapped and audited in Riemann/Deep/Baluyot/MAP.md (Conrey's 4/7 argument);
DI = Deshouillers–Iwaniec, Acta Arith. 43 (1984) 305–312, Lemma 1 (from Invent. Math. 70 (1982) Thm 12); BC = Bettin–Chandee, Adv. Math. (2018) Thm 1.
Exponents: all sizes are powers of T; "relative" = divided by the main term ≍ H = T^θ; ε-losses suppressed. Code: exponents.py (log), kappa_opt.py (log).

## 1. The off-diagonal after the dual step, as a trilinear form
SL (14) (l. 469–475) with E1 Step 1 (x-space form): for h = gq, k = gp, (p,q) = 1, gp, gq ≤ y = T^ν, the oscillatory dual term is
  O = Σ_g Σ_{p,q} b_{gq} b_{gp} (4π/(gpq)) (q/p)^{−z/2} Σ_{n≥1} n^{−z/2} σ_{−α,β}(n) e(−p̄n/q) 𝒯(n),   𝒯(n) = −π∫A(ξ)ξ^{−z/2}B(4π√(nξ/(pq)))dξ.
Key structural fact: 𝒯 depends on (n, p, q) only through x = n/(pq): 𝒯(n) = 𝒲(n/(pq)), 𝒲(x) := −π∫A(ξ)ξ^{−z/2}B(4π√(ξx))dξ.
Reciprocity (p̄ = p⁻¹ mod q, q̄ = q⁻¹ mod p): p̄/q + q̄/p ≡ 1/(pq) (mod 1), so e(−p̄n/q) = e(n q̄/p)·e(−x). Hence
  O = Σ_g (4π/g) Σ_{(p,q)=1} Σ_n β_p γ_q α_n e(n q̄/p) · Φ(n/(pq)),   Φ(x) = e(−x)𝒲(x),
  β_p = b_{gp}p^{z/2−1}, γ_q = b_{gq}q^{−z/2−1}, α_n = n^{−z/2}σ_{−α,β}(n) (|α_n| ≤ d(n)n^{O(1/L)}),
  b_h = μ(h)h^{σ0−1/2}P(log(y/h)/log y) (SL Thm 1.2, l. 158–176; σ0 = 1/2 − R/L), so b_{gq} = μ(g)μ(q)1_{(q,g)=1}·(smooth in q).
Moduli and lengths (dyadic p ≍ P, q ≍ Q, n ≍ N; P, Q ≤ y/g): modulus p (or q, by the same reciprocity), inverted variable q (or p);
  N0 = PQT/H² (E1 Step 4: 𝒲 decays like (N0/N)^{K/2} beyond), N ≤ N0T^ε; X0 = 2√(2π)√(NT/(PQ)) ≤ 2√(2π)·T/H (the Bessel argument).
Weight sizes (E1 Steps 3–4, audited): |𝒲(x)| ≪ (T/H²)^j H (Tx)^{−1/4}, and 𝒲(x) = Σ_± e^{±iX0}·(non-oscillating amplitude of size H X0^{−1/2}):
  the phase ±4π√(ξx) varies by ≍ X0 over each dyadic block in every one of n, p, q (∂_{log n} = −∂_{log p} = −∂_{log q} = X0/2) — not separable.
At θ = 1: X0 ≪ T^ε (no oscillation), N0 = PQ/T: Conrey's/BCR's setting. At θ < 1: X0 up to T^{1−θ}, N0 up to y²T^{1−2θ}.

## 2. Separating the weight costs exactly √X0 (= SL's multiplier norm)
Mellin in the single variable x: Φ(x) = (2πi)⁻¹∫_{(σ)}Φ̃(s)x^{−s}ds, x^{−s} = n^{−s}p^{s}q^{s}: the twist splits into three unimodular sequences.
Partition dyadically in x (x ≍ x_N = N/(PQ)). Stationary phase in log x (phase t log x ± X0(x), second derivative ≍ X0): Φ̃_N(it) lives on
|t| ≍ X0(N) with |Φ̃_N| ≍ H X0^{−1/2}·X0^{−1/2}, so ∫|Φ̃_N(it)|dt ≍ H = sup|Φ_N|·√X0(N): separation costs Λ(N) = √X0(N) (falsified numerically, §6).
SL's own version (15) (l. 485, one line Re u = σ for all N): multiplier norm ≪ H T^{−σ}(T/H)^{2σ}, per block H(N0/N)^σ; at σ = 1/4
  separated scale = H(N0/N)^{1/4} = (true size H(PQ/(NT))^{1/4}) × √X0 (exactly the factor (T/H)^{1/2} at N = N0; E1Audit §2 end).
So every bound that needs the twist separated (BC; DI's bilinear lemma) is charged Λ = √X0 ≍ (NT/(PQ))^{1/4}, and then
  (trivial, E1)   R_triv(N) = N^{3/4}(PQ)^{1/4}T^{−1/4}  (relative, summed over the block; = PQT^{1/2}H^{−3/2} at N = N0),
  (separated)     R_sep(N) = R_triv·Λ = N  (exactly), and a bound with relative gain G over the trivial ℓ¹ bound gives R = N·G.
Consistency (BC, arbitrary coefficients): G = (NPQ)^{−1/2}·BC-bracket ⇒ R = N^{3/5}(PQ)^{−3/10}… ; worst P = Q = y, N = N0 gives
  T^{(17/20)(1−2θ)+(33/20)ν} + T^{1−2θ+(15/8)ν}: exactly SL (newerror) l. 534–538, i.e. SL's ν_B = min{1/2, (17/33)(2θ−1)} (l. 126) is reproduced.

## 3. Conrey's structure: Vaughan's identity in the inverted variable, DI Lemma 1 (type II), Weil (type I)
Fix g, a block (P, Q, N), and choose the modulus (p, say; q by symmetry). γ_q = μ(q)χ0(q)·F(q) (χ0 principal mod g; F = q^{−z/2−1+σ0−1/2}·P(log(y/(gq))/log y),
a polynomial in log q times a power: binomial expansion separates it exactly into O(1) products). Vaughan with any U = U(Q), 1 ≤ U ≤ Q^{1/2}:
  μχ0 = c1 + c2 + c3, c1 = 2(μχ0)1_{≤U} (q ≤ U); c2 = −(μχ0)_{≤U}*(μχ0)_{≤U}*χ0: q = a·k, a ≤ U² (|coef| ≤ d(a)), k smooth;
  c3 = (μχ0)_{>U}*λ, λ = −((μχ0)_{≤U}*χ0)1_{>U}: q = a·b, a, b > U, |λ(b)| ≤ d(b).   (χ0(k) = Σ_{e|(k,g)}μ(e): smooth after O(d(g)) terms.)
Type II (c3; c2 with a > U1): DI Lemma 1 as in Baluyot L6 (p. 103; audited, Riemann/Deep/Baluyot): for |C(a,ℓ)| ≤ 1, a ≍ A, ℓ ≍ N, b ≍ B, v ≍ V,
  Σ_b Σ_v |Σ_a Σ_ℓ C(a,ℓ) e(ℓ\overline{ab}/v)| ≪ (ABNV)^{1/2+ε}{(VB)^{1/2} + (N+A)^{1/4}[VB(N+A)(V+A²) + NA²B²]^{1/4}}.
  Use: v = p (|β_p p^{s}| ≤ P^{−1}, outside), b outside (λ(b)b^{s}… ≤ T^ε), C(a,n) = coef(a)a^{s}F-part·α_n n^{−s}, after the separation of §2.
  Gain G = {…}/(ABNV)^{1/2}; relative error R_II = N·G = (N/A)^{1/2} + (N/(ABV))^{1/2}(N+A)^{1/2}(VB)^{1/4}(V+A²)^{1/4} + N^{3/4}(N+A)^{1/4}V^{−1/2}
  (both orders (A,B) allowed: b and a are interchangeable bilinear variables).
Type I (c2 with a ≤ U1): no separation. Lemma I (incomplete Kloosterman sum with Bessel phase): for fixed (a, p, n), with Φ_N = Φ on x ≍ N/(PQ),
  |Σ_{k≍D} χ0(k) e(n\overline{ak}/p) Φ_N(n/(pak))| ≪ T^ε sup|Φ_N|·(D(n,p)/p + p^{1/2}(1+X0)^{1/2}),  D = Q/a.
  Proof: Poisson mod p; Weil |S(m,r;p)| ≤ d(p)(m,r,p)^{1/2}p^{1/2}; Φ_N = superposition (E1 Step 4 form, ξ-integral) of amplitudes times e^{±iX(ξ,t)},
  X(ξ,t) = 4π√(nξ/(pat)), ∂_t²X = (3/4)X/t² ≍ X0/D²: the r-th Fourier integral is ≪ D(1+X0)^{−1/2}·amp (second-derivative test) and is
  negligible for |r| ≫ T^ε p(1+X0)/D (integration by parts); r = 0 gives the Ramanujan sum. Sum: (1/p)[D(n,p) + p^{1/2}·p(1+X0)D^{−1}·D(1+X0)^{−1/2}].
  Relative: R_I = R_triv(N)·[1/P + P^{1/2}(1+X0)^{1/2}a/Q] (the (n,p) average is ≪ d(p)). At N = N0: R_I ≪ Q T^{1/2}H^{−3/2} + P^{3/2}a·T/H².
  Meaning: the Bessel phase raises the conductor of the incomplete Kloosterman sum from p to p·X0 ≍ pT/H.
Main terms. P = Q = y, N = N0 = y²T^{1−2θ} (exponent n0 = 2ν + 1 − 2θ < ν/2 in the relevant range):
  type II, balanced A = B = y^{1/2} (DI's spectral term T2; N0 < A, V = A²): N0^{1/2}y^{−1/8} = T^{(1−2θ)/2 + 7ν/8}  ← binding, independent of U;
  type II diagonal T1 = (N0/A)^{1/2}: with the larger factor inside, A ≥ y^{1/2}: ≤ T^{(1−2θ)/2 + 3ν/4} (needs only ν < (2/3)(2θ−1));
  type II unbalanced, smaller factor inside (V > A²): T2 = T^{n0 − ν/4 − α/4} ≤ T^{1−2θ+13ν/8} at α = ν/2 (θ = 1: Baluyot's y ≤ T^{8/13});
  type I, a ≤ U1:        y^{3/2}U1·T^{1−2θ} (U1 = y^{1/4}: T^{1−2θ+7ν/4} = (T2)²; U1 → 1: needs only ν < (2/3)(2θ−1)).
  All < 1 ⟺ ν < (4/7)(2θ − 1) =: ν_C(θ). At θ = 1: 4/7 = Conrey (Baluyot's binding terms T^{−1/2}y^{7/8}, T^{−1}y^{7/4}, MAP.md "Where θ < 4/7 comes from").
Every other monomial, every unbalanced block, every N ≤ N0, both orders and moduli, every U on a grid: exponents.py (grid minimax). Result (exponents.log):
  ν_K(θ) := sup{ν: every block < 1 with min(E1 bound, DI/Weil bound)} = max{(3θ−1)/4, (4/7)(2θ−1)} at θ = 1, .99, .9375, .9, .85, .8182, .8, .75 (and .6, .51);
  worst block always P = Q = y, N = N0 (T2/T1/type I). So DI/Weil beats E1 exactly for θ > 9/11 ((4/7)(2θ−1) = (3θ−1)/4 ⟺ θ = 9/11).


## 4. The rest of SL at the new ν (SL's text; E1 §3 table) and the cap ν < 1/2
| SL step | requirement | status at ν < ν_K(θ) |
|---|---|---|
| J/Y dual branch (14) | this file §2–§3 (DI/Weil) or E1 (absolute) per block | ✓ ν < max{(3θ−1)/4, (4/7)(2θ−1)} |
| K-Bessel branch (l. 540–542) | O(H(y²/T)^σ), any fixed σ | ✓ ν < 1/2 |
| ε = +1 primal profile (l. 543–545) | (H/T)^A | ✓ any ν |
| Estermann residues = main term, Lemma diagonal (l. 285–316) | contour to Re s = −δ: HT^{−δ(1−2ν)} = o(H/L) | ✓ ν < 1/2 (Remark l. 336–338: only use of ν < 1/2) |
| Taylor remainder (11) (l. 388–405) | J > (1+ν−θ+δ)/(2θ−1) | ✓ fixed J |
| compact windows (l. 730–758) | Gaussian width T^{θ−η}; ν < ν_K(θ−η) | ✓ ν_K continuous, η → 0 |
| Vaughan pieces, g-sum | P, Q ≤ y/g; exponents of g negative (as SL l. 535) | ✓ |
So E1-M (SL Thm 1.2) holds for 0 < ν < ν_K(θ) := min{1/2, max{(3θ−1)/4, (4/7)(2θ−1)}}.
Comparison: E1 (3θ−1)/4; SL ν_A = min{1/2, 2θ−1, (5θ−2)/6} ≤ (3θ−1)/4; SL ν_B = min{1/2, (17/33)(2θ−1)} < (4/7)(2θ−1).
  New exactly on θ ∈ (9/11, 67/68): (4/7)(2θ−1) on (9/11, 15/16], then the cap 1/2 on [15/16, 67/68) (SL's ν_B reaches 1/2 only at θ = 67/68).
  θ ≥ 67/68: no change (everything at the cap). Lifting the cap needs Conrey's ν ∈ [1/2, 4/7) main-term evaluation and his K-branch
  treatment in the window (θ = 1 style, X0 ≤ 1 there): plausible, not done; it would give ν_K = (4/7)(2θ−1) up to 4/7 at θ = 1.
Conrey's own Gaussian-window proposition (SL l. 190–192: valid for 1 − θ < 1/7 − ν/4, i.e. ν < 4θ − 24/7; Baluyot's factor T^{7δ/2}, Δ = T^{1−δ})
  beats (3θ−1)/4 only for θ > 89/91; the √X0 separation here replaces Conrey's (T/Δ)^{7/2}.

## 5. Why nothing is gained for θ ≤ 9/11 (obstruction, exact)
At the binding block (P = Q = y, N = N0 = y²T/H²) the DI/Weil bounds beat E1's absolute bound iff ν > 4θ/9:
  T2/E1 = N0^{1/2}y^{−1/8}/(y²T^{1/2}H^{−3/2}) = T^{θ/2 − 9ν/8}; at ν = (3θ−1)/4 this is < 1 iff θ > 9/11.
The loss is exactly the Bessel phase X0 ≍ T/H of the short window: Λ = √X0 = (T/H)^{1/2} charged by the separation on DI's spectral term T2
(and on every separated bound), and the conductor p·X0 = pT/H (instead of p) of the type I incomplete Kloosterman sums (Lemma I).
This phase is intrinsic: pre-Voronoi it is (kn/(hm))^{iT}, the frequency T·log(kn/(hm)) of a window of width H centred at T (total ≍ T/H over
the shifts |hm − kn| ≲ hm/H); at θ = 1 it is O(1) and Conrey's 4/7 is recovered.
Input class of the barrier: Mellin (or Fourier) separation of the Bessel weight, then any bilinear/trilinear Kloosterman bound (DI Thm 12/Lemma 1,
BC Thm 1, Watt-type refinements) applied with arbitrary outer coefficients (they all contain the Cauchy–Schwarz diagonal (AN)^{−1/2}), plus Weil
for type I: with DI Lemma 1 the limit is exactly max{(3θ−1)/4, (4/7)(2θ−1)} (exponents.py; the balanced T2 decides; BC gives (17/33)(2θ−1)).
Heuristic targets (not derived; lesson 162): (a) even if every separated DI term lost no √X0 (Kuznetsov with the Bessel weight inside; its
transform lives on spectral parameters ≲ (T/H)^{1/2}), the unseparated type I piece a = 1 of Vaughan (q = k smooth, D = Q) has Lemma I bound
y^{3/2}T^{1−2θ}, so this decomposition stays ≤ (2/3)(2θ−1), which beats (3θ−1)/4 only for θ > 5/7; (b) below 5/7 even that needs
type I sums averaged over n and over the second mollifier variable (Kloosterman sums S(nā, r; p) with r ≍ pX0/D, smooth only after a second
decomposition of p). Near θ = 1/2 every route of this class gives ν → 0 (all three thresholds ∝ 2θ − 1), while E1 keeps 1/8.


## 6. Falsification (own derivation)
Smallest gain case θ = 0.82 (just above 9/11), ν = 0.3657 < (4/7)(0.64) = 0.365714; worst block P = Q = y, N = N0 (n0 = 2ν + 1 − 2θ = 0.0914):
| term | exponent (must be < 0) | value |
|---|---|---|
| E1 absolute bound (P = Q = y) | 2ν + 1/2 − 3θ/2 | +0.0014 ✗ (E1 alone fails: needs ν < 0.365) |
| DI spectral T2, A = B = y^{1/2} (N0 < A ✓: 0.0914 < 0.1829) | (1−2θ)/2 + 7ν/8 | −0.000013 ✓ (margin = 7/8·(ν_C − ν)) |
| DI diagonal T1, larger factor inside | (1−2θ)/2 + 3ν/4 | −0.0457 ✓ |
| DI T2, smaller factor inside, α = ν/2 | 1 − 2θ + 13ν/8 | −0.046 ✓ |
| DI T3 = N^{3/4}(N+A)^{1/4}V^{−1/2} | n0·3/4 + max(n0,α)/4 − ν/2 | −0.069 ✓ (A = y^{1/2}) |
| type I (Lemma I), a ≤ U1 = y^{1/4} | 1 − 2θ + 7ν/4 | −0.0000 (−0.000025) ✓ |
| type I first term Q T^{1/2}H^{−3/2} | ν + 1/2 − 3θ/2 | −0.364 ✓ |
| K-branch, residues/diagonal lemma | ν < 1/2 | ✓ |
| Taylor (11), compact windows | J, D fixed with (2θ'−1)J > 1+ν−θ' | ✓ (J ≥ 2) |
exponents.py scans every block, N ≤ N0, both moduli, both orders, U on a grid: worst exactly at this block (exponents.log θ = 0.8182, 0.85).
θ just above 1/2 (θ = 0.51, ν = (3θ−1)/4 = 0.1325, E1's limit): DI T2 exponent (1−2θ)/2 + 7ν/8 = +0.106 > 0 and type I 1 − 2θ + 7ν/4 = +0.212:
  the Kloosterman bounds are far worse than E1's absolute bound (exponent −0.0000 at ν → 0.1325); ν_C(0.51) = 0.0114. No gain; E1 unchanged.
Numerics (falsify.py/.log, peak < 50 MB):
- F1 (separation cost): Φ = e^{iX0√x}·bump on x ∈ (1,2): (∫|Φ̃(it)|dt/2π)/sup|Φ| over √(1+X0) = 1.25, 0.38, 0.20, 0.12, 0.097, 0.092 at X0 = 0, 10, 40, 160, 640, 2560:
  → constant: Λ ≍ √X0 confirmed (not smaller).
- F2 (Lemma I): p ∈ {1009, 4001}, D ∈ {p^{3/4}, p, 4p}, X0 ∈ {0, 30, 300, 3000}, 12 random c each: max|S|/(D/p + √(p(1+X0))) ≤ 0.19 in all 24 cases ✓.
  Not sharp for single sums: |S| ≈ √(Σw) (random-walk size), i.e. the r-sum of Kloosterman sums after Poisson cancels; exploiting it needs
  averaging over the modulus p (a mollifier variable) — the open route (b) of §5.
- κ(ν) (kappa_opt.py, h_table.py, mpmath dps 30 re-evaluation, sinh-P and its degree-15 Taylor polynomial agree to 12 digits):
  κ(1/2) = 0.3635010, κ(0.56) = 0.4012710 (Conrey's one-piece 4/7 value 0.4088 is the known limit), κ(0.12) = 0.085244 (> E1's frozen-detector κ1 = 0.084478 ✓).


## 7. Review
Reviewer 1 (14 Nov, Claims A–D pasted with exact values): GAP. Points and resolution:
1. "Vaughan's c2 has a = m1m2 ≤ U², not a ≤ U; the range (U, U²] is omitted" — it is not omitted (§3, exponents.py piece 2): q = a·k with the whole
   a = m1m2 as one variable, c(a) = Σ_{m1m2=a; m1,m2≤U} μχ0(m1)μχ0(m2) (|c(a)| ≤ d(a)), k free (χ0(k), smooth). a ≤ U1: Lemma I (D = Q/a);
   a ∈ (U1, U²]: bilinear (a, k), DI Lemma 1 with a inside (C(a,ℓ) = c(a)a^{s}·α_ℓℓ^{−s}ω'(ℓ/N'), independent of b = k and v = p), k outside.
   At P = Q = y (V = y ≥ A² since a ≤ y^{1/2}): T1 = (N0/a)^{1/2} ≤ (N0/U)^{1/2}; T2 = N0^{1/2}a^{1/4}y^{−1/4} (a > N0), ≤ N0^{1/2}y^{−1/8} at a = U² = y^{1/2};
   T2 = N0a^{−1/4}y^{−1/4} (a ≤ N0), ≤ T^{1−2θ+27ν/16}. Max = T^{(1−2θ)/2+7ν/8}: the same 4/7 condition. No truncated divisor weight appears.
2. "R_triv·Λ = N√(1+1/X0)": yes; here X0 ≥ √(T/y²) → ∞ (ν < 1/2), so = N(1+o(1)).
3. "C(a,ℓ) must not depend on b, v; keep (ab,v) = 1": holds (separated twists n^{−s}, a^{s}, b^{s}, p^{s}; the weight P(log(y/(gab))/log y) splits by the
   binomial theorem into O(1) products f(a)g(b); the n-block cutoff ω'(n/N') depends on n only; (p,q) = 1 is (ab, v) = 1).
4. "Lemma I and Claim A need derivative bounds, not only the amplitude size": they hold and are E1's audited Steps 2–4: in the Hankel form
   Φ(x) = Σ_± ∫G(ξ)ξ^{−z/2}X^{−1/2}a_±(X)e^{±i(X−π/4)}dξ, X = 4π√(ξx), ∫|G| ≪ (T/H²)^jH, X^k∂_X^k a_± = O_k(1), x∂_x X = X/2; so for each ξ the
   amplitude is a symbol in x (x^k∂_x^k amp ≪_k |amp|) and the phase has x²∂_x²X = −X/4 ≍ X0 (Claim A: stationary phase in log x; Lemma I: in
   k, ∂_k²X = (3/4)X/k² ≍ X0/D², amplitude of bounded variation after the smooth cutoff ω(k/D)); apply both inside the ξ-integral (Fubini).
5. Arithmetic of the balanced block and of D: confirmed correct by the reviewer.


## 8. Transfer principle (why every separated Kloosterman range is ∝ 2θ − 1) — the barrier
After separation (§2), a block's relative error is R = N·G(N, P, Q; A, B, a) with G the relative gain of the arithmetic bound; G involves only the
lengths, and N ≤ N0 = PQ/T' with T' := H²/T = T^{2θ−1}. At θ = 1, height T' (X0 ≪ 1, no separation needed) the same block has the same
R = N·G with the same N0 = PQ/T'. The separated twists n^{−it}p^{it}q^{it} keep every coefficient class (arbitrary; μ-structured, since Vaughan
applies to μ(q)q^{it}). Lemma I transfers too: R_triv·√X0 = N, so R_I = N·P^{1/2}a/Q + (non-binding R_triv/P). Hence:
  Transfer principle. "Separate (Mellin/Fourier in n/(pq)), then bound the t-slice uniformly in t" applied to SL's short-window J/Y dual at
  (T, H = T^θ) is the θ = 1 off-diagonal problem at height T' = H²/T (t = 0 is one slice). A method with θ = 1 range ν < ν1 gives ν < (2θ−1)ν1.
  Checks: BC (ν1 = 17/33) → SL's ν_B = (17/33)(2θ−1) ✓; DI/Weil (ν1 = 4/7) → (4/7)(2θ−1) ✓ (exponents.py); trivial (ν1 = 1/2) → (2θ−1)/2,
  below E1's unseparated (3θ−1)/4 by exactly the separation cost ((3θ−1)/4 − (2θ−1)/2 = (1−θ)/4).
  Separated methods beat E1 iff ν1(2θ−1) > (3θ−1)/4 ⟺ θ > (4ν1−1)/(8ν1−3): ν1 = 4/7 → 9/11; 17/33 → 35/37; ν1 = 1 → 3/5.
BARRIER (input class: separate the Bessel weight, then any θ = 1 off-diagonal technology — Weil, Kuznetsov/DI Thm 12 and Lemma 1, BC trilinear,
Watt-type refinements, even a hypothetical θ = 1 range ν1 = 1, the length beyond which the θ = 1 dual carries main terms for general coefficients
[BCR's expected limit]): no member improves E1's ν < (3θ−1)/4 for θ ≤ 3/5, and with the known inputs (best ν1 = 4/7, Conrey 1989, unchanged since)
none improves it for θ ≤ 9/11. As θ ↓ 1/2 the effective height T' = T^{2θ−1} → T^0, so every such range → 0 while E1 keeps 1/8.
Escape routes (not in the class): keep the Bessel phase inside the Cauchy–Schwarz/Kuznetsov step (exploit the t-average: S(t) for |t| ≍ X0,
a hybrid Kloosterman–archimedean large sieve), or E1-type absolute bounds sharpened by cancellation in n alone (Lemma I numerics: the
r-sums of Kloosterman sums do cancel, |S| ≈ √D, but using it needs an average over the modulus p). Neither is available; heuristic targets only.
Effect on zeros: none. Wherever ν improves (θ > 9/11), Wang's c(θ) ≥ c(9/11) = 0.5103 exceeds h(θ;κ) and (c+2κ)/3 for every κ ≤ κ(1/2) = 0.3635,
so SL/E1's max{0, c, (c+2κ)/3, h} is c(θ) before and after (h_cols.log); margin ≥ 0.0648 on [9/11, 1] even at κ = 0.409 (dominance.log);
where h(θ;κ) is the maximum (θ ≲ 0.6), no Kloosterman gain exists.

## 9. Prior art (window versions of the off-diagonal; scans 14 Nov 2026)
- Conrey 1989 Gaussian-window proposition: ν < 4θ − 24/7 (SL l. 190–192; Baluyot's T^{7δ/2}, MAP.md); Das–Pujahari arXiv 2104.10243 (2021) Thm 1
  "walk through Conrey's estimates": their bound (TeX l. 487–494) T^{3−5a/2+7ν/8} + T^{5/2−5a/2+7ν/4} gives ν < min{4a − 24/7, 2a − 10/7} = 4a − 24/7
  for a < 1, but Thm 1 states ν < 2(7a−5)/7 = 2a − 10/7 (the larger one) — a misstatement; either way below (4/7)(2a−1) for every a < 1.
  Their Cor 5 (k = 0): ν < min{a − 1/2, 1/4} (elementary diagonal range).
- SL v0.03 (latest; versions 0.01–0.03 checked on Zenodo): ν_A, ν_B = min{1/2,(17/33)(2θ−1)} (BC), E1: (3θ−1)/4 (Steuding's line).
- Tang arXiv 2608.14852: reciprocity for one twisted pair (p, q primes) on windows T^δ, H²/max{p,q} > T^{1+ε}; no mollifier average.
- arXiv monthly listings math.NT 2026-10, 2026-11 (.scratch/lst.py, 28 title hits; none on short-window mollified moments); arXiv API:
  abs:"short intervals" AND abs:mollifier; Kloosterman AND "short interval" AND zeta AND moment; "twisted second moment" AND short;
  "mean square" AND "Dirichlet polynomial" AND "short interval"; scout's 8 arXiv queries (BCR 1411.7764, 1609.02539, 2508.11108, 1403.5786: θ = 1 only);
  Zenodo '"short intervals" mollifier Kloosterman': only SL's records. (4/7)(2θ−1) in windows T^θ not found.

