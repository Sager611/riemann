# explore9 derivations (homological stability → RH tree)

## 1. Where the arithmetic enters (from the TeX sources, e-prints read 2026-11-05)

**EVW 0912.0325 (Annals 2016).** §1.6 eq. (Gro-Lef): |X_n(F_q)| = Σ_j (−1)^j Tr(Frob | H^j_c(X_n ⊗ F̄_q, Q_ℓ)). Assume dim H^i ≤ C^i (eq. (mario)), which follows from homological stability plus an a-priori bound dim H^i ≤ C^n (eq. (moon)). "The Deligne bounds [Weil II] show that the eigenvalues of Frobenius on H^{2n−i}_c are bounded above by q^{n−i/2}", giving |#X_n(F_q)/q^n − 1| ≤ 1/(√q/C − 1) for q > C² (eq. (lef-del)). §8 (Cohen–Lenstra): the same bound (moo) with B(A) = 2C(A); the main term is Tr(Frob | H^{2n}_c) = q^n ⟺ exactly one F_q-rational component (Galois action on components; EVW II 1212.0923 computes it through the lifting invariant).

**Landesman–Levy 2503.03861.** Lemma "weak-point-bound": ‖|Z(F_q)| − q^n‖ < D_q q^{n−1/2} for q^{1/2} > K, from Grothendieck–Lefschetz, Sun's extension of Deligne's bounds to stacks, and the exponential cohomology bound (Lemma "exponential-cohomology-bound"). Main term: Frobenius eigenvalue 1 on H^0 of a geometrically irreducible component. The topological input (stable homology of every component = that of Conf_n after inverting finitely many primes) is RH-free.

**BDPW 2302.07664 + MPPRW 2402.00354.**
- Thm A (topology): Σ_i Σ_λ dim H_i(β_∞, S^λ V)(−z)^i s_{λ'} = Exp(z^{−1} Log(z + Σ_j h_{2j} z^j) − 1). Proof: scanning map, Koszul duality for Gerstenhaber algebras (Berglund). RH-free.
- Thm B (Galois action): H_i(β_∞, V^{⊗k}) is pure Tate of weight −2i + k. Proof: the stable homology is generated under the (motivic, Galois-equivariant) framed E_2 operations by degree-0 classes ("indecomposables" theorem, purely topological); degree 0 is trivial or Tate of weight 2 (symplectic form). Uses Frobenius/Galois equivariance, NOT Weil II.
- Substituting z → q^{−1}, h_j → q^{j/2} h_j gives (up to 1 − q^{−1}) the arithmetic factor R = Π_P (1 + (1+|P|^{−1})^{−1} Σ_{j>0} h_{2j}(x^{deg P})): an Euler product over monic irreducibles, i.e. the function-field CFKRS factor (Andrade–Keating §4).
- §11.3: tr_λ(g) = q^{−2g−1} Σ_d s_⟨λ⟩(Θ̄_d) = Σ_k (−1)^k Tr(Frob | H_k(H_g^{1,0}, V_λ(|λ|/2))) (Grothendieck–Lefschetz for stacks, Behrend). Θ̄_d = q^{−1/2}·Frobenius ∈ USp(2g) is unitary BECAUSE of Weil's RH for y² = d(t); V_λ is pure of weight |λ| for the same reason.
- Jimbo–Miwa / Howe duality (Sp(2g), Sp(2r)): q^{−2g−1} Σ_d L(1/2, χ_d)^r = Σ_{λ ⊆ (r^g)} tr_λ(g) dim V_{λ†}^{Sp(2r)} (pure algebra; needs the eigenvalues of Θ̄_d).
- Thm 11.3.19 (= Thm C): with a uniform stability bound θ, the moment = Q_r(2g+1) + O(4^{g(r+1)} q^{−θ(2g+1)/2}). Proof: stable degrees k ≤ θ cancel against the stable values; unstable degrees are bounded by "the Deligne bounds on Frobenius eigenvalues from [Del80]" (|eigenvalue| ≤ q^{−k/2} on H_k), Lemma "fuks" dim H_k(β_n, V) ≤ binom(n−1, k) dim V (Fuks cells of Conf_n(C)), and Σ_λ dim V_λ dim V_{λ†} = 4^{gr}; stable classes beyond the range by Thm B (weight −2k) and a Cauchy bound.
- MPPRW Prop.: θ(n) = (n − 35)/34 is a uniform stability bound (from their Thm C: relative homology vanishes for d < (n − max(2, l(λ)))/17). BDPW, intro: "Theorem C and [MPPRW] together imply [CFKRS] for q > 4^{34(r+1)}". Conjecture (BDPW 11.3.17): θ(n) = n/4 − c.
- §13 (predictions): tr_λ(g) = Σ_n T_n(λ†) + O(...); T_1 = stable homology; T_n (n ≥ 2) conjecturally = secondary / higher stability (Galatius–Kupers–Randal-Williams); unconditionally tr_λ(g) ~ T_2(λ†) for λ = (3^g). Over Q these are the DGH terms D^{1/2+1/(2n)} (BDPW Conj. 1.2.1, r ≥ 4).

## 2. CFKRS recipe for χ_{8d}, d odd squarefree (Soundararajan's family)
L(s, χ) = X(s) L(1 − s, χ), X(1/2 + a) = (8d/π)^{−a} g(a), g(a) = Γ((1/2 − a)/2)/Γ((1/2 + a)/2) (even χ, root number 1).
Average over odd squarefree d of χ_{8d}(n): 0 unless n is an odd square m², then Π_{p | m} p/(p + 1). Hence
R(β) = Π_{p odd} R_p, R_p = 1 + (p/(p+1)) (H_e(x) − 1), H_e = ½[Π_j (1 − x_j)^{−1} + Π_j (1 + x_j)^{−1}] = Σ_{j ≥ 0} h_{2j}(x), x_j = p^{−1/2−β_j}
(the same local factor as BDPW's R with |P| → p), and A(β) = R(β)/Π_{i≤j} ζ(1 + β_i + β_j) = Π_{i≤j}(1 − 2^{−1−β_i−β_j}) Π_{p odd} R_p Π_{i≤j}(1 − x_i x_j).
Recipe: M_k(α) = Σ_{ε ∈ {±1}^k} Π_{ε_j = −1} X(1/2 + α_j) R(ε α) = Π_j X(1/2 + α_j)^{1/2} Σ_ε K(ε α), K(z) = R(z) Π_j X(1/2 + z_j)^{−1/2}.
CFKRS Lemma (symplectic): Σ_ε K(εα) = ((−1)^{k(k−1)/2} 2^k / k!) (2πi)^{−k} ∮ K(z) Δ(z²)² Π z_j / Π_{i,j}(z_i − α_j)(z_i + α_j) dz. With E(z) = A(z) Π_{i≤j} ζ̃(z_i + z_j) Π_j X(1/2 + z_j)^{−1/2}, ζ̃(s) = s ζ(1 + s) (entire), the poles z_i = −z_j cancel against Δ(z²)², and at α = 0
  Q_k(L) = M_k(0) = ((−1)^{k(k−1)/2}/k!) [Π_j z_j^{2k−1}] { e^{(L/2)Σ z} F(z) Π_{i<j}(z_i − z_j)²(z_i + z_j) },  L = log(8d/π),
F = A Π ζ̃ Π g^{−1/2}. Degree check: k(2k − 1) − 3k(k − 1)/2 = k(k + 1)/2. ✓
Numerics (cfkrs.py): Taylor coefficients of F by an FFT on a polycircle (N = 32, radii 0.1·(1, 0.93, 0.86)); odd p ≤ 2·10^4 exactly; p > 2·10^4 by the p^{−2} term of log A_p,
  log A_p = t²(h_4 − h_2 − h_2²/2 − S_2/2) + t³(h_6 − h_4 + h_2 − h_2 h_4 + h_2² + h_2³/3 − S_3/3) + O(t⁴), t = 1/p, h_m = h_m(p^{−z}), S_m = Σ_{i≤j}(y_i y_j)^m,
  summed with τ_m = Σ_{p > P0} (log p)^m p^{−2}/m! (primes to 10^7, PNT integral beyond). k = 1 check: A(0) = ½ Π_{p odd}(1 − 1/(p(p+1))) = 0.42266532 (carefree constant 0.70444220 × 6/5 × ½) = 2 × 0.21133266. ✓
Independent check (check_cfkrs.py): direct recipe sum at shifts α = δ(1, √2, √3), 60 digits, linear extrapolation in δ (M_k(−α) = Π X(1/2 − α) M_k(α), so the error is O(δ)): k = 2, L = 14: 20.0389043 vs FFT 20.0389047; k = 3, L = 14: 163.80516 vs 163.80592 (O(δ²) remainder). N = 40, r = 0.12 reproduces the k = 3 coefficients to 1e−11. k = 4 (N = 24, r = 0.15, 467 s, 182 MB) and (N = 20, r = 0.18) agree to 1e−11 relative; leading coefficient 6.9159159e−12 = a_4 / 4,838,400 with a_4 = 3.346·10^{−5} (Keating–Snaith factor Π_{j≤4} j!/(2j)!).
Diaconu–Whitehead residue (1804.00690 Thm A, family χ_{2d} = χ_{8d}): (9/(256π)) 2^{1/4}(−181 + 128√2) Γ(1/4)^4 ζ(1/2)^7 Π_{p odd} P(p^{−1/2}) = −0.003420 (theirs: ≈ −0.0034). ✓

## 3. DH mean square (dh_test.py)
f = c_1 L(s, ψ) + c_2 L(s, ψ̄), ψ mod 5, ψ(2) = i, c_{1,2} = (1 ∓ iκ)/2, κ = (√(10 − 2√5) − 2)/(√5 − 1); f(0.808517182457 + 85.699348485378 i) = 6.8·10^{−13} (checked).
Recipe: ∫_0^T |L(1/2+it, χ)|² dt ≈ (φ(q)/q) T (log(qT/2π) + 2γ − 1 + 2 log q/(q − 1)); the cross term L(1/2+it, ψ) L(1/2−it, ψ) has diagonal Σ ψ(n)²/n = L(1, (·/5)) = 2 log φ/√5 from both the direct and the doubly swapped term (factor ε(ψ)²): density L(1, (·/5))(1 + ε(ψ)²). So ∫_0^T |f|² ≈ ((1 + κ²)/2)(4/5) T (log(5T/2π) + 2γ − 1 + (log 5)/2) + 2 Re(c_1 c̄_2 (1 + ε²)) L(1, (·/5)) T, ε(ψ) = τ(ψ)/(i√5) = 0.85065 + 0.52573 i, cross density 0.395675.

## 4. q-threshold (side remark, not recorded)
BDPW bound Σ_{k > θ} q^{−k/2} binom(2g, k) by 2^{2g} q^{−θ/2}. For q > 34² the terms decrease beyond θ = g/17 − 1, so the sum is ≪ binom(2g, θ) q^{−θ/2} ≈ e^{2g H(1/34)} q^{−g/34}, H(1/34) = 0.13268. The unstable part then tends to 0 as g → ∞ iff q^{1/34} > e^{0.26536} 4^r, i.e. q > 8.29·10³ · 4^{34r} (r = 3: 2.1·10^{65}), instead of q > 4^{34(r+1)} (r = 3: 7.6·10^{81}); the stable remainder needs only q > 4^{17r}(1 + o(1)). Bookkeeping on BDPW's own lemma: routine and insignificant; the slope 1/34 and the 4^{gr} = Σ dim V_λ dim V_{λ†} factor dominate. BDPW's conjectured slope 1/4 (θ(2g+1) ≈ g/2) would give, by the same count, binom(2g, g/2) ≈ e^{2g H(1/4)} = 3.079^g, hence q > 3.079⁴ · 4^{4r} = 90 · 256^r (r = 3: 1.5·10⁹; r = 4: 3.9·10^{11}): still no fixed q for all r.


## 5. Side observation: zeros of ζ in quadratic-family moments, and the squarefree count (routine, for triage)
Over Q, the condition "d squarefree" is imposed by Σ_{a² | d} μ(a); in the generating Dirichlet series this produces Σ_a μ(a) a^{−2w}·(local corrections) ≈ 1/ζ(2w), whose poles at w = ρ/2 give terms of size X^{β/2} ≤ X^{1/2} (X^{1/4} under RH). The same factor normalises the DGH multiple Dirichlet series (ζ(2w + 2s − 1) at s = 1/2). Over F_q(t), 1/ζ_{F_q[t]}(2w) = 1 − q^{1−2w} is entire: a genus-0 base has no zeros, which is why the function-field results never meet them. The known moment error terms over squarefree d (Young X^{1/2+ε}, Sono X^{1/2+ε}, Diaconu–Whitehead X^{2/3+δ}) are ≥ X^{1/2}, above these terms.
The simplest case is the squarefree count Q(X) = Σ_{n≤X} μ²(n). Evelyn–Linfoot's argument:
  Q(X) = Σ_{a≤y} μ(a)⌊X/a²⌋ + Σ_{b ≤ X/y²} [M(√(X/b)) − M(y)],
so if M(x) ≪ x^{θ+ε} then Q(X) − 6X/π² ≪ y + X y^{θ−2+ε}, and y = X^{1/(3−θ)} gives
  Q(X) = 6X/π² + O(X^{1/(3−θ)+ε})
(θ = 1/2 recovers Evelyn–Linfoot's 2/5 under RH). The tree's Lean-checked half-plane (no zeros with Re s > B = 0.8733819666) gives M(x) ≪ x^{B+ε}, hence an unconditional power saving
  Q(X) = 6X/π² + O(X^{0.47023}),  and for k-free numbers O(X^{1/(k+1−B)+ε}), e.g. k = 3: X^{0.31982} < X^{1/3}.
Before 2026 the best unconditional bound was X^{1/2} exp(−c (log X)^{3/5}(log log X)^{−1/5}) (Walfisz). Exponential sums (the Montgomery–Vaughan route) should lower 0.4702 toward B/2 = 0.4367; the X^{B/2} terms coming from possible zeros near Re s = B are the limit. With the explicit |M(x)| ≤ 1.84·10^6 x^{0.9382}, the constant is about 5·10³ at exponent 0.485, which beats √X only for X > 10^{240}: no explicit gain.
Status: routine substitution, not reviewed; novelty checked only by an arXiv API query (squarefree / k-free with error term, zero-free, RH; 2022–2026: nothing).
