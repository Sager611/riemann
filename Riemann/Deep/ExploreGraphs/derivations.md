# explore12 derivations: quantum graphs → RH tree (2026-11-05)

Notation: k is the spectral variable of a graph; for ζ we use the same letter, s = 1/2 − ik, so the critical line is Im k = 0, Re s > 1/2 is Im k > 0 and Re s > 1 is Im k > 1/2. Fourier convention: ĥ(x) = ∫ h(k) e^{ikx} dk.

## 0. Objects, normalisations, checks

**Graph.** Directed bonds b with lengths ℓ_b > 0, L = diag(ℓ_b), D bonds, unitary bond-scattering matrix S (S[b',b] ≠ 0 only if b ends where b' starts). Operator: H = −i d/dx on ⊕ L²(0, ℓ_b) with ψ_{b'}(0) = Σ_b S[b',b] ψ_b(ℓ_b). It is self-adjoint iff S is unitary (Exner, "Momentum operators on graphs", arXiv 1205.5941). Spectrum: the zeros of det(I − S e^{ikL}). The Kirchhoff Laplacian gives the same equation with S = σΠ (σ the vertex matrix on half-edges, Π the bond reversal); there S is k-independent iff σ is a Hermitian unitary. Every computation below uses real orthogonal S, so the spectrum is symmetric under k → −k.

**Counting without tracking.** The eigenphases of U(k) = S e^{ikL} increase at total rate Tr L, so
N(k) = (k·TrL + Θ₀ − Σ_j arg_{[0,2π)} μ_j(k)) / 2π
is an exact integer, where μ_j(k) are the eigenvalues of U(k) and Θ₀ is fixed by N(0+) = 0. Eigenvalues are located by bisection on its jumps (g1.py `spectrum`).

**Trace formula** (Roth; Kottos–Smilansky; here for the first-order operator):
d(k) = TrL/2π + (1/π) Re Σ_{m≥1} Σ_{closed walks γ of m steps} (ℓ_γ/m) A_γ e^{ikℓ_γ},
where A_γ is the product of the S-entries along γ. The weight at length x is W(x) = Σ_{ℓ_γ = x} (ℓ_γ/m) Re A_γ. For h with ĥ supported in [−X, X]:
Σ_n h(k_n) = (TrL/2π) ĥ(0) + (1/π) Σ_{x<X} W(x) ĥ(x).

**Zeta** (Weil's explicit formula, same normalisation):
Σ_γ h(γ) = h(i/2) + h(−i/2) + (1/2π) ∫ h(t) (Re ψ(1/4 + it/2) − log π) dt + (1/π) Σ_n W_ζ(log n) ĥ(log n), with W_ζ(log n) = −Λ(n)/√n.

**Test function.** ĥ = Λ_a ⋆ Λ_a with Λ_a(x) = (1 − |x|/a)₊ and a = X/2, so supp ĥ = [−X, X]. Then h(k) = (a²/2π) sinc⁴(ka/2), an entire function decaying like k⁻⁴.

**Checks.**
- The graph trace formula on the Kirchhoff flower with petals log 2, log 3, log 5 agrees to 1.4·10⁻⁸ (g1_out.txt). The offset is the tail of eigenvalues beyond k = 400.
- Weil's formula with 377 zeros up to height 650 agrees to 2·10⁻⁹ (g3_out.txt).

## 1. What makes graph spectra real, and the ζ counterpart of each ingredient

| graph ingredient | consequence on graphs | ζ counterpart | status for ζ |
|---|---|---|---|
| (a) S unitary | H self-adjoint | none known (Hilbert–Pólya) | open; equivalent to RH |
| (b) for Im k > 0, ‖U(k)‖ ≤ e^{−ℓ_min Im k} < 1, with U finite | det(I − U) ≠ 0 in the upper half plane | the Euler product as det(I − U) with U = ⊕ p^{−1/2} e^{ik log p} (§6) | U is a contraction for Im k > −1/2, but det(I − U) exists only for Im k > 1/2 (Re s > 1, Landau's theorem: positive amplitudes); pushing it to Im k > 0 is RH |
| (c) det(I − U) = det(−U)·conj det(I − U) on R | a real "Hardy function", like Hardy's Z | Z(t) = e^{iϑ}ζ(1/2+it) is real | holds, and also holds for Davenport–Heilbronn, so it carries no zero information (§7) |
| (d) Lee–Yang: P(z) = det(I − SZ) has no zeros on D^D or on (C∖D̄)^D | real zeros for every positive length vector (Kurasov–Sarnak 2004.05678) | Bohr lift z_p = p^{−1/2+ik} | for real k this point lies inside the polydisc (\|z_p\| = p^{−1/2}), where a stable function has no zeros; ζ's zeros live beyond the convergence boundary, and the functional equation (Γ-factor) is not a torus symmetry (explore5) |

Ingredient (b) is the whole mechanism: in finite dimension the determinant is entire, so the contraction applies all the way down to the real axis. Unitarity on R then mirrors the zero-free upper half plane to the lower one.

## 2. The sign, as one identity, and the archimedean place as orbits

Every term of Weil's formula has the graph shape with amplitude −e^{−x/2}.
- Prime powers: −(1/π) e^{−x/2} dψ(e^x), since Λ(n)/√n = Λ(n) e^{−x/2} at x = log n.
- Archimedean place, rewriting ψ(1/4 + it/2) by its integral representation:
  (1/2π) ∫ h(t) Re ψ(1/4 + it/2) dt = (1/2π) ∫_0^∞ [ĥ(0) e^{−2x}/x − 2 e^{−x/2} ĥ(x)/(1 − e^{−2x})] dx.
  So the archimedean place is a continuum of "orbits" with amplitude −e^{−x/2} and length density 1/(1 − e^{−2x}): this is ≈ 1/(2x) near 0, which is the log t growth of the zero density, and → 1 at infinity.

Collecting signs:
Σ_γ h(γ) = −{ −[h(i/2) + h(−i/2)] + (−arch)(h) + (1/π) Σ_{p,m} log p · (p^{−1/2})^m ĥ(m log p) }.

The braces hold a self-adjoint-type trace formula:
- orbits of length log p with amplitude A_p = +p^{−1/2} and repetitions A_p^m;
- one bound state k = ±i/2, which enters with the graph sign "−" (spectral side);
- Weyl term −arch.

The zeros are the negative of such a spectrum: an absorption spectrum (Connes). On a graph, a bound state (negative eigenvalue of an attractive δ-coupling) appears on the spectral side, so it enters against the orbit sum. ζ's pole enters with the opposite sign: it pays for the primes (§5).

Known literature workarounds:
- Kuipers–Hummel–Richter (KHR, 1307.6055, PRL 112 (2014) 070406): infinitely many uncoupled butterflies (p, m) with Σ_{d|m} d Tr S_{p,d}^{m/d} = −p^{−m/2}.
- Dyadic orbits 2^j log p with amplitude −p^{−2^{j−1}}, using ∏_j (1 + p^{−2^j s}) = (1 − p^{−s})^{−1}. KHR attribute this idea to an earlier Andreev-reflection proposal.

Both are infinite direct sums.

## 3. Window realisation F_X (finite, self-adjoint, one vertex)

**Construction.** For every p ≤ e^X:
- r_p = ⌊X/log p⌋ directed loops of length log p;
- the block T_p = p^{−1/2} R_{r_p}, where R_r is the cyclic shift on R^{r+1} restricted to the sum-zero subspace (orthogonal, eigenvalues the (r+1)-th roots of unity other than 1).

Then Tr T_p^m = p^{−m/2}((r+1)[r+1 | m] − 1) = −p^{−m/2} for m ≤ r_p.

Let T = ⊕_p T_p, a strict contraction. Take S as its Julia (unitary) dilation: S = [[T, U_J√(1−Σ_J²)], [√(1−Σ_J²)V_J^T, −Σ_J]], from the SVD T = UΣV^T restricted to the d singular values below 1. The d auxiliary loops get length ℓ_aux > X.

**Why it works.**
- A closed walk of length ≤ X never enters an auxiliary loop.
- Inside the short loops, T is block-diagonal by prime.
- So within (0, X) the only orbits are prime powers p^m with m ≤ r_p, and their weight is log p · Tr T_p^m = −log p·p^{−m/2} = W_ζ.
- There are no composite orbits inside the window; mixed walks exist only through the auxiliary loops, at lengths > X.

**Identity.** For every h with supp ĥ ⊂ [−X, X]:
Σ_n h(k_n) − (TrL/2π) ĥ(0) = Σ_γ h(γ) − [h(i/2) + h(−i/2)] − arch(h).

**Numbers (g3_out.txt).** X = log 16.5.
- Short loops: r = (4, 2, 1, 1, 1, 1) for p = 2, …, 13, total length log 720720 = ψ(16) = 13.488.
- Auxiliary loops: 10 of length 2.9. TrL = 42.49. \|SSᵀ − I\| = 3·10⁻¹⁶.
- 4395 eigenvalues in (0, 650] against 377 zeros.
- The identity holds to 1.4·10⁻⁸ at t = 0, 7, γ₁, 17.5, γ₂, 23, γ₃, 40, 60, 100. The constant offset is the eigenvalue tail.

The block traces match at m ≤ r_p and fail exactly at m = r_p + 1, where Tr = +r_p p^{−(r_p+1)/2}.

The auxiliary part could be cut to π(e^X) loops: by Horn's theorem a contraction with eigenvalues of modulus c and singular values (1, …, 1, c^r) exists.

## 4. Newton barrier

**Theorem.**
Hypotheses:
- P ∈ C[z_1, …, z_D] with P(0) = 1 and multidegree d, containing the top monomial z^d. This holds for every graph secular polynomial det(I − SZ), whose top coefficient is det(−S) ≠ 0. For a general Lee–Yang polynomial we assume it.
- f(k) = P(e^{ikℓ_1}, …, e^{ikℓ_D}), with every ℓ_b in {log m : m ≥ 2 an integer}.
- The Dirichlet coefficients of log f (its expansion at Im k = +∞) agree with those of log ζ(1/2 − ik) at all frequencies log n ≤ X. Equivalently, the trace-formula weights of f's zero set equal −Λ(n)/√n on (0, X].

Conclusion: the frequency span Σ_b d_b ℓ_b is at least ψ(e^X). Consequences:
- For a graph, TrL ≥ ψ(e^X), so the spectral density is at least ψ(e^X)/2π.
- For a Lee–Yang f, all zeros are real and their density is ≥ ψ(e^X)/2π ≈ e^X/2π.

**Proof.**
1. Let B_p be the set of variables whose length is a multiple e_b·log p.
2. By unique factorisation, a monomial z^α has frequency m log p only if α is supported on B_p.
3. Put f_p(x) = P(z_b = x^{e_b} for b ∈ B_p, z_b = 0 otherwise). Then log f_p(x) ≡ −log(1 − x/√p) mod x^{M+1}, with M = ⌊X/log p⌋.
4. Hence the polynomial f_p ≡ Σ_{j≤M} (x/√p)^j mod x^{M+1}. All these coefficients are nonzero, so deg f_p ≥ M.
5. Also deg f_p ≤ Σ_{b∈B_p} d_b e_b, so Σ_{b∈B_p} d_b ℓ_b ≥ M log p.
6. The sets B_p are disjoint; summing gives Σ_p ⌊X/log p⌋ log p = ψ(e^X).
7. Since the top monomial is present, the frequency span equals Σ_b d_b ℓ_b. ∎

**Sharpness.** F_X attains equality on its short loops.

**Extremal block.** Matching −c^m for m = 1, …, r with an r × r block forces its eigenvalues to be cω (ω^{r+1} = 1, ω ≠ 1), by Newton's identities. Then Tr = r c^{r+1} at m = r + 1. So r bonds of length log p reproduce at most r prime powers of p.

**Corollaries.**
1. No finite graph reproduces the explicit formula at all lengths. This is the reason KHR's swarm must be infinite with divergent mean density.
2. One petal per prime fails at p³. In g2 the best least-squares fit over O(12), with one petal for each p ≤ 13, leaves residual 0.69, all on N = 2, 4, 8, 16. The p = 3 and composite constraints are met to 10⁻⁸.

**Extensions.**
- *Laplacians with k-dependent vertex conditions* (remark, not a full proof). The orbit amplitudes depend on k, and S(k) → S_∞, a k-independent Hermitian unitary, as k → ∞. If the amplitudes are required to reproduce ζ's k-independent weights, then S_∞ must already do so on the same lengths. The Weyl term is TrL/2π in both cases, so TrL ≥ ψ(e^X) again.
- *Dirichlet L(s, χ).* The per-prime target is 1/(1 − χ(p)x/√p), so the bound is ψ(e^X) − Σ_{p|q} ⌊X/log p⌋ log p.

**Comparison of what the "reflection" costs.**
- A plain Dirichlet polynomial (no reflection) matching to X needs span X: Σ_{n≤e^X} n^{−1/2−s}.
- ζ's own reflection (functional equation) costs conductor t/2π at height t (Riemann–Siegel cutoff √(t/2π)).
- The Lee–Yang/unitary reflection costs conductor ≥ e^{ψ(e^X)} = lcm(1, …, ⌊e^X⌋).

So any graph or Lee–Yang model that agrees with ζ to the depth needed at height t (X ≈ ½ log(t/2π)) carries at least (t/2π)^{1/2}/2π eigenvalues per unit length. ζ carries log(t/2π)/2π zeros per unit length.

## 5. Arbitrary lengths: positivity bound

Spectral positivity with a flat Weyl term gives, for every h ≥ 0 with supp ĥ ⊂ [−X, X]:
(TrL/2π) ĥ(0) ≥ (1/π) Σ_{n<e^X} Λ(n) n^{−1/2} ĥ(log n).

Writing ĥ = u ⋆ ũ with u supported in [−X/2, X/2]:
TrL ≥ λ_max(K_X), where K_X u = Σ_{n<e^X} Λ(n) n^{−1/2} [u(· − log n) + u(· + log n)].

Fejér (u a box) gives TrL ≥ 2 Σ_{n<e^X} Λ(n) n^{−1/2} (1 − log n/X) = (8/X) e^{X/2} (1 + o(1)).

Numbers (g5_out.txt):

| X | λ_max(K_X) | ψ(e^X) (log-integer lengths) |
|---|---|---|
| log 16.5 | 3.78 | 13.49 |
| log 100.5 | 11.15 | 94.05 |
| log 1000.5 | 34.98 | 996.68 |

In ζ's formula the same budget is paid by the pole term h(±i/2) = (1/2π)∫ĥ(x)e^{∓x/2}dx, which has size e^{X/2}. A self-adjoint system can carry it only as a bound state, and §2 shows that enters with the wrong sign.

## 6. The infinite limit is the Euler product

**KHR swarm identity.** log D(k) = −Σ_{(p,m)} l_{p,m} Σ_n Tr(S_{p,m}^n) e^{iknm log p}/n. The coefficient of e^{ikM log p} is
−(1/M) Σ_{d|M} d·l_{p,d} Tr S_{p,d}^{M/d} = p^{−M/2}/M.
So D(k) = ζ(1/2 − ik) exactly on Im k > 1/2.

**Convergence.** On Im k > 1/2 the product converges absolutely: the n = 1 terms behave like Σ_p p^{−1/2−Im k} and the n ≥ 2 terms like Σ_p p^{−2 Im k}. Each factor is nonzero for Im k > 0 (contraction), but the product diverges for Im k ≤ 1/2, because Σ_p p^{−s} has abscissa 1 (Landau).

**Spectrum.** The swarm's spectrum is a dense union of periodic sequences, unrelated to the zeros. The zeros appear only after analytic continuation.

**Known equivalent.** By partial summation with π(x) = li(x) + O(x^{1/2+ε}): RH ⟺ the partial products ∏_{p≤x} (1 − p^{−s})^{−1}, renormalised by exp(−Ei((1−s) log x)), converge on Re s > 1/2 (cf. K. Conrad, "Partial Euler products on the critical line", 2005). This is the swarm product converging on 0 < Im k ≤ 1/2.

## 7. Lee–Yang, the Riemann–Siegel lookalike, fluxes

**Graph Riemann–Siegel lookalike** (Kottos–Smilansky; Band–Harrison–Joyner 1205.4214). det(I − U) = Σ_{n≤D} (−1)^n a_n with a_{D−n} = det(−U) conj(a_n). It is exact, and its zeros are real because of unitarity (monotone eigenphases), not because of the symmetry.

**ζ's frozen analogue.** g_N(τ) = Σ_{n≤N} n^{−1/2} cos(τ log(N/n)): the Riemann–Siegel main term frozen at height 2πN². It is real-valued with the same symmetry, but not real-rooted (g5_out.txt).
- Real-zero deficit on [0, 4000]: 28 % at N = 3; 15–23 % for N = 10, …, 200.
- The positivity barrier (λ_max(K)/2 > log N) certifies non-real zeros only from N ≈ 60, so it is weak. Here λ_max comes from a 600-node discretisation, and real zeros are counted as sign changes on a 0.004 grid.
- Davenport–Heilbronn has the same symmetry.

**Fluxes.** On a graph, a magnetic flux on cycle p multiplies the amplitudes by e^{imΦ_p}; for ζ this is the twist by a character χ(p)^m. Graph spectra are real for every flux. Euler products ∏(1 − f(p) p^{−s})^{−1} with generic f ∈ T^∞ have a natural boundary at Re s = 1/2: they are zero-free to its right almost surely and do not continue. So "real for every flux" has content only on the measure-zero set of characters. The known equivalent is Bagchi's strong recurrence (R36).

## 8. Adversarial tests

**(a) Davenport–Heilbronn.** g4 realises DH's weights −Λ_f(n)/√n on n ≤ 16 by the same construction:
- classes log 2 (×4), log 3 (×2), log 6, log 7, log 11, log 12, log 13, log 14, with nonzero composite weights;
- 12 short loops (length 18.79, against ψ(16) = 13.49 for ζ) and 12 auxiliary loops;
- trace identity to 9·10⁻⁸.

The realisation has a real spectrum, while DH has zeros off the line (0.8085 + 85.70i, explore7). So window realisability, and the reality of any finite realisation, is blind to RH.

**(b) Known equivalents.** Each statement reduces to a known equivalent:
- Positivity transfer on a window is Weil positivity with an extra −TrL·‖g‖² ≤ −ψ(e^X)‖g‖² (N1).
- The limit X → ∞ is the renormalised Euler product (§6).
- Real-rooted LY approximants converging to Ξ would be a Laguerre–Pólya statement (N6), but §4 forces density e^X/2π, so they cannot converge to Ξ.
