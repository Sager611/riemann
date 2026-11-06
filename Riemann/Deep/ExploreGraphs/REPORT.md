NO IDEA: graph spectra are real only because a finite unitary scattering matrix makes the secular function an entire Lee–Yang polynomial, and any such realisation of ζ's explicit formula up to length X needs total bond length ≥ ψ(e^X) (density e^X/2π, not log(t/2π)/2π), tends to the bare Euler product as X → ∞, and works equally for Davenport–Heilbronn.

# Exploration explore12: quantum graphs → RH tree (2026-11-05)

Brief: `.scratch/orca/explore12-brief.md`.

Sources:
- arXiv API via `axq.py`; the abs pages and the e-print of 1307.6055 (read for its construction, then deleted);
- our TREE/NOGO (R7.1, R290, R365; explore1, 3, 5, 10, 11).

Work files: proofs in `derivations.md`; tests `g1.py`–`g5.py` with their outputs `g*_out.txt`; `run_all.sh` reruns everything sequentially and single-threaded in about 12 minutes.

## Survey (one page)

**Trace formulas.** The density is d(k) = L/π + (1/π) Re Σ_γ Σ_m L_γ A_γ^m e^{ikmL_γ}. It is exact and the weights carry no stability factor.
- Sources: Roth 1983; Kottos–Smilansky 1997/99; Gnutzmann–Smilansky review 2006; Berkolaiko–Kuchment 2013.
- The Euler characteristic can be recovered from the spectrum (Kurasov, 2002.02911).
- The first-order operator −i d/dx on directed graphs is self-adjoint for every unitary S (Exner 1205.5941). This is the form used below.

**Graphs and Riemann zeros.**
- Endres–Steiner 0912.3183: the Berry–Keating xp operator on compact graphs. They classify its self-adjoint extensions and give an exact trace formula, but its mean density cannot match ζ's.
- Endres–Steiner 1104.1364: an infinite chain with δ-couplings and lengths π/n. The Dirichlet limit has eigenvalues n² with multiplicity d(n): a k log k Weyl law and a Voronoi-type trace formula. Logarithmic density is possible, but only with short bonds whose orbits accumulate at length 0, like ζ's archimedean place (derivations §2).
- Kuipers–Hummel–Richter 1307.6055 (PRL 112 (2014) 070406):
  - an infinite swarm of uncoupled two-bond "butterflies" (p, m) with Σ_{d|m} d Tr S_{p,d}^{m/d} = −p^{−m/2} reproduces ζ's oscillating density exactly; orbits of equal length together supply the minus sign;
  - the mean density diverges, and the spectrum is unrelated to the zeros;
  - they note that coupled petals create composite orbits of length log pq; they also cite the alternative of orbits of length 2^j log p (an Andreev-reflection proposal).

**Stable polynomials and crystalline measures.**
- Kurasov–Sarnak (2004.05678): secular polynomials are Lee–Yang, so graph spectra are Fourier quasicrystals.
- Alon–Cohen–Vinzant 2303.03201; Alon–Vinzant 2307.13498 (gap laws); Alon–Kummer–Kurasov–Vinzant 2407.11184; Alon–Kummer 2507.16029; Wei 2610.10652 (positive-integer-mass Fourier quasicrystals = Lee–Yang cycles); the Gauss circle problem for Fourier quasicrystals, 2412.05485.
- Gonçalves–Kerr–Klurman–Vedana 2610.26705 (Oct 2026): a real-rooted trigonometric polynomial shares at most 2/3 of ζ's zeros on the line, proved via Kurasov–Sarnak and Alon–Cohen–Vinzant.

**Statistics, expansions, inverse problems.**
- Barra–Gaspard nlin/0011099 (Kronecker flow on the torus); 1301.1126 applies it to the fractional parts of log n.
- Band–Harrison–Joyner 1205.4214 (finite pseudo-orbit expansions: the graph Riemann–Siegel lookalike, exact); Harrison–Kirsten 0911.2509 (spectral ζ-functions of graphs).
- Berkolaiko (nodal surplus = Morse index in the magnetic flux); Gutkin–Smilansky (lengths from the spectrum); isoscattering families 2505.16621, 2610.13512; Aizenman–Warzel (trees); 2606.14931 (scale-invariant couplings, 2026).
- The discrete limit is Ihara–Bartholdi: for equilateral graphs the secular determinant is a weighted Bartholdi zeta.
- Already in our repo: R7.1 (no graph has zeta ζ), R290, R365 (quantum walks), explore5 (Fourier quasicrystals / Lee–Yang).

## (i) What forces real spectra, and the ζ counterpart of each piece (derivations §1–2)

| graph mechanism | ζ counterpart | verdict |
|---|---|---|
| S unitary ⟺ H self-adjoint | Hilbert–Pólya | equivalent to RH |
| U(k) = Se^{ikL} is a contraction for Im k > 0 and U is finite, so det(I − U) is entire and nonzero there | U = ⊕ p^{−1/2}e^{ik log p}: a contraction for Im k > −1/2, but det(I − U) = 1/ζ exists only for Re s > 1 (Landau: positive amplitudes) | extending it to Re s > 1/2 is RH (renormalised Euler product, derivations §6) |
| det(I − U) = det(−U)·conj det(I − U) on R | Hardy's Z is real | also true for Davenport–Heilbronn; no zero information |
| Lee–Yang: real zeros for every length vector | Bohr lift z_p = p^{−1/2+ik} lies inside the polydisc, where stable functions are zero-free | no transfer (explore5) |
| magnetic flux keeps spectra real | twisting by χ(p) | generic fluxes give a natural boundary at Re s = 1/2; the known equivalent is Bagchi's recurrence (R36) |

**The sign.** Weil's formula is the negative of a self-adjoint trace formula: orbits log p with amplitude +p^{−1/2}, a bound state at k = ±i/2 (ζ's pole), and Weyl term −arch. So the zeros form an absorption spectrum.

On a graph, a bound state sits on the spectral side and opposes the orbit sum. ζ's pole has the opposite sign and pays for the primes. In Weil's formula every orbit, prime-power or archimedean, has amplitude −e^{−x/2}; the archimedean length density is 1/(1 − e^{−2x}) ≈ 1/(2x) near 0.

**Can a graph realise lengths log p^k with weights −log p/p^{k/2}?** On a finite window yes (T1); at all lengths no (Newton barrier); in the limit only as the Euler product (T3).

## Transfers and tests

### T1 (N13), best transfer: finite window realisation F_X, and the Newton barrier

**Construction.**
- For each p ≤ e^X: r_p = ⌊X/log p⌋ directed loops of length log p, with block p^{−1/2}·(cyclic shift on R^{r+1} restricted to the sum-zero subspace).
- A unitary dilation with auxiliary loops longer than X.
- Within (0, X) its orbits are exactly the prime powers, with ζ's weights. So for every test function with ĥ supported in [−X, X]:
  Σ_n h(k_n) − (TrL/2π) ĥ(0) = Σ_γ h(γ) − [h(i/2) + h(−i/2)] − arch(h).

This is a finite, connected version of KHR's infinite swarm.

**Test (g3).** X = log 16.5 with the primes up to 13:
- 10 short loops of total length ψ(16) = log 720720 = 13.488, and 10 auxiliary loops of length 2.9;
- 4395 eigenvalues in (0, 650] against 377 zeros;
- the identity holds to 1.4·10⁻⁸ at ten heights t ∈ [0, 100]; the offset is the eigenvalue tail. Weil's formula alone checks to 2·10⁻⁹.

**Kirchhoff baseline (g1).** The flower with petals log 2, log 3, log 5 has weights of the wrong sign (+0.46, +0.73, +0.77, +1.07 at log 2, 3, 4, 5, against −0.49, −0.63, −0.35, −0.72 for ζ) and nonzero composite weights (+0.80 at log 6).

**Barrier (ours, routine, derivations §4).**
- *Theorem.* A graph or Lee–Yang secular function whose bond lengths are logs of integers and that agrees with log ζ(1/2 − ik) up to frequency X has frequency span (total directed length) ≥ ψ(e^X). Hence it has ≥ ψ(e^X)/2π zeros per unit length; F_X attains this on its short loops. The proof is a per-prime Newton/degree count.
- r bonds of length log p reproduce at most r powers of p; the extremal block gives +r p^{−(r+1)/2} at the next power.
- The best fit with one petal per prime (g2, S ∈ O(12)) leaves residual 0.69, all on N = 2, 4, 8, 16.
- For arbitrary bond lengths, spectral positivity alone forces TrL ≥ λ_max(prime form) ≈ (8/X) e^{X/2}. The numbers are 3.78, 11.15, 34.98 at X = log 16.5, log 100.5, log 1000.5, against ψ = 13.49, 94.05, 996.68. ζ meets this budget only through its pole.
- Laplacians with k-dependent vertex conditions inherit the bound through their k → ∞ limit (a remark, derivations §4). For L(s, χ) the bound is ψ(e^X) minus the primes dividing q.
- In conductor terms: the unitary reflection costs conductor ≥ lcm(1, …, e^X), while ζ's functional equation costs t/2π. A model that matches ζ to the Riemann–Siegel depth at height t has ≳ √(t/2π)/2π eigenvalues per unit length, against log(t/2π)/2π zeros.

Verdict: dead as a route; it yields a quantitative barrier (N13, N6).

### T2 (N1): graph positivity as a Weil certificate

Subtracting F_X's identity from Weil's form gives Q_ζ(g) = Σ_n|ĝ(k_n)|² + [pole + arch](g ⋆ g*) − TrL·‖g‖² (normalised).

The free positive part costs −TrL‖g‖² ≤ −ψ(e^X)‖g‖². This dwarfs both the archimedean term (≈ log of the bandwidth) and the indefinite rank-two pole term. So the remainder is negative definite on large subspaces, and the certificate fails beyond tiny windows.

Graph positivity buys exactly the budget that ζ gets from its pole, and pays for it with a flat density. This is a known equivalent (Weil positivity) at a worse constant.

Verdict: dead.

### T3 (N6, N2, N13): Lee–Yang / Riemann–Siegel lookalike / infinite limit

- **Infinite limit.** The swarm determinant equals ζ(1/2 − ik) exactly on Im k > 1/2 and diverges below. Its spectrum is dense and unrelated to the zeros. RH ⟺ the Ei-renormalised product converges on 0 < Im k ≤ 1/2, a known equivalent (derivations §6).
- **Riemann–Siegel lookalike.** The graph version is real-rooted through unitarity (monotone eigenphases), not through its symmetry. ζ's frozen analogue g_N(τ) = Σ_{n≤N} n^{−1/2} cos(τ log(N/n)) has the same symmetry but is not real-rooted (g5):
  - 28 % of its zeros on [0, 4000] are non-real at N = 3, and 15–23 % for N = 10, …, 200;
  - the positivity barrier certifies non-real zeros only from N ≈ 60.
- **Hurwitz route.** Real-rooted Lee–Yang approximants of Ξ would need density log(t/2π)/2π, but T1's barrier forces e^X/2π.

Verdict: dead (known equivalents plus T1's barrier).

## Adversarial tests

**(a) Davenport–Heilbronn (g4).** DH's weights −Λ_f(n)/√n on n ≤ 16 are realised the same way:
- classes log 2 (×4), log 3 (×2), log 6, 7, 11, 12, 13, 14, with nonzero composite weights;
- 12 short loops (length 18.79 against 13.49 for ζ) and 12 auxiliary loops;
- the trace identity holds to 9·10⁻⁸.

The realising graph has a real spectrum, while DH has zeros off the line (0.8085 + 85.70i). So window realisability and the reality of finite realisations say nothing about RH. The Riemann–Siegel symmetry in T3 is shared by DH as well.

**(b) Only known equivalents or zero-free facts?** Yes:
- T2 is Weil positivity with a worse constant.
- T3 is the renormalised Euler product, or the Laguerre–Pólya closure (N6).
- Graph contraction gives "no zeros for Re s > 1", the Euler-product region.
- Nothing reaches the 7/8 half-plane, let alone beyond it.

## Ranking (RH relevance × chance × gain)

| transfer | node | RH relevance | chance | gain | note |
|---|---|---|---|---|---|
| T1 window realisation + Newton barrier | N13, N6 | medium | 0 for RH | small: a clean barrier, ψ(e^X) minimal length | only new item |
| T2 positivity certificates | N1 | high | ≈ 0 | none | Weyl cost ψ(e^X) |
| T3 Lee–Yang / Riemann–Siegel / limit | N6, N2 | medium | 0 | none | known equivalents |

## Suggested NOGO line

- explore12 (Orca exploration: quantum graphs → N13, N1, N6) NO IDEA:
  - graph spectra are real because a finite unitary S makes det(I − Se^{ikL}) an entire Lee–Yang polynomial;
  - a self-adjoint one-vertex graph reproduces ζ's explicit formula exactly on lengths ≤ X (identity to 1.4e−8 at X = log 16.5);
  - but any unitary realisation with log-integer lengths needs total bond length ≥ ψ(e^X), since r loops of length log p reproduce at most r powers of p (general lengths: ≈ 8e^{X/2}/X);
  - its X → ∞ limit is the Euler product (Kuipers–Hummel–Richter 2014; convergent only for Re s > 1);
  - Davenport–Heilbronn is realised the same way.
  [ExploreGraphs](Riemann/Deep/ExploreGraphs/REPORT.md)

## Files

- `derivations.md`: normalisations, the mechanism table, the sign identity, the F_X construction, the Newton barrier with proof, the positivity bound, the swarm identity, adversarial tests.
- `g1.py`: flower graph tools (exact counting function, walk coefficients, test function); Kirchhoff flower log 2, log 3, log 5; trace-formula check. Output `g1_out.txt`.
- `g2.py C`: least-squares fit with one petal per prime (ζ and DH). Output `g2_out.txt`.
- `g3.py`: F_X, its spectrum, zeros to height 650, the window identity. Output `g3_out.txt`; zeros cached in `zeta_zeros.npy`.
- `g4.py`: Davenport–Heilbronn window realisation. Output `g4_out.txt`.
- `g5.py`: frozen Riemann–Siegel cosine polynomials and the positivity bounds. Output `g5_out.txt`.
- `run_all.sh`: sequential rerun.

Novelty of the Newton barrier was checked only against the searches listed in the survey; it is routine and is not claimed as a discovery.

RH relevance: N13 (barrier for graph/Hilbert–Pólya realisations), N6 (Lee–Yang approximants), N1 (no certificate).
