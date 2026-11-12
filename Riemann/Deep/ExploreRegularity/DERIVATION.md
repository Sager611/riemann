# explore43 derivation log (appends ≤ 60 lines)

## §1 Bridges measured (14:21–15:00 UTC)
Field tools: (a) Frieze–Kannan weak regularity, (b) regularity + counting lemma, (c) dense model theorem (Green–Tao, RTTV,
Gowers Hahn–Banach), (d) transference: Green's restriction theorem for majorized functions (Green 2005, Green–Tao 2006).
Node quantities: N8 small gaps D(λ) (record λ > 2/3, Theorem U; pair-LP floor 0.775; DPP cap 1/(1+θ)).

### §1.1 Pair kernel of zeros as a graph (brief's N8 bridge): scale barrier (partA.py, partA.log)
Unfold zeros to density 1. Graph G_s: edge if |x − y| ≤ s. It has ~2sN edges, so its dense cut norm is O(s/N) → 0: every
graphon / Frieze–Kannan limit of G_s is 0 for ζ, GUE, Poisson, a lattice and DH alike (a). In the sparse normalisation
(divide by 2sN) a dense model g (cell counts at resolution s) misstates G_s by
 (i) the atomic term: g ⊗ g counts the N diagonal pairs (x, x), relative size exactly 1/(2s);
 (ii) the cell fluctuations, relative size ≈ √Σ²(s)/s, Σ² the number variance.
Measured (partA.log; ζ 22,462 zeros 100 < t < 2·10⁴; GUE tridiagonal bulk 9,993; Poisson 10⁴; DH on-line 304, t < 420):
| s | atomic 1/(2s) | ζ | GUE | Poisson | DH | lattice |
|---|---|---|---|---|---|---|
| 0.5 | 1.000 | 1.044 | 1.061 | 1.412 | 1.034 | 1.000 |
| 1 | 0.500 | 0.567 | 0.588 | 0.986 | 0.523 | 0 |
| 4 | 0.125 | 0.155 | 0.175 | 0.506 | 0.190 | 0 |
| 16 | 0.031 | 0.038 | 0.049 | 0.257 | 0.078 | 0 |
So a dense model is ε-accurate for G_s only for s ≳ 1/ε (both terms); gap statements live at s = λ < 1, where both errors
exceed 1/2. At s ≫ 1 the model sees only Σ²(s): ζ ≈ GUE (Berry saturation at T ~ 10⁴: Σ²(16) = 0.38), Poisson far off,
DH (low height) less rigid (Σ²(16) = 1.57). Nothing there is a node quantity.
(b)+(c): a counting-lemma or dense-model count of a band-limited kernel φ (φ^ ⊂ [−θ, θ]) is a band-limited n-level
correlation sum; the a-discrete sine DPP, a = 1/(1+θ), shares all of them (Lagarias–Rodgers 1907.03391; explore34
Theorem A) and has no gap below a. Hence no such argument proves D(λ) > 0 for λ < 1/(1+θ): 7/11 at θ = 4/7. Known
(explore34's collapse); the regularity language adds no input. The dense model at scale 1/θ re-derives Montgomery's
F(α): ‖g‖₂² − ‖g‖₁ is the excess ∫φ^ min(|α|, 1), i.e. a re-encoding (R044's CIRCULAR).
=> The live question on N8 is Theorem U's open input (below 2/3): ∫_T^{2T}|R_T|^{2+δ} ≪ T for rough Möbius
polynomials of length T^θ, θ ∈ (1/2, 4/7) (DPCUncond PROOF Step 5a, 6e). Tool (d) is the field's L^q (q > 2) engine.

## §2 Chosen transfer: majorant transference → Theorem U's θ < 1/2 cap (N8) (15:05 UTC)
Tool (Green–Tao transference, tool (d)): do not bound an unbounded weight by a mean-value theorem; dominate it pointwise by a
majorant whose correlations are known, and transfer the moment bound. Theorem U bounds |R_T| only through Montgomery–Vaughan
fourth moments (PROOF 5a, 6e), the sole source of θ < 1/2. Majorant here: the explicit formula makes DPC's rough prime sum an
exact smoothed sum over zeros (DPC 2610.28527 Lemma 5, TeX l.609–671, under RH): with F̂(u) = ρ(u)/u, F(x) = ∫F̂(u)e^{2πixu}du,
  S_Λ(t) := L⁻¹Σ_n Λ(n)n^{−1/2−it}F̂(log n/L) = −Σ_ρ F(L(γ_ρ − t)/2π) + O(T^{−A}),  γ_ρ := −i(ρ − 1/2) = γ − i d_ρ,
(gamma integral ≪ T^{−A} because F̂ vanishes near 0; pole terms ≪ T^{−A}). Off-line zeros enter at complex argument:
|F(x + v − iLd/2π)| ≪_M (1 + |v|)^M(1 + |x|)^{−M} e^{b′Ld}(1 + Ld)^M, b′ = sup supp F̂.
Test (t1.py a, t1a.log): on the 22,491 zeros, T = 10⁴, supp ρ ⊂ (0.05, 0.6)L and (0.05, 1.0)L (67 and 1,229 prime powers),
40 random t: |S_Λ − zero side| ≤ 4.1·10⁻⁹ (|S_Λ| ≈ 1). The majorant identity is exact.

### Lemma A (DPC Lemma 5 without RH). κ > 0 with N(σ, T) ≪ T^{1−κ(σ−1/2)}log T uniformly for σ ≥ 1/2. DPC's R_T
(c_n = μ(n)n^{−1/2}f(log n/L)∏_{q|n}ρ(log q/L), f smooth, f = 0 on [θ, ∞), supp ρ ⊂ (η/2, 2b)): ∫_T^{2T}|R_T|^q ≪ T whenever qθ < κ.
Sketch. (i) DPC: R_T = Σ_{r≤k}(−1)^r∫f̂(v)e_r(t, v)dv, k = ⌊2θ/η⌋, e_r a Newton polynomial in modulated power sums S_j(t, v);
S_j, j ≥ 2, are O(1). (ii) Substitute S_1(t, v) = −Σ_ρ F_v(L(γ_ρ − t)/2π) + O(1) and multiply out. (iii) Do the v-integral first:
a product of r zero kernels is ∫∏F̂(u_i)e^{2πiΣ x_i u_i}e^{LΣ d_i u_i}·f(Σu_i + ·)du, and f(Σu_i) = 0 for Σu_i ≥ θ, so the
joint off-line weight is T^{θ max d_i}, not DPC's crude T^{2rbd} from |e_r| ≪ (1 + |S_1|)^r. Hence, at t = x − iy,
  |R_T(x − iy)| ≪ (1 + n(x))^k [1 + Σ_{d_ρ>0} T^{θ(d_ρ − y)_+} w_ρ(x)],  w_ρ(x) = (1 + Ld_ρ)^M(1 + L|γ − x|)^{−M},
n(x) = Σ_ρ w_ρ(x); the mirror factor R*(x − iy) (s = 1/2 − y + ix) has T^{θy}(1 + n)^k[1 + Σ T^{θd_ρ}w_ρ]. (iv) On-line part:
n has every moment ≪ T (Fujii F2 for windows ≤ 100/L, Minkowski beyond, F1 for |γ − x| > 1): no condition on θ.
Off-line part: bins d ∈ [a/L, (a+1)/L), count σ_aTL ≤ CTLe^{−κa}: ∫(Σ T^{θd}w_ρ)^q ≪ Σ_a σ_a T(a+2)^C e^{qθ(a+1)} ≪ T iff qθ < κ. ∎(sketch)
Test (t1.py b2, t1b2.log; η = 0.05, f = (1 − x/θ)², 20,000 random t ∈ [T, 2T], normalised by (Σc_n²)^{q/2}):
| T | θ | terms | M₄ sampled | MV diagonal Σ|b_n|²/G² | MV off-diagonal factor 1 + T^{2θ−1} | M₆ |
|---|---|---|---|---|---|---|
| 10⁴ | 0.60 | 154 | 1.580 | 1.572 | 7.31 | 3.34 |
| 10⁶ | 0.55 | 1,213 | 1.543 | 1.531 | 4.98 | 3.22 |
| 10⁶ | 0.60 | 2,421 | 1.607 | – | 16.85 | 3.65 |
| 10⁸ | 0.55 | 10,184 | 1.578 | – | 7.31 | 3.50 |
| 10⁸ | 0.60 | 25,574 | 1.738 | – | 40.81 | 4.41 |
The fourth moment equals its diagonal for θ > 1/2; MV's T^{2θ−1} never appears (RH holds here, so this tests only the on-line part).

### Lemma B (pointwise hairpin weight, replaces 6e). For ρ ∈ S_a (depth d), slit point x_ρ = γ ∓ h/2, 0 ≤ y ≤ d: the self
terms give |R_T(x_ρ − iy)| ≲ T^{θ(d−y)}, |R*(x_ρ − iy)| ≲ T^{θ(d+y)}, product T^{2θd} ≤ e^{2θ(a+1)}, times (1 + n(x_ρ))^{2k} and
cross terms from nearby off-line zeros (each reorganises under its own deeper bin). With 6d's Hölder over S_a:
  𝓗_{T,a} ≪ T(a+2)^C e^{2θ(a+1)} σ_a^{1−1/w}   (6e had T(a+2)^C e^{θ(a+1)} σ_a^{E}, E ↑ 1/2: the same threshold).
6f: Σ_a (a+2)^C e^{2θ(a+1)} min(2ε, Ce^{−κa})^{1−1/w} ≪ ε^{(1−1/w) − 2θ/κ}log^C(1/ε) → 0 iff 2θ < κ(1 − 1/w).
Step 7 (E_T term, Hölder (q, r, ·) with q = 2 + δ): needs only (2 + δ)θ < κ.
Global moments + sparse Hölder instead of Lemma B give E = 1 − 2/q with qθ < κ, i.e. θ < κ/3 = 8/21 < 1/2: worse than MV.
So the pointwise slit bound is what removes θ < 1/2.

### Claim (Theorem U′, proof to be written). For every θ < min(1, κ/2), Theorem U holds at length T^θ; hence D(λ) > 0 for every
λ > 1/(1 + min(1, κ/2)). Jutila κ = 1 − δ: 2/3 (no gain). Conrey κ = 8/7 − δ: λ > 7/11 = 0.63636. Density hypothesis near the
line (κ = 2): λ > 1/2, DPC's RH value.

## §3 Novelty check of §2, at-scale tests, verdict (15:10–15:20 UTC)
- Our files: dpcbeyond (Riemann/Deep/DPCBeyond, commit 331db28, 12 Nov 2026) already proves this: Lemma Z = Lemma A (coupled cutoff
  f(Σlog q_i/L) inside the r-fold integral; joint off-line cost X^{max d}), Lemma M (∫|R_T|^q ≪ T for q < κ/θ), Step 6e′ = Lemma B
  (pointwise slits, e^{2θa} per bin), Step 7′ (any q > 2), Theorem U′: λ > max{2/(2 + κ), 1/2} = 7/11 with κ = 8/7; falsified in
  R454 (Lemma Z on zeros, T = 2000, 6000: ratio 1.47, 1.53), referees CORRECT/GAP-repairable. §2 is an independent re-derivation
  (found before reading DPCBeyond): same mechanism, same threshold 2θ < κ, same corollary (density hypothesis ⇒ λ > 1/2). Not new.
  Extra here: global q-th moments + generalised Hölder (p₁, p₂, p₃) on the slits give E = 1 − 2/q with qθ < κ, so θ < κ/3 = 8/21:
  the pointwise slit bound is necessary for any gain over Montgomery–Vaughan (dpcbeyond: Lemma Z gives ∫|R_T|⁴ ≪ T only for θ < κ/4).
- At-scale check of Lemma Z (t2.py, t2.log, t2_eta01.log; θ = 0.55, f = 1 − S(x/θ) smooth, ψ, ω as DPC; RH range so W = On):
  max |R_T(t)|/(1 + On(t))^k, On = Σ_ρ min(1, (L|γ − t|)^{−A}):
| height | T | terms | η = 0.2 (k = 2), A = 2 / 4 | η = 0.1 (k = 5), A = 2 / 4 | max |R_T| (η = 0.2 / 0.1) |
|---|---|---|---|---|---|
| t ∈ [10⁴, 1.99·10⁴] (495,000 pts) | 10⁴ | 35 / 56 | 0.987 / 1.208 | 1.123 / 1.852 | 1.28 / 2.15 |
| 34 Gram windows, t ≈ 10⁸ (4,080 pts) | 6.7·10⁷ | 2,342 / 3,926 | 0.905 / 1.134 | 0.950 / 1.737 | 1.26 / 2.46 |
| 34 Gram windows, t ≈ 9.8·10¹⁰ | 6.5·10¹⁰ | 68,573 / 125,854 | 0.898 / 1.163 | 0.873 / 1.916 | 1.36 / 2.42 |
  The constant does not drift over seven decades of T. Weak test: |R_T| ≤ 2.5 here, and off-line zeros (the lemma's content) cannot
  occur in the RH-verified range; R454's planted-mechanism test covers that part.
- Adversarial. The conclusion D(λ) > 0 is RH-blind: D(7/11) = 0.18–0.21 for ζ, 0.208 for the sine process, 0.118 for DH's on-line
  zeros (explore34 table). The inputs are ζ-specific: Euler product (main term at Re s = 3/2) and a uniform density slope κ > 0;
  DH has N(σ, T) ≫ T for 1/2 < σ < 1 (Voronin), so κ_DH = 0 and Theorem U′ says nothing for DH. No RH information (as for all N8
  gap results); the gain over 2/3 is Conrey's mollifier length 4/7 > 1/2 imported through κ = 2·(4/7).
Verdict for the field: (a)–(c) give no input beyond band-limited correlations; (d) reproduces dpcbeyond. NO IDEA (new).
