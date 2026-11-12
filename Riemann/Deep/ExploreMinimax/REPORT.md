BARRIER: by Le Cam's two-point method, no positive-margin argument from primes ≤ X rules out a near-line off-line zero pair at height t once log X < ½log(t/2π) (X below the Riemann–Siegel length √(t/2π)): moving ≤ 100 neighbouring zeros by < 1 spacing hides the pair from every Weil window c ≤ X to relative operator norm ~10⁻⁶ when log X ≤ 0.6·½log(t/2π) (→ 0 with more neighbours), never when log X ≥ ½log(t/2π) (heights 10³–10¹¹); the exact global version ζ(s)(1 + b p^{−s} + p^{1−2s}) shares every datum below p with ζ, is an Euler product with ζ's gamma factor (conductor p²) and has zeros on Re s = 0.9026 for p = b = 7; the law depends on the zero density only, so it is DH-blind and carries no RH information (node N1/N7)

# ExploreMinimax — minimax lower bounds (Le Cam, Assouad, Fano) → N1, N7 (explore36, 12 Nov 2026)

Brief: [.scratch/orca/explore36-brief.md](../../../.scratch/orca/explore36-brief.md). Proofs and conventions: [DERIVATION.md](DERIVATION.md) (D1–D6).

## 1. Transfer
- Tool: Le Cam's two-point method in its noise-free form. If two hypotheses give data within η, every procedure whose
  conclusion survives data perturbations of size η gives the same answer for both. Every certificate with a positive
  margin is such a procedure. So a statement true for one hypothesis and false for the other cannot be derived.
  Assouad = the same over a hypercube of independent flips (D6). Fano needs a noise model; with exact data it is
  counting, so it was not used.
- Bridge (D1, D2): the prime data D_X = {Λ(n) : n < X} are, by the explicit formula, exactly the zero symbol
  G_Z(y) = Σ_ρ e^{iγ_ρ y} on |y| < L = log X, equivalently the Weil form Q_Z on every window c ≤ X (N1's window c uses
  primes ≤ c). Le Cam distance at window c: η_c(Z₀, Z₁) = sup_{‖f‖=1}|Q_{Z₁}(f) − Q_{Z₀}(f)| = ‖T_{ΔG}‖ on L²(−a, a), 2a = log c
  (Q_Z is unbounded; the difference is a bounded convolution operator).
- Hypotheses:
  - Z₀ = ζ's zeros (ExploreHelio's zeros to 2·10⁴; gram_win -zeros windows at t ≈ 1.0·10⁸ and 9.8·10¹⁰, from ExploreDPP).
  - Z₁ = Z₀ with two consecutive ordinates merged into the off-line pair ½ ± δ + iγ̄ (γ̄ their mean), and the M nearest
    other ordinates moved by ε_k, real, |ε_k| ≤ 1 mean spacing. Moves ≤ 0.45 spacing keep the order (same results).
  - ε is chosen by Levenberg–Marquardt on the Hilbert–Schmidt norm. Reported: ρ = η_c(Z₀, Z₁^ε)/η_c(Z₀, Z₁^0), the
    distance after compensation relative to the bare moved pair (operator norm, Nyström).
  - Every Z₁ keeps the symmetry and Re s ≤ ½ + δ < B₂ (δ ≤ 0.37). With moves ≤ 0.45 spacing it keeps the order of the
    ordinates and N(t) up to ±1 at the pair.
- Heuristic (D4, Landau density; a conjecture, not a theorem): local compensation works approximately for L < L* = πD =
  ½log(t/2π), i.e. X < X* = √(t/2π), the Riemann–Siegel main-sum length. explore30's off-line lattices are the exact
  lattice case.

## 2. Statements and tests
**S1 (local floor, numerical).** Median ρ over 15 pairs per height, M = 100, δ = 1/(2L*); in brackets the share of pairs
with ρ < 10⁻³ (logs/stats_*_cap1.txt):

| height t | L* | X* = √(t/2π) | log X = 0.6 L* | 0.8 L* | 0.9 L* | 1.0 L* | 1.1 L* |
|---|---|---|---|---|---|---|---|
| 1.0·10³ | 2.536 | 12.6 | 4.4·10⁻⁶ [0.93] | 4.8·10⁻³ [0.33] | 2.4·10⁻² [0.07] | 0.44 [0] | 0.78 [0] |
| 1.0·10⁴ | 3.686 | 39.9 | 3.2·10⁻⁶ [0.93] | 4.3·10⁻³ [0.33] | 0.16 [0] | 0.60 [0] | 0.80 [0] |
| 1.0·10⁸ | 8.291 | 3,989 | 2.1·10⁻⁶ [0.80] | 1.9·10⁻³ [0.40] | 4.7·10⁻² [0] | 0.40 [0] | 0.66 [0] |
| 9.8·10¹⁰ | 11.73 | 1.25·10⁵ | 4.1·10⁻⁶ [1.00] | 4.1·10⁻³ [0.33] | 0.15 [0] | 0.47 [0] | 0.58 [0] |

At L* the 10th percentile is ≥ 0.12 at every height. Dependence on the number of neighbours M at t ≈ 10⁸
(logs/mdep_z1e8.txt, logs/mdep2_z1e8.txt):

| case | M = 50 | 100 | 200 | 400 | 800 |
|---|---|---|---|---|---|
| centre pair, 0.6 L* | 6.4·10⁻⁶ | 2.2·10⁻⁶ | 8.1·10⁻⁷ | 2.7·10⁻⁷ | – |
| centre pair, 0.8 L* | 1.3·10⁻⁴ | 8.6·10⁻⁶ | 3.7·10⁻⁶ | 1.5·10⁻⁶ | – |
| centre pair, 0.9 L* | 7.9·10⁻² | 3.1·10⁻² | 1.0·10⁻² | 5.7·10⁻³ | – |
| pair +150, 0.9 L* | – | 3.2·10⁻³ | 1.5·10⁻³ | 2.0·10⁻⁵ | 2.3·10⁻⁵ |
| pair +150, 1.0 L* | – | 0.52 | 0.48 | 0.44 | 0.43 |
| pair +150, 1.1 L* | – | 0.754 | 0.754 | 0.754 | 0.754 |

Below L*, ρ falls like M^{−1.5}; at and above L* it is flat. So, by the two-point lemma (D2): no argument with a positive
margin that uses only primes ≤ X excludes a near-line off-line pair at height t when log X ≤ 0.9·½log(t/2π), i.e.
t ≥ 2πX^{2.2} (the M-scans point to the limit 2πX²). The alternative only moves neighbouring zeros by < 1 spacing, which
nothing in the prime data ≤ X forbids. The construction is float64 numerics;
exact invisibility needs infinitely many moves (D3), which is conjectured, not proved.

**S2 (the separation δ).** Moves of at most one spacing also need δ·log X ≲ 1.5. Median ρ, 9 pairs at t ≈ 10⁸, M = 100
(logs/stats_z1e8_delta*.txt):

| δ | 0.3 L* | 0.5 L* | 0.7 L* | 0.8 L* | 0.9 L* | 1.0 L* |
|---|---|---|---|---|---|---|
| 0.01 | 3.6·10⁻⁷ | 2.0·10⁻⁶ | 3.9·10⁻⁵ | 1.2·10⁻² | 0.30 | 0.50 |
| 0.2 | 1.1·10⁻⁶ | 3.0·10⁻⁶ | 1.0·10⁻² (δL = 1.16) | 0.27 (δL = 1.33) | 0.54 | 0.69 |
| 0.5 | 2.0·10⁻⁶ (δL = 1.24) | 0.16 (δL = 2.07) | 0.69 | 0.81 | 0.88 | 0.91 |

δ = 0.5 lies beyond B₂ = 0.8734 and is only a probe of the δ-dependence. Combined, empirically: hidden when
log X < min(½log(t/2π), ≈ 1.5/δ). That is the brief's "X against the separation".

**S3 (Assouad, blocks).** K disjoint pairs moved off the line (δ = 0.06) in a block of 160 zeros at t ≈ 10⁸, every other
zero of the block movable (logs/multi_z1e8.txt). Median ρ:

| share of zeros off the line | 0.5 L* | 0.6 L* | 0.7 L* | 0.8 L* |
|---|---|---|---|---|
| 1 % | 3.0·10⁻⁷ | 2.7·10⁻⁶ | 3.3·10⁻⁵ | 1.1·10⁻³ |
| 5 % | 3.7·10⁻⁵ | 1.4·10⁻³ | 4.3·10⁻² | 0.35 |
| 10 % | 2.0·10⁻⁵ | 1.8·10⁻³ | 5.6·10⁻² | 0.30 |
| 20 % | 1.2·10⁻⁵ | 1.2·10⁻² | 0.32 | 0.66 |

Fewer compensators lower the threshold, roughly to (1 − share)·L*. By Assouad (D6), at log X ≤ 0.5 L* any η-robust count of
the off-line pairs in such a block from primes ≤ X has worst-case error ≥ K/2, and any η-robust guess of which pairs are
off-line has expected Hamming error ≥ K/2.

**S4 (exact global pair, D5).** For a prime p and a real b ≠ −(p + 1), take F_{p,b}(s) = ζ(s)(1 + b p^{−s} + p^{1−2s}).
- It is an Euler product with ζ's gamma factor, conductor p² and root number 1.
- It has a_n = 1 and Λ_F(n) = Λ(n) for every n < p.
- Its extra zeros are on the line iff |b| ≤ 2√p; otherwise they lie on Re s = ½ ± log|w|/log p.
- fpb_check.py (logs/fpb_check.txt):

| p | b | Re s of the extra zeros | zero sum, support 0.9 log p | zero sum, support 1.5 log p |
|---|---|---|---|---|
| 7 | 5 | 0.5 | 3.891820 | 3.771686 |
| 7 | 7 | 0.097412, 0.902588 | 3.891820 | 3.723633 |
| 11 | 8 | 0.236685, 0.763315 | 4.795791 | 4.606841 |
| 11 | 11.5 | 0.021508, 0.978492 | 4.795791 | 4.524175 |

The zero side equals the prime side (Poisson) to 10⁻¹⁶ in every case. Below support log p the sum is 2 log p·g(0) for
every b; above it, it differs. So F_{7,5} (RH ⟺ RH for ζ) and F_{7,7} (zeros at Re s = 0.9026 > B₂) agree in every
explicit-formula datum below 7, and their Weil forms on all windows c < p are Q_ζ + 2 log p‖f‖². Hence no argument that
uses only an Euler product, a degree-1 functional equation of conductor q and the coefficients below √q can prove RH,
a zero-free region or the 7/8 half-plane. For q = 1 no such pair exists (Hamburger: a Dirichlet series with ζ's functional
equation, finite order and finitely many poles is cζ), so S4 says nothing about ζ itself.
S1 is its local analogue: at height t, ζ has analytic conductor t/2π, and the threshold is √(t/2π) again.

## 3. Against the node lines (iii)
- N1 (window records c ≤ 47 for all f, 71 even). Window c has the data of primes ≤ c. Against S1-type alternatives:
  - A near-line pair (δ ≲ 1.5/log c) is hidden above height ≈ 2πc², i.e. 13,900 at c = 47.
  - Weil positivity sees isolated moved pairs (δ = 0.3) only below ≈ 13.9c, i.e. 640 at c = 47 (explore19's
    c* ≈ 0.072m + 0.73).
  - So positivity uses about 1/20 of the height range in which its own data still fix the zeros against local
    compensation. No record moves. Positivity is not a Le Cam test: its margins (2⁻⁸⁰³ at c = 47) are far below the
    η of S1, so S1 does not bound positivity's reach; explore19's law does.
- N7 (RH verified to 3·10¹², Platt–Trudgian): the Riemann–Siegel main sum at T has n ≤ √(T/2π) = 6.9·10⁵ = X*(T).
  Verification by Riemann–Siegel already works at the information floor. explore30's off-line lattice (step 0.467 at
  3·10¹², invisible for X < 6.9·10⁵) is the lattice case of S1.
- No new bound on any node quantity.

## 4. Adversarial test (iv)
Davenport–Heilbronn (conductor 5, ExploreHelio's 304 on-line zeros and 9 real off-line pairs below 420), each off-line
pair against the on-line alternative (two zeros one spacing apart plus moved neighbours, M = 40), next to ζ with the
nearest pair moved off the line by the same δ (logs/dh.txt):

| | ρ at 0.6 L* | 0.8 L* | 1.0 L* | 1.1 L* |
|---|---|---|---|---|
| DH, 9 pairs, L*_DH = ½log(5t/2π) | 7.5·10⁻⁶ – 1.2·10⁻³ | 3.5·10⁻³ – 0.11 | 0.38 – 0.67 | 0.65 – 0.99 |
| ζ moved pairs, same t and δ, L* = ½log(t/2π) | 4.4·10⁻⁶ – 4.4·10⁻⁵ | 5.6·10⁻³ – 0.12 | 0.45 – 0.76 | 0.56 – 0.88 |

DH's first pair (85.699, δ = 0.3085) is hidden from its own coefficient data below X*_DH = √(5t/2π) = 8.3. DH's Weil
positivity fails only at c = 30.7. The law is the same for ζ and DH relative to their densities (gamma factor and
conductor), so the statement does not tell ζ from DH: no RH information. S4 tells RH-true from RH-false only through data
at n ≥ p, by construction.

## 5. Where it breaks
1. No noise. η > 0 is the margin of a robust argument: the two-point lemma is a robustness implication (D2), and the exact
   prime data of Z₀ and Z₁ differ. Margin-free indistinguishability needs infinitely many moved zeros (D3) and is not proved;
   the sharp threshold at L* is a numerically supported conjecture.
2. Above L* only local compensation was tried: ≤ 800 neighbours, moves ≤ 1 spacing, on-line. Far or off-line
   compensators were not tried, so "visible above L*" holds only for this family.
3. Large δ: with bounded moves the floor is the δ·log X ≲ 1.5 zone near the line, not the strip. Far-off pairs are
   excluded by short data, as classical zero detection says.
4. Fano/Assouad add nothing beyond S3: the data are deterministic.

## 6. Novelty
- arXiv math.NT listings 2026-10 (854 titles) and 2026-11 (280): keyword scan (explicit formula, Weil, minimax, Le Cam,
  detect, Turing, information, off-line, Beurling, Davenport, Selberg class, pair correlation, zero-free, zero density).
  Nothing on information floors for prime data. Nearest titles: "On Landau's local method for obtaining zero-free
  regions" and "The uncertainty operator, the Weil ground state and the zeros of zeta", both 2026-11, neither on this.
- web_search (arXiv) was throttled (HTTP 429, retried once).
- S4 is folklore in type: ζ times a Dirichlet polynomial with a functional equation, excluded by the Selberg-class axiom
  θ < ½.
- explore30 has the lattice case; explore19 has the positivity reach.
- New here: the Le Cam distance as a Weil-window operator norm, the √(t/2π) threshold with ζ's actual zeros at
  10³–10¹¹, the δ·log X ≲ 1.5 zone, the block (Assouad) test, the DH comparison and the 20× gap to positivity.
- No discovery claimed: a barrier on information, with no node quantity moved.

## Suggested NOGO line
- explore36 (Orca, exploration: minimax lower bounds — Le Cam two-point, Assouad hypercubes, Fano → N1, N7) BARRIER: primes ≤ X fix the zeros only through G(y) = Σ_ρ e^{iγ_ρ y} on |y| < log X (Weil windows c ≤ X); merging a zero pair at height t off the line (δ ≤ 0.2) and moving ≤ 100 neighbours by < 1 spacing changes every such window by ~10⁻⁶ of the bare pair's effect at log X ≤ 0.6·½log(t/2π) (falling like M^{−1.5}), never below 0.12 at log X ≥ ½log(t/2π) (t = 10³ … 10¹¹), so no positive-margin argument from primes ≤ X excludes near-line off-line zeros above 2πX², the Riemann–Siegel length (bounded moves also need δ log X ≲ 1.5; 20 % of a block hides at 0.5·L*); exact global pair ζ(s)(1 + bp^{−s} + p^{1−2s}): Euler product, ζ's gamma factor, conductor p², identical data below p, zeros at Re s = 0.9026 (p = b = 7); Weil positivity reaches only ≈ 13.9c (explore19) against 2πc² here; density-only, Davenport–Heilbronn obeys it with conductor 5: no RH information [ExploreMinimax](Riemann/Deep/ExploreMinimax/REPORT.md)

RH relevance: N1 (window records see about 1/20 of the height their data could cover); N7 (zero verification needs prime data up to √(t/2π), where Riemann–Siegel already sits); no RH information.

## Review
One `reviewer` pass (verdict only, claims inline): C1, C2 correct; C4, C5 correct with full hypotheses stated; C3 wrong for
b = −(p + 1) (pole cancellation), now excluded; C6's "iff" overstated, now a conjecture; D6's "≥ K on some vertex, K/2 on
average" wrong, now worst case ≥ K/2 and Hamming error ≥ K/2. All fixed in DERIVATION.md and here.

## Files and runs
| file | what | run (single-threaded, perl alarm) |
|---|---|---|
| lecam1.py | one pair: compensation (Levenberg–Marquardt, HS norm), operator norms; args set δ M-list fracs op [offset] | 8–259 s, 36–186 MB (M = 800) |
| lecam_stats.py | percentiles over pairs; args set npairs M fracs δ cap | 54–157 s, ≤ 44 MB |
| lecam_multi.py | blocks of K pairs (Assouad); args set B K-list fracs | 28 s, 39 MB |
| lecam_dh.py | DH off-line pairs vs on-line alternatives, ζ moved pairs at the same heights | 62 s, 46 MB |
| fpb_check.py | S4: zeros of A, zero side = prime side, two-bump sign | 4 s, 127 MB |
| logs/*.txt | outputs quoted above | – |

Run: `uv run --no-project --with numpy python3 lecam1.py z1e8 rel 100,200 0.6,0.8,1.0 op` from this directory. Inputs are
read from ../ExploreHelio/results and ../ExploreDPP/logs.

Process note: the first fpb_check.py run recomputed a 4001-node Gauss–Legendre rule per call. It ran 600 s to its alarm,
peaked at 287 MB (below the 300 MB cap) and printed nothing. It was fixed with precomputed nodes (4 s, 127 MB). Every
other script self-aborts above 280 MB.
