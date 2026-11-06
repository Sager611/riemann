NO IDEA: on the critical line a function with a functional equation is a class-D system, so its only protected invariant is the window parity (Hardy/Selberg sign changes); the integer winding of a companion is the argument principle (its defect count equals the number of off-line zeros for every deformation parameter, and Levinson's method is its bulk estimate), integer protection needs a definite Krein metric, which is RH, and the one place the Euler product acts (a zero-free bulk symbol) is Gonek's finite-product theorem at finite cutoff and RH itself at infinite cutoff (nodes N8, N2, N7, N13, N1).

# explore18: topological phases of matter (tenfold way, Kitaev/SSH zero modes, bulk-boundary index theorems) -> RH tree

Orca exploration, 6 Nov 2026. Derivations: [derivations.md](derivations.md) (D1-D8). Raw numbers: t1_*, t2_*, t3_* _out.txt.

Short answer. The dictionary works and is exact, but every entry is a known theorem: Kitaev's Z2 = Hardy's sign change,
the Toeplitz/chiral winding = the argument principle with Levinson's companion, the Krein-definite (integer-protected)
case = de Branges/Hermite-Biehler = RH, the Euler product = a zero-free 'bulk symbol' = Gonek's finite Euler products.
All finite-window invariants treat Davenport-Heilbronn (DH) exactly like zeta; only density invariants separate them,
and those are Bohr-Landau (1914), blind to o(T) sets of off-line zeros.

## Survey (arXiv API, 34 queries; web_search arxiv, 6 queries; Nov 2026)
- No paper with mathematical content connects bulk-boundary correspondence, the tenfold way, Kitaev/SSH zero modes,
  Toeplitz/edge indices or quantum Hall indices with zeta zeros or RH (queries: topological insulator, bulk-boundary,
  Kitaev, edge states, index theorem, spectral flow, Berry phase, Majorana, Krein signature, PT-symmetric, Fredholm,
  winding number + critical line, long-range Kitaev + polylogarithm, chiral symmetry, quantum Hall, Altland-Zirnbauer +
  L-functions, zero modes, Hermite-Biehler, Pfaffian, gap labelling + Dirichlet, Maslov, Pontryagin, negative squares).
- Hits, all peripheral: Sierra 1404.4252, 1601.01797 (Dirac fermion in Rindler space with prime-located potentials;
  spectrum = zeros only under RH-type conditions); Bender-Brody-Mueller 1608.03679, Bellissard's comment 1704.02644,
  2607.19067 (seen; PT symmetry needs a positive metric); 2108.07852, 2503.09644 (seen; Majorana equation in Rindler
  space); 1607.06276 (a chiral 'Riemann lattice' after Dyson's quasicrystal; conduction only, explore5 territory);
  1005.4299 (Chern-Simons zero-mode holonomy giving zeta values; interpretive); 2408.13292 ('topological proof of RH',
  saddle point, no content); 2609.27962 (Morse function on a phase torus for elliptic-curve zero statistics);
  2611.06705 (3-manifold decision problems made unconditional by quasi-RH: a consequence, not a tool); 1910.14368
  (Connes-Consani: Maslov phases in xp); 2501.15436 (Witten index of non-Fredholm Toeplitz operators: half-integer
  winding = argument principle with half residues); math-ph/0506016 (Kellendonk-Schulz-Baldes: rotation numbers and gap
  labels); long-range Kitaev chains 2609.28735, 2009.04111 (polylogarithm symbols with zeta(alpha) at k = 0, real alpha).
- Classical inputs used (not re-derived): Altland-Zirnbauer 1997, Kitaev 0901.2686, Schnyder-Ryu-Furusaki-Ludwig
  0803.2786; Gohberg-Krein Toeplitz index; Avron-Seiler-Simon charge deficiency; Prodan-Schulz-Baldes 2016 (1510.08744)
  and Kellendonk-Richter-Schulz-Baldes edge index; Katz-Sarnak symmetry types; Krein-Langer Pontryagin spaces; de Branges;
  Levinson 1974; Levinson-Montgomery 1974; Lagarias 1999 (Re xi'/xi > 0 on Re s > 1/2 iff RH); Bohr-Jessen mean motion;
  Bohr-Landau 1914; Davenport-Heilbronn 1936; Gonek 0704.3448 (finite Euler products: RH up to finitely many exceptions);
  Bombieri-Hejhal 1995; Lamzouri 2610.26841.
- In the repo already, not repeated: spectral-flow counts (R308, R335), LeClair-Mussardo (R285), symmetric-rectangle
  windings (R077), Thouless/cocycles (R061, R297), K-theory regulators (R22, R330), KdV/Toda (R062), non-Fredholm
  Toeplitz symbols (explore13 T3), Pontryagin index (R315, R325), Hermite-Biehler (R29, R34, R51), Levinson and Speiser
  (R5, R14, R47, R72), Gram blocks (R59).

## Transfers, tests, verdicts
Test functions: xi; mv = xi with its 6th and 7th zeros (37.5862, 40.9187) moved to 1/2 +- 0.15 + 39.2524i (functional
equation and reality kept, no Euler product); DH = completed Davenport-Heilbronn function (functional equation,
off-line zeros 0.808517 + 85.699348i, 0.650830 + 114.163343i below 120).

### T1 (N8): Kitaev's Z2 invariant = Hardy's sign change - reproduces Hardy/Selberg
Class D (particle-hole symmetry): one level at zero is protected, two can gap out. Here the functional equation pairs
rho with 1 - conj rho at the same height, so a simple on-line zero leaves the line only by meeting another one (D1). The
window invariant nu(I) = sign A(a) sign A(b) = (-1)^{N_0(I)} is the Kitaev invariant; N_0 >= #odd windows.
Test ([t1_z2.py](t1_z2.py); zeta via a numpy Riemann-Siegel formula checked against mpmath to 1e-3 at t = 20, 1e-5 for
t >= 100): zeta on [10, 2000]: N_0 = 1517 = N(2000). Share of zeros certified by odd windows chosen without knowing the
zeros: Gram intervals 0.8945; equal windows of 0.25, 0.5, 1, 2 mean spacings: 0.9974, 0.9591, 0.6809, 0.2498.
DH on [1, 120]: N_0 = 64, two zeros with Re s > 1/2 (rectangle count), N = 68; Gram intervals 64/68 = 0.9412, equal
windows 0.9412, 0.9118, 0.6471, 0.1618: the same mechanism with the same yields, and the four off-line zeros leave no
trace in any window.
Adversarial: (a) DH: parity certifies its on-line zeros exactly as for zeta and proves nothing about the off-line ones;
(b) mv: gamma_6 and gamma_7 sit in adjacent Gram intervals [35.4672, 38.9992], [38.9992, 42.3636]; moving the pair
removes two odd windows, which only the bulk count N(T) exposes (Turing's method, how RH is checked to 3*10^12).
Known? Yes: Hardy 1914 (parity), Selberg 1942 (parity on windows of length c/log T with mollified moments).

### T2 (N8, N2): chiral winding of a companion = argument principle; Levinson is its bulk estimate - reproduces known
Class AIII: zero modes carry a chirality; edge modes >= |winding| (Toeplitz index). Transfer (D2): any E with
Re E = (positive) * A on the line gives N_0 >= |W| - 1, W = (1/pi) Delta arg E; the bulk side is W = N - 2 kappa + b,
kappa = zeros of E in Re s > 1/2 (defects), b a horizontal-edge term. For the exact companion E_c = F + cF' one has
E_c = A - icA' on the line, so no zero can cross the line as c varies, and c -> 0 gives
kappa(c) = #{zeros of F with Re s > 1/2} for every c > 0: the integer invariant is the argument principle.
Levinson's G = zeta + zeta'/L equals xi + (2/L) xi' up to O(1/L) (D2); his loss is in bounding kappa by Littlewood's
lemma (convex weights charge the good zeros just left of the line) and the mollified mean value, not in the index step.
Test ([t2_companion.py](t2_companion.py), argument principle on [1/2, 3] x window, adaptive mpmath):

| function, window | N_0 | zeros of F, Re s > 1/2 | kappa: c = 0.05, 0.2, 1, 5, Levinson E_L | W (same order) |
|---|---|---|---|---|
| xi, [10, 60] | 13 | 0 | 0, 0, 0, 0, 0 | 12.995, 12.981, 12.920, 12.912, 12.792 |
| mv, [10, 60] | 11 | 1 | 1, 1, 1, 1, 1 | 10.995, 10.981, 10.921, 10.913, 10.792 |
| DH, [83, 88] | 2 | 1 | 1, 1, 1, 1, 1 | 1.834, 1.549, 1.197, 1.043, 1.296 |
| DH, [112, 116.5] | 1 | 1 | 1, 1, 1, 1, 1 | 1.041, 1.099, 1.047, 1.010, 1.117 |

Adversarial: (a) DH and (b) mv: kappa = 1 per off-line pair for every c and for Levinson's companion, so the bound
N_0 >= N - 2 kappa is an identity and never 'proves' a false line zero; it proves a true one only when kappa is bounded
from bulk data, i.e. Levinson-Conrey (and that bound stops below 100% for zeta and for DH-type combinations alike).
Known? Yes: Levinson 1974; kappa(c) = 0 for all c iff RH, whose positivity form is Lagarias 1999 (Re xi'/xi > 0 on
Re s > 1/2); Levinson-Montgomery 1974 (c -> infinity, zeros of the derivative).

### T3 (N7, N8): the Euler product as a zero-free bulk symbol - Gonek's finite products; circular at infinite cutoff
Transfer (D4): models F = G B + G* B* (G a Gamma factor, B -> 1 as Re s -> infinity, F = 2 Re(G B) on the line) with
companion E = G B, so W = theta/pi - 2 kappa_B + O(osc arg B). A finite Euler product is invertible on Re s > 0
(kappa_B = 0, arg bounded): the line must carry the whole zeta-like count. Dirichlet polynomials and DH's sum of two
Euler products have zeros right of the line (Bohr-Jessen mean motion < 0), the bulk defects.
Test ([t3_bulk.py](t3_bulk.py), [t3b_where.py](t3b_where.py); t in [1, 250]; zeta has N(250) = 108):

| bulk symbol B | N_0 | W | kappa_B | off-line zeros |
|---|---|---|---|---|
| Euler P_X, X = 3, 7, 13, 29, 59, 101 | 108 each | 107.87-108.11 | 0 | none with 0.505 < Re s < 0.8 (X = 3, 29); 8, 2, 24, 100, 236, 396 near Re s = 0, 1 next to the poles of P_X |
| zeta partial sums, N = 2, 5, 10, 30 | 108, 92, 78, 102 | 107.5, 63.6, 42.2, 9.9 | 0, 22, 33, 49 | 2, 20, 40, 34 |
| DH partial sums, N = 5, 10, 30 | 171, 155, 132 | 171.0, 138.9, 118.5 | 0, 16, 26 | 0, 18, 42 |
| DH sum of two Euler products, X = 7, 13, 29, 59 | 157, 155, 157, 161 | 148.8, 150.7, 142.8, 122.8 | 11, 10, 14, 24 | 7 with 0.505 < Re s < 0.8 (X = 29) |

(DH-type count to 250: 171. X = 101 for the DH products not run: time.) Single Euler products put every zeta-like zero
on the line (N_0 = 108 = N(250), index tight to 0.13) and have no off-line zero with 0.505 < Re s < 0.8 (checked for
X = 3, 29); the DH combination loses 5-8 pairs from the line up to 250. The index of a non-exact companion is loose (N = 30: W = 9.9
against N_0 = 102): bulk zeros overcount the defects, Levinson's loss in miniature.
Adversarial: (a) DH: correctly flagged (its bulk symbol has zeros), so the mechanism is DH-failing, but only through
zero-freeness of B on Re s > 1/2, which for zeta (X = infinity) is RH, and a bounded arg on the line, which is false
for zeta (S(t) is unbounded); (b) mv: its bulk symbol zeta*Q vanishes at 0.65 + 39.25i (kappa_B = 1), flagged trivially.
Known? Yes: Gonek 0704.3448 (a related family built from finite Euler products satisfies RH up to finitely many
exceptions), Davenport-Heilbronn 1936 and
Bohr-Jessen (sums of Euler products have zeros right of the line), explore5 (off-line zeros of the DH model).

### Barrier notes (no numerics needed)
- N13, symmetry class (D6): a Hilbert-Polya H for {+-gamma} needs a chiral or particle-hole symmetry; xp with complex
  conjugation is class D. The 0D class-D invariant (Pfaffian sign) is the root-number parity: the forced central zero
  of Katz-Sarnak SO(odd); trivial for zeta. No tenfold-way symmetry forces a real spectrum; reality is self-adjointness
  or a positive PT metric (BBM, Bellissard), i.e. the definite Krein metric of D3, equivalent to RH. DH has the same
  symmetries.
- N1, gap closing (D7): protection needs a uniform gap; Weil's window form has lambda_min -> 0 even under RH (R84).
- Density invariants (D5): Prodan-Schulz-Baldes index pairings and gap labels are densities; their zero analogue (mean
  motion) is 0 for zeta on every Re s > 1/2 (Bohr-Landau) and nonzero for DH: they separate the two only at the level
  of zero-density estimates.
- N8, isolated zeros (D8): a zero alone in a symmetric disc is on the line and survives perturbations smaller than the
  gap; this is already the Bombieri-Hejhal/Lamzouri mechanism behind our 0.6724 and 0.3449 results.

## Ranking (RH relevance x chance x gain)

| transfer | node | relevance | chance | gain | verdict |
|---|---|---|---|---|---|
| T2 companion winding | N8, N2 | high: it is Levinson's frame | ~0: the index step is already lossless | none | known |
| T3 zero-free bulk | N7, N8 | high: shows where the Euler product must act | 0 at X = infinity (RH itself) | none | known (Gonek, DH) |
| T1 class-D parity | N8 | medium | 0 | none | Hardy, Selberg |
| Krein-definite protection | N13, N1 | high | 0 (equivalent to RH) | none | R315, R325 |
| density invariants | N7 | low (blind to o(T) sets) | 0 | none | Bohr-Landau |

## Suggested NOGO line
- explore18 (Orca, exploration: topological phases of matter, tenfold way, Kitaev/SSH zero modes, bulk-boundary and Toeplitz indices, Krein signature, gap labels → N8, N2, N7, N13, N1) NO IDEA: Xi on the line is class D, so the protected invariant is window parity (Hardy/Selberg; Gram windows certify 0.8945 of zeta's and 0.9412 of DH's zeros alike); the chiral winding of a companion F + cF' has defects kappa(c) = number of off-line zeros for every c (xi 0, moved pair 1, DH 1 per pair, also Levinson's companion), i.e. it is the argument principle and Levinson's method is its bulk estimate; integer protection needs a definite Krein metric (de Branges, ≡ RH); single Euler products are a zero-free bulk (edge count 108 = N(250) for every cutoff X ≤ 101) while DH's sum of two Euler products has bulk zeros, but at infinite cutoff zero-freeness is RH [ExploreTopo](Riemann/Deep/ExploreTopo/REPORT.md)

## Files
- [REPORT.md](REPORT.md), [derivations.md](derivations.md)
- [topo_common.py](topo_common.py): xi, mv, DH and derivatives; adaptive argument principle.
- [t1_z2.py](t1_z2.py) -> [t1_xi_out.txt](t1_xi_out.txt), [t1_dh_out.txt](t1_dh_out.txt): Z2 window parity.
- [t2_companion.py](t2_companion.py) -> t2_xi_out.txt, t2_mv_out.txt, t2_dh1_out.txt, t2_dh2_out.txt: companion winding.
- [t3_bulk.py](t3_bulk.py) -> t3_euler_out.txt, t3_zpart_out.txt, t3_dhpart_out.txt, t3_dheuler_out.txt; [t3b_where.py](t3b_where.py) -> t3b_where_out.txt.
- Run: uv run --no-project --with mpmath --with numpy python3 <script> <arg> (single-threaded, under perl alarm).

RH relevance: N8 (critical-line counting), with barrier notes for N2, N7, N13, N1.
