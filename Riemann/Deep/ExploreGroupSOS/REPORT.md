IDEA: certify Weil windows dually by a Krein positive-definite extension whose atoms sit at logarithms of S-units (S = primes < c), with the high-frequency envelope proved on the S-unit torus (group-ring SOS or branch-and-bound on 𝕋^{π(c)}); at c = 10 and 16 the node scan of a lattice (Chekmachev-type) certificate drops 75× and 244× below Chekmachev's own scans, with no more bad nodes than his matrix order (33 vs 33, 83 vs 133) and no precision gain (node N1).

# explore8: computer-assisted positivity certificates in group algebras (property (T) by sums of squares) -> RH tree

Orca exploration, 2026-11-05, 17:50–20:50 UTC. Remote field; bridge checked against N1 (Weil positivity on windows) and N5.
Notation: window f ∈ L²(−L/2, L/2), L = ln c (half-width a = L/2); Q_c the Weil form (= PROOF.md sec. 0, = Chekmachev's Q_c);
w_n = Λ(n)/√n, σ_c = Σ_{n<c} w_n; λ_min(c) = inf Q_c(f)/‖f‖² ≈ Θ_c = (2πΦ(a)/(2πc − 9/2))⁴ ≈ e^{−4πc} up to powers of c (Chekmachev's scale).

## 1. Survey (arXiv API and abs pages; queries in .scratch/explore8/arxq.py)
- Ozawa 1312.5431: Γ has (T) ⟺ Δ² − κΔ = Σ ξᵢ*ξᵢ in ℝ[Γ] (finite certificate). The engine is an order-unit lemma: Δ dominates every
  self-adjoint element of the augmentation ideal supported in a ball (Netzer–Thom 1110.5619: the SOS cone of the augmentation ideal has an
  interior point ⟺ H¹ = 0; Mizerka–Nowak 2301.07590: families of order units in powers of the augmentation ideal).
- Numerical SOS + exact rounding: Netzer–Thom 1411.2488 (SL₃(ℤ): Δ² − Δ/6, gap 1/72, double-precision SDP rounded to rationals),
  Fujiwara–Kabaya 1703.04555, Kaluba–Nowak 1703.09680, Kaluba–Nowak–Ozawa 1712.07167 (Aut(F₅)), Kaluba–Kielak–Nowak 1812.03456
  (Aut(F_n), n ≥ 6: symmetrization over S_n (Wedderburn blocks) and an induction Δ_n² = Sq + Adj + Op whose gap GROWS with n),
  Kaluba–Kielak 2306.12358 (Chevalley groups via root gradings), Mizerka 2404.10287 and Mizerka–Szymański 2504.16625 (inducing gaps for
  cohomological Laplacians), Nitsche 2009.05134 (strong SDP duality; the dual objects are harmonic cocycles, i.e. positive-definite
  functions; Aut(F₄) has (T)), Nitsche 2312.13917 (computer-free proof for large n).
- Positive-definite functions and extension: Krein (ℝ: every p.d. function on an interval extends), Rudin (fails on ℤ², ℝ²),
  Bakonyi–Timotin math/0512042 and 0808.1273 (free and amenable groups), Jorgensen–Pedersen–Tian 1507.02547, Burton–Juschenko 1912.12365.
- Weil side: Connes–Consani 2006.13771 (archimedean place; positivity for supp ⊂ [2^{−1/2}, 2^{1/2}]), Connes–Consani–Moscovici
  2310.18423 (semilocal prolate operator; semilocal Sonin spaces stable as S grows) and 2403.01247 (moment problem and Jacobi matrix of
  |∏_{v∈S} local factors|² on the critical line, S = {p, ∞}; q-series, integrality); Suzuki 2606.09096 (Weil form via Krein's screw
  function and de Branges spaces), 2206.03682, 2610.32704 (assumes RH). Chekmachev (Zenodo 23848597) periodizes on a circle of length
  exactly 2L ("extension by zero"); he uses no extension, no Fejér–Riesz and no SOS (pdf searched).
- Nothing found that applies group-ring SOS rounding, Krein extensions or torus Positivstellensätze to Weil windows.
- Novelty (light, for T2): arXiv listings math.NT, math.CA, math.OA, math.FA for 2026-10 and 2026-11 (title scan: Weil, Krein,
  positive definite, sums of squares, Fejér, explicit formula, Kazhdan, Positivstellensatz): only 2610.32704 (Suzuki, assumes RH) is on
  the Weil form; 2611.02967 (operator Fejér–Riesz on unions of intervals) is unrelated. Not a DISCOVERIES entry: a method idea, no record.

Mechanism map. SOS in ℝ[Γ] certifies that an ELEMENT is positive; Weil positivity asks that a FUNCTIONAL be positive definite. The two are
dual: the dual certificate of a window form is a positive-definite extension (Krein/Bochner), the dual of an SOS is a p.d. function
(Nitsche's cocycles). Three transfers follow, plus a remark.

## 2. Transfers, tests, verdicts

### T1. Ozawa rounding (numerical SOS, exact residual absorbed by an order unit) -> "certify with less precision". DEAD.
- Test (prec_test.py): even Galerkin matrix Q_N (closed forms of ChekAudit's galerkin_fast), Cholesky of Q_N − (λ_min/2)I at b bits,
  exact residual E in high precision, certificate iff ‖E‖_F < λ_min/2.

  | c | λ_min(Q_N) | fails at | certifies at | log₂(1/λ_min) |
  |---|---|---|---|---|
  | 5 (N = 60) | 1.12e-17 | 53 bits (double) | 64 bits | 56 |
  | 10 (N = 60) | 2.20e-43 | ≤ 128 bits | 160 bits | 142 |
  | 16 (N = 80) | 1.29e-74 | ≤ 224 bits | 256 bits | 245 |

- Lemma (trivial, decisive): a verifier that accepts every form within distance δ of its input accepts Q − δI, which is indefinite once
  δ > λ_min. So every certificate type (primal Cholesky, Chekmachev's lattice, dual extension, SOS) must resolve Q to ≈ λ_min ≈ Θ_c along
  the near-null direction: ≳ log₂(1/Θ_c) ≈ 17c bits (Chekmachev states "approximately 17c bits", pdf sec. 6). Ozawa rounding works for
  (T) because the gap is macroscopic (SL₃(ℤ): κ = 1/6 against ‖Δ‖ ≤ 24, relative gap ~10^{-2}); the Weil gap is e^{−4πc}: 10^{-17} at
  c = 5, 10^{-124} at c = 24.5, 10^{-242} at c = 47. The order unit of the Weil problem would be a form with the same near-null space as
  Q, i.e. Σ_γ |F(γ)|² over the zeros: no explicit "augmentation ideal" exists (R84's λ_min → 0 is intrinsic).

### T2. Dual certificate = Krein positive-definite extension; envelope by group-ring SOS on the S-unit torus. ALIVE (engineering, N1).
Identity (weil_lat.py; checked against galerkin_fast at c = 5 for P = 2L, 3L, 5L, agreement 5e-7 at 2e6 nodes): for every period
P ≥ 2L and every even distribution E on the circle ℝ/Pℤ supported in the gap L ≤ |x| ≤ P/2,
  Q_c(f) = (1/P) Σ_k ω_k |F(2πk/P)|²,  ω_k = m(r_k) + Ê(r_k),  m(r) = pole(r) + Re ψ(1/4 + ir/2) − log π − T(r) − 2Σ_{n<c} w_n cos(r log n),
with T(r) = −2Σ_j c^{−a_j}(a_j cos rL − r sin rL)/(a_j² + r²), a_j = 2j + 1/2 (Chekmachev's tail at r = πk/L). If all ω_k ≥ t then
Q_c ≥ t‖f‖². Chekmachev's identity is the case P = 2L, E = 0; it has negative ω_k (his certificate treats them with two node systems
and a finite matrix) and sporadic high-frequency negative ω_k wherever all cos(r log p) ≈ 1 (handled by a scan to k ≈ 2L e^{1+2σ_c}: 5.1e10 nodes at c = 47, 1.7e13 at 71).
- T2a, pure dual certificate (lp_ext.py, LP over E in a hat basis): c = 2, P = 6L, 60 hats, Σ|e| ≤ 100: every weight ≥ t* = 5.9e-5 on
  420 nodes (λ_min(2) = 1.33e-3). A matrix-free certificate exists at c = 2. With P = 2L the best E (one atom at x = L) gives t* < 0.
  At c ≥ 3 a double-precision LP cannot resolve t* ≤ λ_min (5.6e-8 at c = 3); a high-precision LP was not built.
- T2b, the envelope. Beyond the low band, ω_k ≥ ψ-part − |T| − |pole| − (max over the torus 𝕋^{π(c)} of the prime polynomial). With
  E = 0 that maximum is 2σ_c. Atoms of E at x = k log p (p^k ≥ c, wrapped mod P into the gap) add one-variable cosines; the best
  per-prime bound is a Fejér–Riesz/Fejér–Egerváry problem (β_p = w_p/cos(π/(K+2)) for one fixed term and K free ones), solved by LP.
  The lower limit for ANY extension is ‖C_p‖ on L²(window) (state h = f*f̃ restricted to Γ), computed as a Perron eigenvalue:
  4.901 at c = 24.5, equal to our Schur constant 4.903 (REPORT5): the weighted Schur test IS a binomial SOS in the crossed product
  L^∞ ⋊ Γ and is exact for positive operators (Collatz–Wielandt).

  | c | 2σ_c (E = 0) | β_sep, 8 free powers | ‖C_p‖_window | log₁₀ nodes: E = 0, P = 2L | β_sep | ‖C_p‖ |
  |---|---|---|---|---|---|---|
  | 10 | 7.08 | 3.95 | 2.51 | 4.2 | 3.7 | 2.2 |
  | 24.5 | 14.32 | 7.83 | 4.90 | 7.5 | 5.6 | 3.4 |
  | 47 | 21.61 | 12.04 | 7.30 | 10.7 (Chekmachev: 5.1e10) | 7.5 | 4.5 |
  | 71 | 27.70 | 16.12 | 9.21 | 13.4 (Chekmachev: 1.7e13) | 9.3 | 5.4 |

- T2c, actual lattice counts with wrapped per-prime atoms (scan_wrap.py; nodes = range that must be scanned, neg = nodes with ω_k < 0):

  | c | P = 2L, E = 0: nodes / neg / last neg r | P = 1.25·2L + atoms | P = 1.5·2L + atoms | P = 2·2L + atoms |
  |---|---|---|---|---|
  | 10 | 5466 / 50 / 5548 | 1375 / 41 / 1379 | 748 / 58 / 543 | 712 / 79 / 452 |
  | 16 | 116211 / 125 / 77404 | 12468 / 158 / 8608 | 5670 / 164 / 3907 | 4430 / 256 / 2322 |
  | 24.5 | 10665084 / 516 / 3955354 | 488255 / 731 / 302065 | 122580 / 823 / 57552 | 77011 / 1253 / 34074 |

  The scan shrinks 4–139×, growing with c (it is e^{2σ_c − β_c} times P/2L), but the bad nodes grow 0.8–2.4× (denser lattice; the atoms
  oscillate). In a Chekmachev certificate the bad nodes set the order of the high-precision matrix (7170 at c = 47, ≥ 84% of it sporadic
  high-frequency nodes), so the trade pays only where the scan dominates (c = 71: 1.7e13 GPU nodes).
- T2d, non-separable extension: atoms at S-units v = ∏ p^{a_p} (a in a small exponent box) whose log lands in the gap mod P; the envelope
  is the minimum of a trigonometric polynomial on the torus 𝕋^{π(c)} = dual of the S-unit group. Optimal envelope by LP with cutting
  planes on sampled torus points (torus_lp.py; sampled, NOT a certificate), then the exact lattice scan with the found E:

  | c, P | no extension: envelope / nodes / neg / last neg r | prime-power atoms | S-unit atoms (box) |
  |---|---|---|---|
  | 10, 2L (Chekmachev lattice) | 7.08 / 5466 / 50 / 5548 | – | – |
  | 10, 1.25·2L | 7.08 / 6832 / 50 / 4650 | 5.46 / 1375 / 41 / 1379 | 3.74 / 257 / 33 / 171 (33 atoms, box 3,2,1,1) |
  | 10, 1.5·2L | 7.08 / 8198 / 64 / 5548 | 4.72 / 794 / 59 / 543 | 3.38 / 220 / 31 / 109 (53 atoms) |
  | 16, 2L (Chekmachev lattice) | 9.94 / 116211 / 125 / 77404 | – | – |
  | 16, 1.25·2L | 9.94 / 145263 / 152 / 77404 | 7.72 / 15812 / 155 / 8608 | 5.41 / 1757 / 83 / 1070 (122 atoms, box 2,1,1,1,1,1) |

  Chekmachev's own certificates (ChekAudit cert.py logs): c = 10: scans of 14797 odd + 5444 even nodes, matrix order 33 (384 bits);
  c = 16: 313858 + 115462 nodes, order 133 (512 bits). Against these, the S-unit extension needs 268 nodes at c = 10 (envelope
  VERIFIED: branch-and-bound on 𝕋⁴ with a second-order Taylor bound, 4.6e6 boxes, Φ ≥ −3.7856 everywhere; bb.py, rescan.py) and
  1757 nodes at c = 16 (envelope sampled only), i.e. 75× and 244× fewer, while its negative weights (33, 83) do not exceed his
  matrix orders (33, 133). Coupling the primes (atoms such as 10, 14, 15, 3/2, …) is what helps: prime-power atoms alone leave the
  bad nodes unchanged. The gains grow with c because the scan length is ∝ P·e^{envelope} and the envelope floor is ‖C_p‖ (2.51 at
  c = 10, 3.64 at 16, 7.30 at 47) against 2σ (7.08, 9.94, 21.6): at c = 47 that floor would mean ~10^{4.5} instead of 5.1e10
  nodes, and the sporadic bad nodes beyond r ≈ 2πe^{8.3} (≥ 84% of Chekmachev's order-7170 matrix: 6047 selected odd nodes beyond
  k = ⌈5c²/3⌉) would disappear.
  Certifying the envelope is the open part at larger c: a group-ring SOS (Fejér–Riesz in π(c) variables; Dritschel: strictly
  positive trigonometric polynomials are SOS) or branch-and-bound on 𝕋^{π(c)}, which is cheap for π(c) = 4 but not for π(c) = 15
  (c = 47): couple the small primes, treat the large ones prime by prime. Tested: box-limited SOS on ℤ⁴ (torus_sos2.py, Gram support
  81 elements) gives 6.35 with S-unit atoms against 11.46 with prime powers only and 13.17 with none (degree too low to reach
  2σ = 7.08 even without atoms; larger Gram supports needed > 1 GB and were stopped).
- Verdict: a different certificate design for N1 (dual, Krein): with prime-power atoms it saves only in the cheap part (scan); with
  S-unit atoms it also cuts the negative weights by a third against the zero-slack lattice (50 → 33, 125 → 83; Chekmachev's matrix
  orders 33, 133), the expensive high-precision part, with the potential to remove the sporadic ones altogether at c ≥ 47. Precision is unchanged (T1). Next steps: (i) certify the torus envelope in interval
  arithmetic (done in float64 at c = 10) and at c = 16, 23 (π(c) = 6, 9); (ii) optimize E jointly for the low band and the envelope (one LP at
  ≈ log₂(1/Θ_c) bits); (iii) redo Chekmachev's c = 16 certificate on the P = 1.25·2L lattice with S-unit atoms and compare orders.

### T3. Semilocal positivity Q_S (archimedean place + all powers of the primes in S), "SOS on the subgroup generated by S". DEAD.
- Test (semi.py, semi_thr.py: even sector, Galerkin N = 50 at 80 digits; a negative value is a valid upper bound for λ_min):

  | S | first Galerkin failure (half-width, ±3e-5) | next prime p enters at (ln p)/2 | λ_min upper bound past it |
  |---|---|---|---|
  | {∞} | 0.4101 | 2: 0.3466 | −7.4e-3 at 0.45 |
  | {∞, 2} | 0.5720 | 3: 0.5493 | −1.3e-3 at 0.60 |
  | {∞, 2, 3} | 0.8125 | 5: 0.8047 | −0.12 at 0.90 |
  | {∞, 2, 3, 5} | 0.9771 | 7: 0.9730 | −3.6e-3 at 1.00 |

  Each Q_S dies 0.004–0.06 past the point where the first missing prime enters; the overshoot shrinks with λ_min. No semilocal positivity
  statement exists beyond the range where Q_S = Q. (Side datum: the archimedean form alone stays positive to 0.410 > (ln 2)/2, beyond
  Connes–Consani's window; even sector, numerics only.)

### T4 (remark). KKN-type induction and Wedderburn symmetrization. DEAD.
KKN's induction works because the (T) gap grows with n; the Weil gap shrinks with the window (λ_min ≈ e^{−4πc}), so gluing windows from
smaller ones puts the whole difficulty into the cross terms (Adj), which are the RH-hard part. The only symmetry of a window on ℝ is
x → −x (already used: two parity blocks); the Wedderburn decomposition on the idele class group gives the Dirichlet-character forms.

## 3. Adversarial tests
- Davenport–Heilbronn. T1: not a mechanism (pure numerics). T2: DH-blind, as a certificate method must be. By Krein, an extension
  exists iff the window form is positive, so on windows where the DH form is negative (its zero 0.8085 + 85.70i) no extension exists and
  nothing false can be certified. The per-prime version uses a structure DH lacks (−ζ'/ζ lives on prime powers, so the torus
  polynomial splits into one-variable Fejér–Riesz problems; DH's log-derivative lives on all n < c and does not split), but the S-unit
  version works for DH too (every n < c is an S-unit). Either way it buys certificate cost, never positivity: no RH lever.
  T3: needs an Euler product to truncate; DH has none, and ζ's truncations fail anyway.
- Precision (window c = 5, half-width 0.805):

  | certificate | working precision | size | margin |
  |---|---|---|---|
  | ours, Legendre-Galerkin (WeilRecord run2/cert_5.log) | 320 bits | N = 140 (2 × 70), Ω = 48, Schur 8000 cells | 2.8e-21 (τ = 1e-20) |
  | Chekmachev (our cert.py, ChekAudit) | 256 bits (80-bit scan) | order-9 matrix + 201 nodes, 0.5 s | ε₅ = 7.4e-18 |
  | Ozawa rounding of Galerkin Q_60 (T1) | 64 bits; 53 fails | 61 × 61, Galerkin space only (no tail bound) | 5.6e-18 − 1.4e-18 |
  | dual extension (T2) | ≥ log₂(‖terms‖/t*) ≈ 60–70 bits predicted; not run (double-precision LP cannot) | M atoms + ~100 nodes, no matrix | ≤ λ_min ≈ 9.2e-18 |
  | (c = 2, run) dual extension | ≈ 21 bits (Σ|e| = 100, t* = 5.9e-5) | 60 hats, 420 nodes | 5.9e-5 vs λ_min 1.33e-3 |
  | (at 1.6, for scale) ours REPORT5 | 480 bits | N = 4700 | 3.5e-124 |

  No proposed certificate beats log₂(1/λ_min) bits; T2 changes size (scan), not precision.

## 4. Ranking (RH relevance × chance × gain)
1. T2 (Krein extension with S-unit atoms, envelope by torus SOS): relevance N1 records only (window records do not approach RH);
   chance of a working certified version: moderate (the envelope certificate in π(c) variables is the risk); gain: scans 75× and 244× shorter
   (measured, c = 10, 16; ~10^6× at c = 47 if the envelope reaches its floor ‖C_p‖) and no more bad nodes, i.e. cheaper certificates beyond c = 71 (half-width 2.13). Score: low × moderate × moderate.
2. T3, T1, T4: dead.
N5 (Nyman–Beurling): no transfer found (its Gram matrix is a minimization, not a positivity certificate; d_N → 0 has no SOS form).

## 5. Suggested NOGO line
- explore8 (Orca, exploration → N1): property-(T) SOS methods: Ozawa rounding needs a macroscopic gap, Weil's is ≈ e^{−4πc}, so every
  window certificate needs ≳ log₂(1/λ_min) ≈ 17c bits; semilocal forms Q_S die just after the next prime enters (0.410, 0.572, 0.812,
  0.977); KKN induction fails (the gap shrinks with the window); prime-power-only Krein extensions shorten lattice scans but not their
  bad nodes (S-unit extensions alive, see report) [ExploreGroupSOS](Riemann/Deep/ExploreGroupSOS/REPORT.md)

## 6. Files
Riemann/Deep/ExploreGroupSOS/: REPORT.md; results.log (all printed numbers); scripts weil_lat.py (lattice symbol), test_identity.py,
lp_ext.py (T2a), env.py (T2b), scan_wrap.py and scan_ext.py (T2c), torus_lp.py, bb.py, rescan.py and torus_sos2.py (T2d), prec_test.py (T1), semi.py and semi_thr.py
(T3), lam_small.py. Run: uv run --no-project --with mpmath --with numpy --with scipy --with python-flint python3 <script> (torus_sos2:
--with cvxpy --with clarabel); semi.py, prec_test.py and lam_small.py import ../ChekAudit/evidence/galerkin_fast.py.
torus_lp.py writes tlp_*.npz (E found) into .scratch/explore8/ for bb.py and rescan.py. State: .scratch/explore8/STATE.md.
Resource note: one process at a time, single-threaded; three dense LP/SDP runs exceeded the 1 GB memory rule (1.6–2.0 GB, ≤ 15 min
each: torus_sos2 box (2,2,1,1), torus_lp box (4,3,2,2) at c = 10 and the c = 16 S-unit run); all others ≤ 0.9 GB. Side check: Chekmachev's Thm 1.1(b) (ε_c ≥ 0.97 Θ_c) is stated for c ≥ 3 only; our Galerkin bound
λ_min(2) ≤ 1.3296e-3 < 0.97 Θ₂ is therefore no conflict.

RH relevance: N1 (certificate engineering only; no RH mechanism)
