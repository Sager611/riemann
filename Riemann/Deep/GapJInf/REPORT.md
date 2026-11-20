NO RESULT: no certified ceiling (the per-gap class's dual objects are jet-balanced mixtures of single-gap pieces, not functions; periodic-zero-set witnesses miss by ≥ 0.04 in h, and an all-order witness needs ν's density); in floats an exact all-order per-gap optimiser gives κ_J = 3.341, 3.359, 3.367, 3.371, 3.373, 3.373 (J = 4..9, record amplifier; 3.374 after one amplifier step), increments halving then vanishing, so the J → ∞ limit is ≈ 3.375, not far above 3.358 (node N8)

# gapjinf: the J → ∞ limit of the per-gap method (node N8, Bui–Milinovich large gaps)

## Statement (exact extremal problem)
Units: gaps rescaled to length ≤ π; ν_κ = BM's main-term spectral measure ν (variable ξ = v − S) dilated by κ, so
m_j(ν_κ) := ∫ξ^{2j}dν_κ = κ^{2j}c_j(v). Only even moments enter, so ν_κ may be replaced by its symmetrisation.
- **Order-J class 𝒜_J.** (a, B): a ∈ ℝ^{J+1}, B a quadratic form on jets (f', …, f^{(J−1)}), with
  Q(f) = Σ_{j≤J} a_j‖f^{(j)}‖²_{L²(0,g)} + B(jet f(g)) − B(jet f(0)) ≥ 0 ∀ g ∈ (0, π], f ∈ H^J(0, g), f(0) = f(g) = 0.
  Success at κ: Σ a_j m_j(ν_κ) < 0 (then Λ > κ). Reflection f ↦ f(g − ·) kills the p + q even entries of B, and
  asymmetric end forms B1(jet(0)) − B2(jet(g)) with B2 ⪰ B1 are dominated by the telescoping B2 form: (V_J) of
  ExploreShape (φ_pq, p odd, q even) is the whole quadratic per-gap class. Symbol: P_a(x) = Σ a_j x^j (x = ξ²).
- **Dual cone.** W_J = closed cone of w(g, f) = ((‖f^{(j)}‖²)_{j≤J}, ([f^{(p)}f^{(q)}]_0^g)_{p+q odd}); 𝒜_J = W_J*.
  C_J := {m : (m, 0) ∈ W_J} (mixtures of single-gap pieces whose odd jet products balance: E[jet⊗jet] at entry =
  at exit). By the bipolar theorem: **order J succeeds at κ ⟺ m_{≤J}(ν_κ) ∉ C_J.** proj C_{J+1} ⊆ C_J, so κ_J ↑.
- **J → ∞ problem.** κ_∞ = sup_J κ_J = sup{κ : ∃J, m_{≤J}(ν_κ) ∉ C_J} (amplifier and v re-optimised at each κ).
- **Lemma A (balanced mixtures are Stieltjes).** m ∈ C_J ⇒ [m_{i+k}] (2k ≤ J) and [m_{i+k+1}] (2k+1 ≤ J) are PSD.
  Proof: ‖f^{(i+k)}‖² = (−1)^{k−i}⟨f^{(2i)}, f^{(2k)}⟩ + Σ[f^{(p)}f^{(q)}]_0^g, p + q = 2(i+k) − 1 odd, p, q ≤ J − 1;
  the brackets cancel in a balanced mixture, so Σα_iα_k m_{i+k} = E‖Σ(−1)^iα_i f^{(2i)}‖² ≥ 0. So C_∞ is a cone of
  positive measures σ on x ≥ 0: the J → ∞ method asks whether ν_κ is a "balanced-mixture spectral measure".
- **Lemma B (one gap length: explicit validity).** At the moment-cone level, balanced mixtures of pieces of one length
  π/ω (ω ≥ 1) are the spectral measures of F = sin(ωx)H(x) (converse: reflection doubling equalises the entry/exit
  jet Gram matrices; a stationary Gaussian chain of pieces then matches jets exactly; Bloch decomposition). Bloch fibres: on each lattice φ + 2ωℤ, F's content is d_n = v_n − v_{n+1}
  (v = H's content on the lattice φ − ω + 2ωℤ), so the realisable fibre weights form
  Cone_0 = cone{d⊙d̄ : Σd_n = 0} = {w ≥ 0 : √w_n ≤ Σ_{k≠n}√w_k ∀n} (dual: diag(c) ⪰ 0 on 1^⊥).
  Hence every valid symbol satisfies, for all φ ∈ (−1, 1), ω ≥ 1:
  P ≥ 0 on [1, ∞),   −P(φ²) ≤ [Σ_{k≠0} 1/P((φ + 2kω)²)]^{−1}.                                           (F)
  Check: Wirtinger P = x − 1 gives Σ_{k≠0}1/((φ+2k)²−1) = 1/(1−φ²): **equality for every φ** (ω = 1).
  So Wirtinger is extremal fibre by fibre; going below it at a low frequency φ requires P large at the partners
  φ ± 2ω, |φ ± 2ω| ≥ 2 − |φ|.
- **Lemma C (periodic zero patterns: genuine dual witnesses).** Z = {z_1..z_m} + Pℤ with all gaps ≤ π, θ ∈ ℝ,
  d finitely supported with Σ_n d_n e^{2πinz_i/P} = 0 (i ≤ m): F(x) = Σ_n d_n e^{i(θ+2πn/P)x} vanishes on Z
  (each Bloch fibre must, as Z is P-periodic), Re F, Im F are real with zero gaps ≤ π, and
  m_j(F) = Σ_n |d_n|²(θ + 2πn/P)^{2j} (almost-periodic Parseval). So for every valid (a, B) of any order:
  Σ_n |d_n|² P_a((θ+2πn/P)²) ≥ 0. Hence the dictionary cone 𝒟 = cone{Σ_n|d_n|²δ_{θ+2πn/P}} lies inside the
  closure of genuine spectral measures, and **m_{≤J}(ν_κ) ∈ proj_J 𝒟 ⇒ no inequality Σ a_j m_j ≥ 0 valid for all
  functions with zero gaps ≤ π (per-gap or not) certifies κ**. No gluing is needed (cone = sums of atoms).
- Mixed patterns matter: gaps (π, π, π/2) (P = 5π/2, m = 3): Π sin(π(x − z_i)/P) has spectral atoms at 0.4 and 1.2,
  weight ratio 0.382 (`logs/pattern_pi_pi_halfpi.log`): low-frequency content with a partner at distance 0.8, not 2.

## Method
- **Moments** (`mom.py`, `ampmom.py`): ExploreShape's central_forms with a vectorised Gm (agrees with joint.central_forms to
  3·10⁻¹³, d = 12, K = 4), d = 80, K = 24, κ ∈ {3.358, 3.40, 3.45}; amplifier = explore53's J = 5 record amplifier
  (`../ExploreShape/logs/amp_d80_k3.358_J5.json`, v = 1.361505). Check: its rational J = 5 form gives H/c0 = 7.5895·10⁻⁴
  (explore53: 7.5895·10⁻⁴). Uncentred moments fail the Stieltjes Hankel test at m_10 (m_12 < 0: cancellation in
  Σ binom·(X − v0)^i·(τ-moments), `logs/hankchk_k3.358.log`). `moments_c.py` centres the τ-moments ((τ1+τ2−1.2)^m by a
  shifted recursion, weight (X − v0 + 1.2)^i): Hankel tests then pass to m_14 with smoothly decaying eigenvalue ratios
  and the order-14 Hamburger Hankel is PD (`logs/hankchk_c_k3.358.log`); m_0..m_5 agree with the uncentred values to
  10⁻⁸, m_7, m_8, m_9 differ by 2·10⁻⁵, 4·10⁻⁵, 1.1·10⁻³ (uncentred ones wrong), `logs/ampmom_c_k3.358.log`.
- **Exact per-gap forms of any order** (`pgq.py`, `pgschur.py`): u-representation. With b_k the sine coefficients,
  α_l, β_l = f^{(2l)}(0), f^{(2l)}(g), e_k^{(l)} = (2/g)(α_l − (−1)^kβ_l), ε_{k,l} = (−1)^l e_k^{(l)} ω_k^{−2l−1},
  u_{k,n} = b_k − Σ_{l≤n} ε_{k,l}:
  ‖f^{(j)}‖² = (g/2)Σ_k ω_k^{2j}u_{k,⌊(j−1)/2⌋}² + [j odd] (β−α)_{(j−1)/2}²/g   (u_{k,−1} := b_k),
  [f^{(2m+1)}f^{(2l)}]_0^g = (β−α)_m(β−α)_l/g − (g/2)Σ_k (−1)^{m+l} ω_k^{2(m+l+1)} u_{k,m} ε_{k,l}.
  Each mode couples only to σ = α + β (k odd) or δ = β − α (k even); eliminating every b_k (k ≤ 3000, tail ∝ 1/k²)
  leaves two L×L matrices R_odd(g), R_even(g), L = ⌊(J−1)/2⌋: **valid ⟺ p(ω_k²) ≥ 0 (B_k = 0 where p(ω_k²) = 0)
  and R_odd, R_even ⪰ 0 ∀g ≤ π** (reviewer: singular modes need the generalised Schur complement; the
  infinite-dimensional completion is asserted, not proved here).
  Checks: all ‖f^{(j)}‖², [f^{(p)}f^{(q)}] of a degree-5 polynomial piece reproduced to 8 digits (`logs/pgq_test.log`);
  for J = 4, R = (2/g)·(S_odd, S_even + (a3+φ32)/2) of explore53's exact4 (Lemma 2), ratio 0.63662 = 2/π at g = π
  (`logs/pgschur_test.log`): an independent re-derivation of Lemma 2 extended to every J.
- **Optimiser** (`pgcut.py`, `runpg.py`): LP max −a0m0 s.t. Σ_{j≥1}a_jm_j = 1 and ⟨(a, φ), w(piece)⟩ ≥ 0 over cut pieces
  (w = (‖f^{(j)}‖², [f^{(p)}f^{(q)}])); separation = bottom eigenvector of R on 43 gap lengths, plus p ≥ 0 on [1, ∞).
  h_J > 1 ⟺ order J proves Λ > κ (explore53's normalisation). The LP dual is a balanced mixture of the cut pieces.
- **Genuine-function dictionary** (`dictlp.py`, `rundictlp.py`): cuts from Lemma C (Bloch fibres of 12–156 periodic
  zero patterns), an upper bound h_dict ≥ h_J.

## Evidence (logs in `logs/`)
- **Lemma B is far from the true class** (`logs/dictlp_sep_k3.358.log`, `.scratch/gapjinf/j3test.py` output below):
  single-gap-length atoms alone (12 values ω ∈ [1, 2]) allow h = 1.0868 at J = 3 (κ = 3.358, J5 amplifier), against the
  true h_3 = 0.97689. The sep-optimal J = 3 symbol (−0.8126, 1, −0.2173, 0.02987) is valid at each gap length with its
  own φ12 (valid ranges: g = π [0.061, 0.40]; 0.99π [0.064, 0.074]; 0.9π [0.070, 0.203]; 0.2π [0.196, 0.40]) but the
  ranges have empty intersection: the **common telescoping form across gap lengths is the binding constraint**, and its
  dual is mixing pieces of different lengths with only the total jet covariance balanced.
- **Genuine periodic patterns (Lemma C) do not close the gap**: h_dict(J = 3) = 1.0868 (12 one-gap patterns),
  1.0175 (36 patterns with 2–4 gaps), 1.0173 (156 patterns, 2–7 random gaps), against 0.97689
  (`logs/dictlp_mix_k3.358.log`, `logs/dictlp_mix2_k3.358.log`). Balanced-jet mixtures (the per-gap dual) are not
  approximated by functions with periodic zero sets: the ceiling of the per-gap class is not a "genuine function" ceiling.
- **Small gaps bind.** Cutting planes on gaps ≥ 0.05π only gave forms that are invalid at small g (degree-d polynomial
  pieces on tiny gaps need, e.g., a2 + φ12 ≥ 0 at J = 3, a_d + φ_{d−1,d} ≥ 0 in general): h = 0.99013, 1.00222, 1.00762,
  1.01031 (J = 4..7), each ≈ 0.0014 too high (`logs/pgS3_k3.358_J3-5.log`, `logs/pgS_k3.358_J6-7.log`,
  `logs/valform_k3.358.log`). Fixed by a second oracle (explore53's Legendre bubble basis, degree 30, gaps 10⁻⁴π..0.3π,
  basis function of degree d scaled by g^{min(d,J)−1/2}, which keeps every entry O(1) as g → 0; unscaled, the J = 7
  check returns a spurious −0.51 at g = 10⁻⁴π: this is the "float conditioning of 7th derivatives" that stopped
  explore53 at J = 7).

## Results (float): h_J of the whole order-J per-gap class, one amplifier
explore53's J = 5 record amplifier (r = 2, ϑ = 0.2499, midpoint mollifier, degree 80, v = 1.361505), its moments at
κ = 3.358 (`logs/pgU_k3.358.log`, forms `logs/form_k3.3580_J*.json`, checks `logs/valformU_k3.358.log`):

| J | 1 (Wirtinger) | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| h_J | 0.9730744 | = h_1 (Lemma 0) | 0.9746621 | 0.9884780 | 1.0008557 | 1.0063009 | 1.0090279 | 1.0103483 | 1.0104648 | ≤ 1.0111 (LP not converged) |
| h_J − h_{J−1} | | 0 | 0.0016 | 0.0138 | 0.0124 | 0.0054 | 0.0027 | 0.0013 | 0.00012 | ≤ 0.0007 |
| poly-basis check | exact | | −7.9·10⁻⁹ | −3.2·10⁻⁹ | −5.6·10⁻¹⁰ | −9.2·10⁻⁹ | −1.7·10⁻⁹ | −1.6·10⁻¹⁰ | −9.2·10⁻¹¹ | −2.4·10⁻⁷ |
| κ_J (h = 1) | (3.319)* | (3.319)* | (3.321)* | 3.3413 | 3.3592 | 3.3671 | 3.3711 | 3.3730 | 3.3732 | ≤ 3.3741 |

J ≤ 7 from uncentred moments (m_≤7 accurate to 2·10⁻⁵), J = 8, 9, 10 from centred moments (`logs/pgC_k3.358.log`, forms
`logs/formC_k3.3580_J*.json`; the uncentred J8, J9 runs gave 1.0103548, 1.0103551: the 0.1% error in m_9 had hidden the
J = 9 gain). κ_J: h at κ = 3.40 for J = 5, 7 (0.9720452, 0.9799693, `logs/pgU_k3.40.log`) gives dh/dκ = −0.686, −0.692;
κ_J = 3.358 + (h_J − 1)/0.69 (J = 5: 3.3592 = explore53's 3.3592 at degree 80). Same amplifier at every κ.
*Linear extrapolation over 0.04 with an amplifier tuned for J = 5; GapBar's Wirtinger ceiling with its own optimal
amplifier is 3.32658.

- The LP is a relaxation (finitely many cut pieces), so each h_J above is an upper bound for the order-J class up to the
  separation accuracy, and the validated form gives the matching lower bound up to the residual: h_J is bracketed to
  ~10⁻⁸ in floats. J = 5: our optimum 1.0008557 vs explore53's form 1.00083 (independent code: agreement 3·10⁻⁵).
- From J = 5 on the increments halve (0.0054, 0.0027, 0.0013), then J = 9 adds 0.00012 and J = 10 at most 0.0007.
- **Amplifier step** (`reopt.py`, form fixed, top eigenvector of (C0, Σ a_j C_2j(v)), v re-optimised): J = 8: 1.0103483 →
  1.0107946; J = 9: 1.0104648 → 1.0109566 (v = 1.360455), i.e. κ_9 ≥ 3.3739 (`logs/reopt_J8.log`, `logs/reopt_J9.log`).
- **J → ∞ (float, extrapolated): κ_∞ ≈ 3.375 (between 3.3739 and ≈ 3.377)** for this amplifier class: about +0.016 over
  the order-5 record (3.358 certified, 3.3592 float) and +0.007 over order 6. Trade-off for a record: J = 8 (three
  boundary jets f'', f⁗, f⁽⁶⁾, as J = 7) reaches ≈ 3.373–3.374, 95% of the remaining gain; J ≥ 9 adds < 0.001.
  explore53 owns the records; nothing here is certified.

## Checks
- Moments: fast Gm = ExploreShape joint.central_forms (3·10⁻¹³); explore53's rational J = 5 form on our moments
  H/c0 = 7.5895·10⁻⁴ (theirs 7.5895·10⁻⁴); Stieltjes Hankel test passes to m_9 (`logs/hankchk_k3.358.log`).
- Exact reduction: polynomial piece reproduced (`logs/pgq_test.log`); J = 4 R-matrices = (2/g)× explore53's exact4
  S_odd, S_even + (a3+φ32)/2 (`logs/pgschur_test.log`); our first (unfixed) J = 4 LP form fails exact4 by S_odd(π) =
  −5.6·10⁻⁴ (mode truncation, found and fixed: all modes ≤ 3000 eliminated exactly, `chk4.py`).
- Two independent optimisers agree: our J = 5 optimum for the J5 record amplifier, h5 = 1.0008557, against explore53's
  own J = 5 form (polynomial-basis SDP, different code) 1.00083 on the same moments; J = 3: h3 − h1 = 0.0016
  (explore53: ≈ 0.0015–0.0018); J = 2 = Wirtinger to 5·10⁻⁶ (Lemma 0).
- Every reported form passes an independent check on explore53's polynomial basis (degree 30, 700 gaps in
  [10⁻⁶π, π], degree-adapted scaling): min relative eigenvalue ≥ −10⁻⁸ (`logs/valformU_*.log`). Not a certificate.

## Review (reviewer subagent, verdict on the formulas pasted inline; `.scratch/gapjinf/review1.md` not written, verdict quoted)
- (1) u-representation CORRECT. (2) reduction to R_odd, R_even: GAP — singular modes (p(ω_k²) = 0, e.g. Wirtinger at
  k = 1, g = π) need B_k = 0 and the generalised Schur complement; the infinite-dimensional completion (mode-wise
  minimisation = infimum over H^J) is not proved. (3) Lemma A CORRECT (necessary conditions only). (4) Lemma B CORRECT at
  the moment-cone level; the converse needs reflection doubling + gluing (sketched above); the Wirtinger equality is an
  elementary resolvent identity. (5) Lemma C CORRECT. (6) LP bracket: GAP — sampled validity checks are not uniform
  bounds; an a0-repair needs a proved bound. **Overall: ROUTINE.** Accepted: all h_J here are floats, not certificates.

## Novelty
- arXiv listings (lst.py, math.NT/math.CA, this and last month), arXiv API (abs:zeros AND "spectral measure" AND gap;
  Wirtinger AND zeta AND gaps), Zenodo "Wirtinger": no source on per-gap inequalities of order ≥ 3, their J → ∞ limit, or
  spectral measures of functions with bounded zero gaps (`.scratch/gapjinf/nov_lst.md`, `nov_api.md`; hits: 2610.22927,
  2610.26038, 2611.14486, 0906.5458 — gap results by other methods). In repo: only explore53 (NOGO.md:694, TREE.md:530).
- Uncertainty-principle gap literature (orchestrator's pointer) does not decide the limit: Eremenko–Novikov (math/0301060)
  prove the reverse implication ("Suppose that Fourier transform of a function f is zero on the interval [-a,a]. We prove
  that the lower density of sign changes of f is at least a/pi ..."); Mitkovski–Poltoratski (0903.4499) and Poltoratski
  (0908.2079) determine G_X = sup{a : ∃ μ ≢ 0, supp μ ⊂ X, μ̂ = 0 on [0, a]} = 2πD_*(Λ) for separated Λ: gaps in the
  Fourier support of measures on a set, not mean values of functions vanishing on a set (`.scratch/gapjinf/lit_report.md`).
  Beurling–Malliavin only forces the spectrum of F beyond |ξ| = 1, weaker than Wirtinger.

## Caveat
- Nothing here is certified. h_J are float LP/eigenvalue computations; validity of the forms is sampled (gap grids,
  degree-30 polynomial basis), not proved; moments are float (reliable to m_9; J ≤ 9 only).
- One amplifier (explore53's J = 5 record amplifier at κ = 3.358), not re-optimised per J; the brief's "amplifier
  re-optimised at each κ" is not done beyond one eigenvector step (see below). Re-optimisation can only raise κ_J.
- No ceiling for the whole class: the dual of the per-gap class consists of balanced mixtures of single-gap pieces
  (Lemma A cone C_J); genuine functions with periodic zero sets (Lemma C) are far from it (h_dict(3) ≥ 1.017 against
  h_3 = 0.9747), and a balanced mixture matching ν's moments of every order needs ν's density, which no code in the
  repo computes (only moments). So "no per-gap inequality of any order certifies κ ≥ κ_cert" is not proved.
- κ_∞ is an extrapolation from J ≤ 9 (halving increments, J9 = J8). A late change at J ≥ 10 (new telescoping terms
  with p = 9) is not excluded.

RH relevance: N8 (large gaps between zeros, Λ → ∞). Quantifies how much the per-gap (local quadratic) step can still
give over the record: about +0.015 in κ, then nothing; a record beyond ≈ 3.38 needs non-local information (sign or
phase coupling between gaps, or other mean values than ∫|F^{(j)}|²).

## Files (Riemann/Deep/GapJInf/)
- `mom.py` (moments, vectorised Gm; CSH env = centring), `moments_c.py` (copy of ExploreShape/moments.py with centred
  τ-moments), `ampmom.py` (moments of one amplifier, explore53 check), `hankchk.py`.
- `pgq.py` (u-representation, Hurwitz tails, polynomial-piece test), `pgschur.py` (mode elimination, R_odd/R_even),
  `pgcut.py` + `runpg.py` (cutting-plane LP for the true order-J class; second oracle on small gaps), `valform.py`,
  `valform2.py` (independent polynomial-basis checks; `validate53.py` = copy of ExploreShape/validate.py), `chk4.py`,
  `chk4b.py`, `pgschur_test.py` (exact4 cross-checks), `reopt.py` (amplifier step), `extrap.py` (crossings, tails).
- `dictsdp.py`, `dictlp.py`, `rundict.py`, `rundictlp.py` (Lemma C dictionaries).
- Logs: `logs/pgU_k3.358.log`, `logs/pgC_k3.358.log`, `logs/pgU_k3.40.log`, `logs/valform*.log`, `logs/form*_k*_J*.json`
  (all forms, full precision, with their moments), `logs/dictlp_*.log`, `logs/mom*_*.log`, `logs/hankchk*.log`,
  `logs/ampmom*.log`, `logs/pgq_test.log`, `logs/pgschur_test.log`; superseded (invalid at small gaps):
  `logs/pg_k3.358_J1-5.log`, `logs/pgS*_k3.358_*.log`, `logs/valform_k3.358.log`. Peak memory ≤ 200 MB (memcap).

## Suggested NOGO line
- gapjinf (Orca, J → ∞ of explore53's per-gap method) NOT A BARRIER (uncertified): exact all-order per-gap optimiser
  (mode elimination + cutting planes, small-gap oracle, centred moments) gives κ_J ≈ 3.341, 3.359, 3.367, 3.371, 3.373,
  3.373 (J = 4..9, record amplifier), so the per-gap class tops out near 3.375; genuine-function (periodic zero set)
  witnesses are far from its dual, so no certified ceiling [GapJInf](Riemann/Deep/GapJInf/REPORT.md)
