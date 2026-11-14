# B2Below / MECHANISM: how PR 137 goes below the contact B_2 (node N7)

Paths: `src/` = .scratch/b2below/src/Submissions/marked_anchor_induction/Qrh/ (PR 137 head 9d81a4c0). Notes appended as found.

## 1. Which inequality improves: the count window of the return, via the marked-moment slot capacity
- Old (native) gate of the marked inverse moment (and of the ordinary two-factor face): free width r, total slot width S,
  `r + 2S < 1 ∧ 2r + 8S < 3` (src/MarkedAnchor/Interface/Native.lean:51, :69; Contract.lean:96). At r ∈ [0.7335, 0.7775] the first binds: S < (1−r)/2.
- New (PR 137): `S ≤ markCap r − ν`, `markCap w = (1−w)/2 + g0 − g1·w`, (g0,g1) = (552429/10^7, 7227847/10^8), for
  r ∈ [anchorLeft, anchorTop] = [1467/2000, 311/400] (src/MarkedAnchor/Interface/FiberImproved.lean:47–55, :62–74 inverse, :78–87 ordinary at doubled width 2m).
  Gain over native G(w) = g0 − g1 w: +0.002227 at 0.7335, +0.000861 at 0.7524, 0 at g0/g1 = 0.76431, −0.000953 at 0.7775 (used only where it helps: the count is paid at the crossing).
- Consumer: ONE improved count window (src/MarkedAnchor/Return/ImprovedWindowQ.lean:12–20, :34–44; ImprovedScheme.lean:56–68), δ ∈ [189/512, 15/32], x(=p) ∈ [15/32, 1/2]:
  count exponent = max of inverse face 1 − δ(a + k·r) at r = cross, ordinary face 1 − δ(a + k·2(t−cross)), long piece 1 − δ + (5/6 − δ)(t − 1),
  a = x(1+2g0), k = 1 − x(1+2g1) ⇔ face exponent 1 − δ(r + 2x·markCap(r)). Native (g0=g1=0) = the audited R_d = 1 − δ(x + 2(1−x)m_d),
  m_d = (5/6 − δx)/(5/2 − δ − 2δx) (Hexagon2 model; checked: three-way balance cross = 2t/3 = 2m_d reproduces it exactly).
- So the improved inequality is the CLASS COUNT R_d (the number of zero-detected rows at amplitude δ), not E, G_N, the direct pairing (11/3) or the dual class sum.
  New count at the contact: value = N/J (ImprovedWindowQ.lean:37–40), N/J = 0.62630 vs R_2 = (281 − √1081)/396 = 0.626566 at (δ_2, 1/2).
  Linearized: ΔR = −6Vδx·G(c)/(3V + 2δ(1−x)) = −0.3159·G(c), c = crossing ≈ 0.7524, V = 5/6 − δ.
- Return: PR103 scalar system unchanged otherwise (PositiveGap/Return/Tuning/Data.lean: θ = 1 − c, c = 316633/2500000, ℓ = 4c − 1/3, y = h − 2ℓ, h = 81447/100000;
  Hexagon2 model: y = l_y = (1+b−ℓ)/2, h = (1+b+3ℓ)/2 ⇒ b = 2h − 1 − 3ℓ = 0.1091015 (PR103 contact: b_2 = 0.0999533, h_2 = 0.8096863)).

## 2. Which phases the proof keeps (inverse half, src/MarkedAnchor/Inverse/Poisson/*)
- Moment statement: Σ_{u ∈ F.rows} |inverse poly(u)·slot product(u)|² ≤ C U^{1+ε} (FiberImproved.lean:62–74): a mean square over the family rows u, invariant under u ↦ e^{iφ_u}·(row data). No cross-row (rows u) phase enters the statement.
- Proof: canonical energy over Gauss rows k, Σ'_k W(‖k‖²/Y)|Σ_U c_U χ_U(k)|² (Interface.lean, PairPoisson.lean). Second Poisson in k, EXACT (Grouping.lean "EXACTLY");
  column pairs U₁ = g∪a, U₂ = g∪b keep their phases: Phase.lean χ_A(d)·conj χ_A(h)·χ_A(g)^4·χ_A(f)^4 = conj χ_A(y); GaussRoot.lean: primitive Gauss root·conj(α)γ₂ = μ(a)·fixed ray character;
  result: Σ over physical rows y = f²(e g′)g′h, Ny ≤ Z^{H+κ}, H = 2(N′+V) − K′ = 1 − 6u (Interface.lean prefactor_identities, K′ = K+4u, N′ = N−u), of |c_g|²·QCP(a,b)·χ_a(−1)·conj P_a(y)·P_b(y).
- Phases dropped: triangle inequality over labels f, (g, E) and dual frequencies k (FDecomp.lean), and Cauchy–Schwarz ½(‖chA‖² + ‖chB‖²) on each physical row y (InBall.lean), so the child is again a mean square over physical rows y, paid by the anchor at radius H (induction).
- Kept: rows k (Gauss rows of the moment's OWN dual, exact Poisson) and the column (prime-product) phases up to the physical row; dropped: rows u (statement), physical rows y (child squares).
  Type: Heath-Brown-style recursive Poisson/duality bound for a selected sextic energy (a large-sieve-type mean value), not a signed sum over rows of a zero-detected class.
- BeyondB2 map: dual lever S1/S2 (signed class sum over rows u): untouched; direct lever S7 (theta-row pairing, constant 11/3): untouched (ℓ = 4c − 1/3, y = h − 2ℓ unchanged, Data.lean);
  the gain sits in the moment inputs that fix F_d, which BeyondB2 §1 filed under "positive by nature (moments, counts)" and treated as fixed.

## 3. Equivalent levers (levers.py, exact closed form of BeyondB2 B(σ); infimum.py on the audited Hexagon2 model)
| quantity | PR103 system | PR 137 line | PR 137 θ |
|---|---|---|---|
| infimum B* (b re-optimized) | 0.873381967 (= B_2 to 5e-10 grid) | 0.8733466805 (b = 0.11039, h = 0.81511) | 2183367/2500000 = 0.8733468 at h = 0.81447: B*(b) = 0.8733467682, margin 1.55e-7 |
| count at contact R_d | 0.6265692817 | 0.6262976996 (ΔR = −2.716e-4; crossing 0.7524124) | |
| σ_eq (B(σ) = B*) | 0 | 3.270e-4 | 3.259e-4 |
| λ_eq = B_2 − B* | 0 | 3.529e-5 | 3.517e-5 |
| R310 F_d multiplier | 1 | 1 + 7.27e-4 (× dB/dε = −0.0482 → −3.5e-5 ✓) | |
- Linear lever at the contact: ΔR = −6Vδx·G(c)/(3V+2δ(1−x)) = −0.31591·G(c); dB/dG(c) = −(3h/4)·0.31591/(3+4δ_2) = −0.0408. 0.873 would need G(c) ≈ 0.0094.
- Full-system entries (BeyondB2 fmodel list) at θ, h = 0.81447: active 1.55e-7, next 1/250 − E 2.68e-4, all others ≥ 1.1e-3: the contact is the only active margin.

## 4. Reach: the anchor induction converges; ceiling 0.8733459 (chain_limit.py, ceiling.py, perturb.py)
- The chain is an explicit Euler scheme. One stage from anchor (R, z): Y(r) = (r/R)z + (1/k − 2)(1 − r/R)/6, k ≥ thr (gate 1: G(k) = H·F0 > 0, F0 = 0 at k = thr; High/Algebra.lean:59, :74, :219).
  Limit STEP → 0, pads → 0: r z′ = z − T(r, z), T = (1/thr − 2)/6, z(7/9) = 0, i.e. z(r) = r ∫_r^{7/9} T(ρ, z(ρ)) ρ^{−2} dρ.
  With wb = R: thr = (5R/4 − 5/6)/(R − 1/2 − z), so T ≈ (7 − 9ρ)/(3(15ρ − 10)) at z = 0: T = 0 at ρ = 7/9 (thr = 1/2, start pinned), T = 0.059 at 0.7524.
- Why wb = R is optimal (reviewer GAP → CLOSED): High/Algebra.lean Stg.Hyp needs l(R) ≤ z_R and 0 ≤ a1 ≤ 2/5 (ChainCheck LineGates 0 < a1 < 2/5), so wb − R = (z − l(R))/(1/2 − a1) ≥ 0;
  d thr/d wb = (5/4)((1−R)/2 − z)/(3wb/2 − (1+R)/2 − z)² > 0; the zero-pad chord through (R, z) attains wb = R (Lean chain: wb − R ≈ 1e-6). Finite stages only approach it (strict gates): a supremum.
- Rigorous comparison (every admissible schedule's gain ≤ z_∞; reviewer-2 repair with w = (33ρ − 14)/48 − z_∞): REPORT.md §4.
- Why small steps win: k(R, z) at wb = R has ∂k/∂R = (5/24 − 5z/4)/(R − 1/2 − z)² > 0, so each stage pays with k at its right end; refinement re-evaluates k lower. One stage 0.7775 → 0.70 gives z(0.7524) = 1.6e-5.
- Convergence of the gain at the contact crossing c = 0.7524 (PR137 line there: 8.61e-4):
| schedule | pads | z(0.7524) | z(0.7335) |
|---|---|---|---|
| PR137: 176 × 1/4000 from 311/400 | PR137 (k + 1/5000 ↑10⁻⁶, 4e-8, 2e-7, ↓10⁻¹⁰, gate 1 3e-5) | 8.617e-4 | 3.0764e-3 (exact rationals: chain.py) |
| 1/40000 | PR137 | 8.350e-4 (pads accumulate) | 3.0353e-3 |
| 1/1000, 1/4000, 1/40000, 1/400000 | 0 | 8.436, 8.705, 8.786, 8.794 e-4 | 3.0342, —, 3.1105, 3.1122 e-3 |
| 1/400000 from 7/9 | 0 | 8.795e-4 | 3.1123e-3 |
- Endpoint with the count paid by the exact gain curve (balance solved per (δ,x), b re-optimized; full-system entries all ≥ 2.8e-4 except the contact):
| variant (freedom) | B* |
|---|---|
| PR137 line (Lean) | 0.8733466805 |
| PR137 chain curve instead of the line | 0.8733466253 (line costs 5.5e-11) |
| **continuum, zero pads, whole window (CEILING)** | **0.8733458652** |
| same, PR137's window δ ∈ [189/512, 15/32] (its left edge then binds) | 0.8733459793 |
| keep gate 1's 3e-5 / keep k pad 1/5000 / STEP 1/4000 / start 0.775 / start 7/9 | 0.8733459858 / 0.8733460769 / 0.8733462487 / 0.8733462229 / 0.8733458617 |
| [lemma constant] Poisson radius 1 − 5u (u = (1−H)/5) | 0.8733383642 |
| [lemma constant] F0's 5/6 → 5/6 + 0.001; 5/4 → 5/4 − 0.001 | 0.8733416818; 0.8733426899 |
- Binding constraint: gate 1 (HIGH-row long-branch margin, k ≥ thr) with wb ≥ R and the start at thr = 1/2 (R = 7/9). The other gates have slack ≥ 1e-4 at the contact crossing (gate 2, S < 137/1000, binds only near r = 0.7328, where also k → 9/25).
- Iterating does not help: the chain already feeds every earlier anchor to the HIGH child and the base anchor has z(7/9) = 0 (no gain above 7/9 since k < 1/2 ⇔ 1/k − 2 > 0); θ enters no gate.
- Secondary cap of the window shape: native_top 1 − δ(x + (1−x)·311/400) and plain_top (hplain1) stay ≥ 0.6217 at the contact, so even unlimited gain on [0.7335, 0.7775] keeps R ≥ 0.6217 (ΔR ≥ −0.0049) unless Top also moves.
- To reach 0.873 by this lever: G(c) ≈ 0.0094 (dB/dG = −0.0408), 10.7 × z_∞(c). The strongest lemma-constant change tested (radius 1 − 5u) gives only 1.2 × z_∞(c).

## 5. Falsification of the ceiling against the tracker (read 14 Nov 17:25 UTC; the session was paused 12 Nov 22:36 → 14 Nov 17:23)
- QuasiRiemannTracker PR 141 (Joseph Bebel, "Dirichlet nonvanishing at 0.8733458561067185768 via continuum marked anchors", opened 13 Nov 04:04, MERGED 13 Nov 08:28 UTC, head 03a1e5ad):
  "extends Sil's admitted marked-anchor induction with a continuum supplier, averaged floor and exact return"; QRH-Decimal 0.873345856106718576801184796.
- PR 143 (Bebel, OPEN, "Exact continuum contact", depends on PR 142): θ_∞ = the unique root in (13/15, 7/8) of 3888·19^5·(15t − 13)^4·(1 − t)^5 − 7·(7 − 6t)^4·(9t − 7)^5 = 0.
- Our ceiling, recomputed exactly (ode.py): zero pads and wb = R make the continuum limit the LINEAR ODE r z′ = z − T, T = (14 − 18r − 12z)/(30(3r − 2)), z(7/9) = 0
  (integrating factor r^{−4/5}(3r − 2)^{−1/5}: the source of the 5th powers in PR 143's polynomial). RK4 gain curve + count balance + golden search in δ (x = 1/2) and b:
  B* = 0.8733458561067, b = 0.11079928, δ = 0.42582999, crossing 0.75232203, count 0.6265272858; PR 143 root 0.8733458561067185768; difference 1.1e-16.
- So the ceiling is exact and is attained: PR 141 sits at a rational just above θ_∞, PR 143 proves θ_∞ itself. Our float ceiling.py value 0.8733458652 (12 Nov 22:26) was 9.1e-9 high (δ-grid + step convergence).
- No submission goes below θ_∞; the marked-anchor mechanism is exhausted at θ_∞ = B_2 − 3.611·10⁻⁵ unless a lemma constant changes (§4 table).

## 6. PR 141 = this ODE; PR 143's polynomial derived (REPORT.md §6; poly143.py, theta_inf.py)
- PR 141 Curve.lean: flow(r, z) = A z + b, A = (15r − 8)/(r(15r − 10)), b = (3r − 7/3)/(r(15r − 10)); ideal = P − (P(7/9)/q0(7/9))·q0, P = 11r/16 − 7/24, q0 = (r⁴(r − 2/3))^{1/5}: identical to ode.py's ODE and to z_∞ = (33r − 14)/48 − (35/144)(9r/7)^{4/5}(9r − 6)^{1/5}.
- Contact: δ = (45w − 14)/(6(9w + 1)) (b-coefficient 0, R_d = 5/9 + δ/6); balance δ = thr(w, ideal(w)) ⇒ ideal(w) = −(45w² − 2w − 24)/(2(45w − 14)); θ = (99w − 14)/(12(9w − 1)).
- (P − ideal)⁵ = (5/16)⁵·7·w⁴(w − 2/3) with P − ideal = (95/48)(27w² − 12w − 4)/(45w − 14) ⇒ −48600000(15θ − 13)·Q₉(θ) = 0, Q₉ = PR 143's polynomial (coefficients identical).
- θ_∞ = 0.873345856106718576801184796065, w* = 0.752322053814, δ* = 0.425830047445; λ_eq 3.61105e-5, σ_eq 3.34630e-4.
