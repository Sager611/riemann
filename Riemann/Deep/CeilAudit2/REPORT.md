ROBUST: B1 (θ_∞ = 0.8733458561 is the ceiling of PR 137/141's own lemmas: the HIGH line, thresholds and windows cannot enter the binding class; the barrier's excluded constant, the Gauss extension cost 4 behind "radius 1 − 6u", is ~40× more potent than B2Below stated, and a V14-type twist of it, outside the class, gives 0.8733350 in floats); B2 (B♭ is pinned in NF3's class, and the record that passed it changed one listed input only: BFlat's own map with M2 prime cost 2 is B_2 exactly, with δ_2 and b_2, to 40 digits); B3 (13/15 = direct line at the count-free generic-row edge, robust to every row-wise input up to GRH and to signed full-family sums, which equal their diagonal in a finite model); B4 (BeyondB2's signed large sieve: a Vaughan Type II split with any large-sieve constant ≥ U loses Z^0.089, worse than LS1's Z^0.033).

(ceilaudit2, Orca, 20 Nov 2026, 03:24–04:05 UTC)

# CeilAudit2 — falsification of our own N7 barriers (zero-free regions)

Scope: B1 b2below (marked-anchor ceiling θ_∞), B2 bflat (B♭ = 0.874886647), B3 the 13/15 floor (R109, QRHImprove REPORT3, R291),
B4 BeyondB2's signed-large-sieve loss (TREE N7 "phase levers σ, λ remain unproved"). Lessons 191/235/309: a barrier must cover every
freedom the method's own lemmas admit, name the step it fixes and its exact input class. All runs single-threaded under
bash .scratch/memcap.sh 290 <s>; footprints from /usr/bin/time -l in logs/ (memcap samples every 2–3 s and prints 0 for short runs).
Tracker read through raw.githubusercontent.com (PR 137 = submissions/2026-11-12-marked-anchor-induction; PR 141 = 2026-11-13-continuum-anchor).

## Notes (written as I go)
- B1 source read (raw.githubusercontent, PR 137 = submissions/2026-11-12-marked-anchor-induction, files Inverse/Poisson/Interface.lean,
  High/Algebra.lean). Stage data (R, z, k, a0, a1, H): r = HR, u = (1−H)/6, s = H·S_R + u/k, S_R = (1−R)/2 + z; new capacity S_r = s + u
  (chain.py = ChainData), so z_r = Hz + u(1/k + 1) − (1−H)/2 = Hz + (1/k − 2)u. The 6 in u: prefactor_identities, K' = K + 4u (Gauss rows
  after the fourth-power Gauss selection), N' = N − u, H = 2(N'+V) − K' = 1 − 6u; payment (K−F+5u) + (H+u) = F.
- Gate 1 = Stg.G_k: G(k) = H·F0, G(ν) = H + u − H/6 − (5/6 − ν)(3gν/2) − ν(r+s) (long branch, sixth-power restoration); at ν = k,
  g = r for EVERY line through (R, z) (g ν = r + u(1/k − 1/ν)/(1/2 − a1)), so the HIGH line (a0, a1 ∈ [0, 2/5]) cannot enter the binding
  class; u cancels in G(k) (k·u/k = u), so thr = (5R/4 − 5/6)/(R − 1/2 − z) is independent of the Poisson radius.
- Consequence for the radius freedom: with K' = K + (4−ρ)u (Gauss cost 4−ρ per unit slot width), H = 1 − (6−ρ)u, the payment identity
  still closes, and z_r = Hz + u(1/k − 2 + ρ/2), u = (1−H)/(6−ρ). B2Below's "radius 1 − 5u" row kept (1/k − 2) and so understated it.
- Nearest proved source of such a ρ: Liu V14 (Hexagon 2610.00655v3) Thm tw:thm:norm, origin A_max + 4W < 1 → A_max + (4−ρ)W < 1,
  ρ < ρ_max(p) = 6(p−6)(2p−3)/(p²−193p+198), paid by moving radicals in returned sums; B_2's chain does not use V14 (Hexagon2 §1).
- B1 float test (code/b1_twist.py, logs/b1_twist.log, peak 37 MB, 36 s): continuum ODE r z' = z − T_ρ, T_ρ = (1/thr − 2 + ρ/2)/(6 − ρ),
  start where T_ρ(r0, 0) = 0, count window + return as B2Below/ode.py. ρ = 0 reproduces θ_∞ (−5e−15); 'udiv' mode reproduces B2Below's
  radius-1−5u row (0.8733383543). Consistent bookkeeping: ρ = 0.0830 (V14's ρ_max(6κ)) 0.8733349830 (−1.09e−5); ρ = 0.1304 (ρ_max(2))
  0.8733279269 (−1.79e−5); ρ = 0.1452 (ρ_max(3)) 0.8733256007 (−2.03e−5); ρ = 1 (radius 1 − 5u) 0.8730268822 (−3.19e−4; B2Below: −7.5e−6).
  Start moves to r0 = 0.7818–0.7849 (gain ρu/2 > 0 even at the native anchor). So the radius constant is ~40× more sensitive than
  B2Below's table says; still outside PR 137's lemma class (a twisted child with a moving radical is not covered by ChildBound, which is
  uniform only over ray characters mod 4).
- B1 cross-check (code/b1_chain.py, logs/b1_chain.log, 34 MB, 94 s): zero-pad Euler chain (step 1/40000, k = thr, wb = R) with the
  twisted bookkeeping agrees with the ODE (B* +4e−8 above it at every ρ, the usual finite-step bias): ρ = 0 0.873345895 vs 0.873345856;
  ρ = 24/289 0.873335027 vs 0.873334983; ρ = 0.14516 0.873325649 vs 0.873325601. θ = 0.873 needs ρ = 1.042 (cost 2.96 instead of 4).
- Where the 4 lives (Select/Scale.lean selection_exponents, l.29 and l.323–339): K' = 2(N+V) − 1 + κ + 4u', because the selected slot product
  enters the Gauss frequency as a FOURTH power, k ↦ k·(∏p_i)⁴, ‖·‖² ≤ Z^{2(N+V)−1+κ}·T⁴, T⁴ = Z^{4u'}; Hc' = 1 − 6u' − κ + κP; payment
  K' − N' − V + Hc' = N + V − u' + κP. With 4 → 4 − ρ the payment is unchanged (the ρ terms cancel), child widths (Pmax/H, (Smax−u)/H) keep
  their form, so S_r = H S_R + u/k + u holds with H = 1 − (6−ρ)u. This is exactly V14's "extension into the Gauss frequency costs 4 per unit
  width" (Hexagon REPORT l.5), which V14 Thm tw:thm:norm relaxes to 4 − ρ for p = 6κ ∈ [112/25, 9/2], ρ_max ∈ [24/289, 5662/67351].
- B2 float test (code/b2_bflat.py, logs/b2_bflat.log, 8 s, 21 MB footprint): BFlat's own stage-map pipeline (bflat_sens.py definitions,
  executed read-only), ONLY the M2 prime cost 6κ replaced by a constant c: c = 9/2 0.8748869894 (= B♭ up to 3e−7), 4 0.8747451, 3 0.8743125,
  5/2 0.8739542, **c = 2: 0.8733819666132, δ_c = 0.426082 = B_2 = (777 − √1081)/852 and δ_2 = (61 − √1081)/66 to all printed digits**;
  7/4 0.8729322 (valid, rectmax −3e−9); 3/2 0.87224 but another class binds (rectmax +4.8e−4 at (0.291, 0.21)); 1 → 13/15 (invalid).
  At c = 2: M1 cost 2 → 19/10 gives 0.8732001; amplification 5/6 → 4/5 gives 0.8732298. θ_∞ ≈ M2 cost 1.98 (or M1 cost ≈ 1.98).
  So B♭ was passed by exactly ONE of BFlat's five listed needs (an M2 prime cost below 6κ: 2611.00082's second induction n1+n2+2z ≤ M,
  Hexagon2 §1), not by the amplification exponent its NOGO line names; BFlat's per-unit sensitivity (+0.00019) misranked it because the
  cost moved 2.5 units.
- B3 float test (code/b3_family.py, logs/b3_family.log, 19 MB, 8 s): BeyondB2's finite model (characters mod q, exact L by FFT), FULL family
  (R = 1, g = 0: the generic-row end), T = L(1/2 − θ, χ)/L(0.75, ηχ̄). Signed sum / diagonal (#rows)/L(1.25 − θ, η) = 0.987–1.025 for
  q = 12289 … 163841, θ ∈ {0, 0.047}, η ∈ {χ_7, 1}: the family first moment IS its diagonal; r_full = |ΣT|/Σ|T| = 0.11–0.18 (η = 1, flat
  in q: no power saving); after removing the diagonal r_off = 1.8e−4 … 6e−3, r_off·√q = 0.07–0.68 (square-root cancellation).
  So the signed-sum freedom at the generic-row end reduces C0's sum to its family diagonal, which has the same exponent as the absolute
  bound; it moves the floor only if the actual (14.7) summand's diagonal vanishes (e.g. by the detector's probe subtraction) — not
  decidable from our files; and the direct side's reflected energy (5ℓ − 1)₊/4 forces 13/15 independently (REPORT3), so both would be needed.
- B4 test (code/b4_typeII.py, logs/b4_typeII.log, 19 MB): BeyondB2's "LS1 loses Z^0.0333". Freedom: Vaughan/Heath-Brown split of μ
  into Type II pieces, Cauchy–Schwarz and a family mean value with ANY large-sieve constant Δ(N) ≥ max(U, N) (even the conjectural diagonal
  sextic one). Exponents in Z at the contact (δ_2, a = (1+δ_2)/2 = 0.713041, h = 0.809686, R = 0.626569): absolute a + hR = 1.220366;
  LS1 1 + hR/2 = 1.253662 (loss 0.0332965, = BeyondB2 C1); Type II ≥ 1/2 + h = 1.309686 (loss 0.0893: needs h(1 − R) < δ/2, here
  0.302 > 0.213); subset large sieve (Δ = N + #𝓑·U^{1/2}) 1.412168. On the whole class domain (native count) Type II never beats the
  absolute route (min loss 0.0085 at δ = 1/8, x = 0). ROBUST: the family size U = Z^h > Z^{1/2} caps every mean-value signed tool.

## B1 — marked-anchor ceiling θ_∞ = 0.873345856106718576801 (b2below, NOGO l.661) — ROBUST in its class
**Step fixed.** Gate 1 of PR 137's HIGH child in the continuum limit: G(k) = H·F0 ≥ 0, i.e. k ≥ thr = (5R/4 − 5/6)/(R − 1/2 − z) at
wb = R, with the stage gain Y = Hz + (1/k − 2)u, u = (1 − H)/6, started at thr = 1/2 (R = 7/9); plus the return's zero-b contact
(B2Below §4, §7).
**Input class.** PR 137's stage lemma (High/Algebra.lean Stg.Hyp; Select/Scale.lean selection_exponents with Gauss extension 4u′ and
physical radius Hc′ = 1 − 6u′; ChildBound uniform over ray characters mod 4), PR 141's continuum supplier/averaged floor/exact return,
the count envelope 𝒞 of §7; any finite chain, pads ≥ 0, lines, h, window.
**Freedoms the lemmas admit.**
| freedom | verdict | reason |
|---|---|---|
| stages, pads, threshold, start, window, h, count curve | in class | B2Below §3 table, §7 obstruction |
| HIGH line slope a1 ∈ [0, 2/5] (Stg.Hyp) | adds nothing | g(ν) = r + u(1/k − 1/ν)/(1/2 − a1), so g(k) = r for every line through (R, z); G(k) = H·F0 with u cancelling (reviewer A) |
| several thresholds / cutoffs | adds nothing | the new capacity is one uniform bound (worst LOW amplitude k); witnesses below R are not yet anchors, and the continuum limit already sits at width r |
| phases kept by Phase.lean/GaussRoot.lean | outside | FDecomp (triangle) and InBall (½(‖A‖² + ‖B‖²) per physical row) drop them; a bilinear a ≠ b child is not a mean square, so no anchor of the class consumes it |
| amplification power (5/6, 5/4 = (3/2)(5/6)) | forced | sixth-power cofactor cores (BFlat §3.4) |
| exact evaluation for a bound | already exact | second Poisson exact (Grouping.lean); limit ODE exact (PR 141 Curve.lean, PR 143 root) |
| Gauss extension cost 4 (k ↦ k(∏p_i)⁴, K′ = … + 4u′) | named outside | the barrier lists "Poisson radius 1 − 6u" as excluded; tested below |
**Test** (code/b1_twist.py, b1_chain.py; logs/b1_twist.log, b1_chain.log; 20–24 MB footprint, 36 s + 94 s). With K′ = 2(N+V) − 1 + κ + (4 − ρ)u′:
Hc′ = 1 − (6 − ρ)u′ − κ + κP, payment K′ − N′ − V + Hc′ = N + V − u′ + κP unchanged, S_r = H·S_R + u/k + u, so
Y = Hz + u(1/k − 2 + ρ/2), u = (1 − H)/(6 − ρ); limit r z′ = z − (1/thr − 2 + ρ/2)/(6 − ρ), start r0 = (28 − 10ρ)/(36 − 15ρ); thr unchanged.
| ρ | r0 | z(0.7523) | θ (ODE) | θ − θ_∞ | finite chain (step 1/40000) |
|---|---|---|---|---|---|
| 0 | 7/9 | 8.872e−4 | 0.8733458561067 | −5e−15 | 0.873345895 |
| 24/289 = ρ_max(9/2) (V14) | 0.781760 | 1.140e−3 | 0.8733349830 | −1.09e−5 | 0.873335027 |
| 0.1304 (V14's formula at p = 2) | 0.784163 | 1.303e−3 | 0.8733279269 | −1.79e−5 | — |
| 0.1452 (formula at p = 3) | 0.784931 | 1.356e−3 | 0.8733256007 | −2.03e−5 | 0.873325649 |
| 1 (radius 1 − 5u) | 0.857143 | 7.737e−3 | 0.8730268822 | −3.19e−4 | — |
| 1.042 | | | 0.873 | | |
Gates (code/b1_gates.py, logs/b1_gates.log; 7 MB): PR 137's 14 stage gates (chain.py = ChainData) along the zero-pad twisted chain at
PR 137's step 1/4000: ρ = 24/289 all hold (min 2.3e−4, gate 2 at the left end 0.7335; gates 12, 14 8.2e−3 at the new top 0.7815);
ρ = 0.1452 only gate 2 fails, by 1.1e−5 at r = 0.7337, far from the count crossing 0.752.
**Correction to B2Below.** Its row "radius 1 − 5u → 0.8733384" kept the coefficient 1/k − 2 (chain_limit.py udiv); the consistent coefficient
is 1/k − 3/2 (the omitted gain is (1 − H)/10; reviewer D). The radius is the most potent stage constant: 1 − 5u gives 8.7 z_∞(c), not 1.2 z_∞(c),
and 0.873 needs extension cost 2.96 instead of 4 (dθ/dρ ≈ −1.3·10⁻⁴ near ρ = 0). Every other entry of B2Below's lever map is unaffected.
**Lead (not a freedom of PR 137's lemmas).** Liu V14 (Hexagon 2610.00655v3) uses the same identity G(n, a⁴h) = χ_n(a)²·G(n, h) (l.580–581, 9096–9097)
to keep a part s of the unit primes as χ_n(a)² instead of extending the frequency by a⁴: K − 4s, Cauchy loss s over free twisted tuples,
moving radical M_* + 5s; net origin saving 3s = ρW, ρ < ρ_max(p) = 6(p−6)(2p−3)/(p² − 193p + 198) ∈ [24/289, 5662/67351] for p = 6κ ∈ [112/25, 9/2],
fixed by its returned-sum bound z < (1/p + α)G (Lemma tw:lem:supremum). B_2's chain does not use V14 (Hexagon2 §1).
A certificate would need: (i) a twisted selection in Select/Scale (part of the slot product kept as χ_n(a)², Cauchy over free twisted tuples);
(ii) anchors uniform over moving cubic radicals (ChildBound is uniform only over ray characters mod 4, IsChildRay); (iii) a returned-sum bound
that fixes PR 137's own ρ_max (V14 Cor tw:cor:return-bound analogue; V14's formula is proved only for p ∈ [112/25, 9/2]); (iv) the chain with
H = 1 − (6 − ρ)u and every gate rechecked (all have ≥ 1e−4 slack in the ρ = 0 limit), and the contact algebra of PR 141/143 with the new linear ODE.
**Verdict.** ROBUST: no freedom of PR 137/141's own lemmas passes θ_∞; the barrier's named exclusion (radius constant) is where the lever is.
**Checks.** Reviewer (logs/reviewer_b1.txt): CORRECT on A–D (bookkeeping; θ floats not re-run by it). Novelty: tracker PRs ≤ 143 (PR 143,
θ_∞ itself, merged 15 Nov 22:38 UTC; TREE N7 still says "open"); arXiv math.NT 2026-10/11 and Zenodo "quasi-Riemann" (lst.py): no twisted
selection or bound below θ_∞; our files: no prior note of the radius correction.

## B2 — bflat: quasi-RH half-plane below B♭ = 0.874886647 (NOGO l.253) — ROBUST in its class; passed by one listed input
**Step fixed.** (A3) at the critical class endpoint (δ_c, x) = (0.4693, 1/2) (double root in δ): count R = 0.590332 set jointly by M1
(Möbius second moment, prime cost 2), M2 (plain fourth moment, prime cost 6κ) and sixth-power amplification (5/6); numerator by Hölder
with M2 ((3R+1)/4 − θ/2, θ = g/(6κ)); stage map Φ(B) = F♭(2B − 1) with κ licensed by β_* ≤ B.
**Input class.** NF3 (2610.26835) Lemma 17.1 (r + 2z ≤ 1), Thm 12.16 (n1 + n2 + 6κz ≤ M; arbitrary lengths at z = 0), Lemma 17.6
(amplification), contour lemma (θ ≤ 1/10), detectors of product, two-piece or mollified form (QRHImprove REPORT, REPORT3).
**Freedoms the lemmas admit.** Thm 12.16's arbitrary lengths at z = 0 and unequal n1, n2 (the face depends on n1 + n2 only), several
cutoffs, Hölder mixtures, a trivial second factor, Heath-Brown splits, the 2/27 capacity (κ ≥ 3/4; equals 1/(18κ) at κ = 3/4), κ
licensing (F♭ flat: 0.874312 even at κ = 1/2), M1's side condition 2r + 8z ≤ 3, twisted M2 (prime product Q inside the fourth moment;
already used as the shift θ): all in class and none better (BFlat NOTE §3; QRHImprove REPORT3). No unused own-lemma freedom found.
**Test** (code/b2_bflat.py, logs/b2_bflat.log; 8 s, 21 MB): BFlat's own pipeline (bflat_sens.py definitions, executed read-only), only
M2's prime cost 6κ replaced by a constant c:
| c | 6κ (NF3) | 9/2 | 4 | 3 | 5/2 | **2** | 7/4 | 3/2 |
|---|---|---|---|---|---|---|---|---|
| fixed point | 0.8748866470 | 0.8748869894 | 0.8747451 | 0.8743125 | 0.8739542 | **0.8733819666132** | 0.8729322 | (0.87224: another class binds) |
At c = 2 the fixed point (B, δ_c, b) equals (B_2, δ_2, b_2) = ((777 − √1081)/852, (61 − √1081)/66, −11/71 + 595√1081/76751) to 10⁻⁴²
(code/b2_exact.py, logs/b2_exact.log, 40-digit findroot; 1278B² − 2331B + 1061 = 0 there):
Liu's B_2 (Hexagon 2611.00082, the Lean-checked PR 103/129 contact) is BFlat's stage map with M2 cost 2 (2611.00082's second moment
induction n1 + n2 + 2z ≤ M, Hexagon2 §1), nothing else changed. At c = 2: M1 cost 19/10 → 0.8732001, amplification 4/5 → 0.8732298;
θ_∞ corresponds to M2 cost ≈ 1.98. Every record between (V14's B_*, B_♮, B′, B_ex, B_tw; Chang–Wang 0.874462; B_2) came from new
M2-type norms, an enlarged M2 count (V14 Prop ec:count, extra width W_cap for m > 3/8) or a numerator contour, i.e. from BFlat's own
"what a record needs" list.
**Verdict.** ROBUST in its class. The NOGO line's parenthetical "(amplification exponent < 5/6)" named the wrong input: B♭ fell to the
M2 prime cost, whose per-unit sensitivity (+0.00019) looked smallest but which moved by 2.5 units. Lesson-235 wording below.

## B3 — the 13/15 floor (R109; QRHImprove REPORT3; R291 "cubic-theta product-moment method cannot pass 13/15") — ROBUST
**Step fixed.** Two constraints meet at ℓ = 1/5, b → 0: (G) the count-free generic-row end C0 = −1/4 + b/6 + 5ℓ/4 < 0 (Liu (14.7) at
δ → 0: a = 1/2, R = 1, K = 1, g = 0), which caps ℓ < 1/5 − 2b/15, against the direct line B = 11/12 − ℓ/4; so B > 13/15 + b/30.
For ℓ > 1/5 the direct side's reflected-energy excess (5ℓ − 1)₊/4 (LP-attained over Lemma 14.3, REPORT3) forces 13/15 again.
**Input class.** Rows summed in absolute value; any row-wise input up to GRH/Lindelöf for every row L-function; counts and moments of
any strength; product, two-piece and mollified detectors; NF3/OpenAI Lemmas 14.2–14.3 on the direct side.
**Freedoms.**
| freedom | verdict | reason |
|---|---|---|
| stronger counts, moments, large sieves, PR 137-type mean squares | provably nothing | they act through R, but 1 − δ ≤ R ≤ 1 gives R = 1 at δ → 0; C0 holds even under GRH per row |
| contour/reflection shift θ | nothing | θ ∝ g = 0 at the generic-row end |
| detector restructuring | nothing | REPORT3 (a), (b) |
| direct-pairing phase saving λ (BeyondB2 S7) | would move it one-for-one | floor 13/15 + b/30 − λ, since C0 involves only (b, ℓ); no input gives λ > 0 |
| trilinear step on the direct side (REPORT3 sketch, 28/33) | not alone | (G) still caps ℓ at 1/5 |
| signed full-family sum at the generic-row end | tested | would relax (G); the original excess still forces 13/15, the trilinear one 28/33 |
**Test** (code/b3_family.py, logs/b3_family.log; 8 s, 19 MB): BeyondB2's finite model (all characters mod prime q, exact L by FFT),
full family, T = L(1/2 − θ, χ)/L(0.75, ηχ̄). For q = 12289, 40961, 65537, 163841, θ ∈ {0, 0.047}, η ∈ {χ_7, 1}: ΣT / diagonal
((#rows)/L(1.25 − θ, η)) = 0.987–1.025; r_full = |ΣT|/Σ|T| = 0.11–0.18 for η = 1, flat in q; after removing the diagonal r_off =
1.8e−4 … 6.1e−3 with r_off·√q = 0.07–0.68. The family first moment is its diagonal: signs buy no exponent at R = 1 unless the actual
(14.7) summand (theta and Gauss coefficients, after the detector's probe subtraction) has a vanishing diagonal, which our files do not
decide. Proxies: characters mod q, not sextic rows.
**Verdict.** ROBUST in its class, doubly forced; below 13/15 needs λ > 0, or a vanishing generic-row diagonal plus REPORT3's trilinear step.

## B4 — BeyondB2: "the one provable signed tool, a large sieve in the Möbius variable (LS1), loses Z^0.0333 at the contact" — ROBUST
**Step fixed.** The dual class sum Σ_{u∈𝓑} b_u F_u, F_u = Σ_{m~Z} μ(m)ηχ̄_u(m)w(m/Z), at the contact class of B_2 (δ_2, 1/2);
LS1 = duality plus the large sieve in m: |Σ b_u F_u| ≤ Z‖b‖₂, against the absolute route Z^a‖b‖₁ (a = (1+δ)/2).
**Freedom tested** (a different inequality in the same step): Vaughan/Heath-Brown split of μ into Type II pieces Σ α_a β_b χ_u(ab),
Cauchy–Schwarz and a family mean value with ANY constant Δ(N) ≥ max(U, N) (even the conjectural diagonal sextic large sieve), or a
subset large sieve Δ = N + #𝓑·U^{1/2}.
**Test** (code/b4_typeII.py, logs/b4_typeII.log; 19 MB). Z-exponents at the contact (a = 0.713041, h = 0.809686, R = 0.626569):
absolute a + hR = 1.220366; LS1 1 + hR/2 = 1.253662 (loss 0.0332965 = BeyondB2 C1); Type II ≥ 1/2 + h = 1.309686 (loss 0.0893);
subset large sieve 1.412168 (loss 0.1918). On the whole class domain [1/8, 5/6] × [0, 1/2] (native count) Type II never beats the
absolute route (min loss 0.0085 at δ = 1/8, x = 0); it would need h(1 − R) < δ/2 (contact: 0.302 > 0.213).
**Verdict.** ROBUST: the family size U = Z^h, h > 1/2, caps every mean-value signed tool; a dual saving needs cancellation that is
not a mean square (BeyondB2 §3: phase-aligned subclasses carry the same magnitude data).

## Suggested NOGO lines
- b2below (N7, restated with its class): θ_∞ = 0.8733458561067 is the ceiling of PR 137/141's lemma class (stage lemma with Gauss extension
  4u′ and radius 1 − 6u; HIGH gate G(k) = H·F0, independent of the radius; count envelope 𝒞); HIGH line, thresholds and cutoffs add nothing;
  the radius is the potent constant: extension cost 4 − ρ gives Y = Hz + u(1/k − 2 + ρ/2), θ = 0.8733350 at ρ = 24/289 and 0.873 at
  ρ = 1.04 (B2Below's 1 − 5u row understated it 42×); lead: a V14-type twist χ_n(a)² in the selection [CeilAudit2](Riemann/Deep/CeilAudit2/REPORT.md)
- bflat (N7, restated): B♭ = 0.874886647 is pinned in NF3's class (M1 cost 2, M2 Thm 12.16 cost 6κ, sixth-power amplification 5/6,
  Hölder–M2 numerator, contour θ ≤ 1/10); passed only by new M2-type inputs: BFlat's own map at M2 cost 2 is B_2 = (777 − √1081)/852
  exactly [CeilAudit2](Riemann/Deep/CeilAudit2/REPORT.md)
- 13/15 floor (N7, restated): the direct line 11/12 − ℓ/4 at the count-free generic-row edge C0 = −1/4 + b/6 + 5ℓ/4 = 0 (ℓ = 1/5), and
  independently the reflected excess (5ℓ − 1)₊/4; robust to every row-wise input up to GRH, every count, moment or detector, and to signed
  full-family sums whose diagonal survives (finite model: sum = diagonal within 2.5%); moves one-for-one with a direct-pairing saving λ, or
  with a vanishing generic-row diagonal plus REPORT3's trilinear step (then ≥ 28/33) [CeilAudit2](Riemann/Deep/CeilAudit2/REPORT.md)
- beyondb2 signed tools (N7): at the B_2 contact every mean-value signed tool loses to the absolute route: LS1 Z^0.0333, Vaughan Type II with
  any large-sieve constant ≥ U Z^0.0893, subset large sieve Z^0.1918 (family size U = Z^0.81 > Z^{1/2}) [CeilAudit2](Riemann/Deep/CeilAudit2/REPORT.md)

## For the orchestrator
- Tracker: PR 143 (exact continuum contact, θ_∞ itself) was merged 15 Nov 22:38 UTC; TREE N7 still calls it open. Newest PR is 143.
- Not a discovery (corrections and identifications inside our audits); no Lean; nothing edited outside this folder and .scratch/ceilaudit2/.

RH relevance: N7

## Files
| file | role |
|---|---|
| code/b1_twist.py | B1: continuum ODE with Gauss extension cost 4 − ρ; count window + return via B2Below (system.py, infimum.py, read-only) |
| code/b1_chain.py | B1: finite zero-pad chain vs ODE; ρ for θ = 0.873 |
| code/b1_gates.py | B1: PR 137's 14 stage gates along the twisted chain |
| code/b2_bflat.py | B2: BFlat's stage map (bflat_sens.py definitions, read-only) with constant M2 costs; lever map at cost 2 |
| code/b3_family.py | B3: full-family signed sum vs diagonal in BeyondB2's finite model (finite_model.py rows, read-only) |
| code/b4_typeII.py | B4: LS1, Type II and subset-LS exponents at the contact and on the class domain |
| logs/*.log, *.err | outputs; /usr/bin/time -l footprints 7–24 MB (all < 300 MB); reviewer_b1.txt |
Run: uv run --no-project --with mpmath --with numpy [--with sympy] python3 code/<file> (each < 2 min, single thread).
