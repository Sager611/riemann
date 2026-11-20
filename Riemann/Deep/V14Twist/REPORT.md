COUNTEREXAMPLE: Liu's V14 twist loses at every stage of PR 137/141's marked-anchor chain. Its own price, the Cauchy count U^τ over free twisted tuples (V l.9103), comes out of the child's excess. At equal radius a twisted stage therefore loses (τ/3)(1/k − 2) > 0 whenever k = thr < 1/2, which is every stage below R = 7/9. Honest continuum θ = 0.8733462063 > θ_∞; CeilAudit2's 0.8733350 left out the count. Passing θ_∞ needs a new twisted-family bound with any saving Z^{γτ}, γ > 0.

# V14Twist: import Liu's V14 unit twist into PR 137/141's marked-anchor chain (node N7)

(v14twist, Orca, 20 Nov 2026, 04:13–04:50 UTC.) Sources: Hexagon 2610.00655v3 TeX (curl https://hexagonmath.org/src/2610.00655v3, zip;
number_field_new_v1-4.tex, 11726 lines, "V l."). Tracker GTimG/QuasiRiemannTracker main, PR 137 = submissions/2026-11-12-marked-anchor-induction
(145 files of MarkedAnchor/{Inverse,Interface,Anchor,Driver}, High/Algebra.lean, *Twist.lean, read via raw.githubusercontent.com).

## Statement
1. **What V14 twists, in PR 137's terms.** PR 137's "fourth-power selection" is V14's frequency extension. The selected tuple d (width u′) is absorbed as
   k ↦ k·P_d⁴ through the injection (d,k) ↦ kP_d⁴ (Select/Separate.lean selected_energy_bound; Scale.lean hmap, hmult). V14 keeps part τ of it as χ_n(a)²,
   using G(n, a⁴h) = χ_n(a)²G(n,h) (V l.9096). The twisted tuple is still removed from the column. It is not extended, and its Cauchy count is Z^τ (V l.9103).
2. **Hypotheses a twisted stage needs** (exact list):
   (H1) unit locus: automatic in PR 137 (a column meeting img d gives χ_{U′}(kP_d⁴) = 0, Separate.lean row_full_pool);
   (H2) child bound with a moving radical: ChildBoundInj for column characters ψ·χ̃_P, χ̃_P(A) = χ_A(P)⁴, N P ≤ Z^τ, at radius Hc′ = 1 − 6u′ + 4τ (+κ terms),
        bound Z^{Hc′(1+x)} with x·Hc′ ≤ u′ − τ. This follows from the payment K′ − N′ − V + Hc′ = N′ + V + κP (K′ cancels) plus the count τ.
        NOT available: RestoreChildAt (Assembly.lean l.61–82) allows only base characters that factor modulo a fixed q chosen before Z, and IsChildRay
        allows only ray characters mod 4;
   (H3) stage algebra with excess e = u′ − τ: highPriceF(a0, a1, k, r, H·S_R + e/k, H) ≤ e, with the selected width u′ in S_r = s + u′. Gate 1 is unchanged
        (G k = H·F0), so thr is unchanged;
   (H4) diagonal: N + V + 4u′ − 3τ + κ ≤ 1. This is V14's relaxed origin A_max + (4−ρ)W < 1 with 3τ = ρu′. It is the ONLY place where V14's saving
        appears, and it does not bind in the continuum limit;
   (H5) V14's ρ < ρ_max(p), p = 6κ ∈ [112/25, 9/2], is the closure of V14's returned sums (E_α ≤ 0, V l.9276; E_0 = 7ρ/18 > 0). It has no counterpart here,
        because PR 137's children are anchors, not M2 sums, and its κ is a slack, not a prime-length constant. PR 137's own closure condition is item 3 below
        (k ≥ 1/2, false on the chain), so its ρ_max is 0.
3. **The decisive step (falsification).** Under (H2)–(H3) the stage is z_r = Hz + e(1/k − 2) with 1 − H = 6e + 2τ. A twisted width is a pure radius loss:
   the stage has the same additive gain as an untwisted stage with u = e, but a radius 2τ smaller. Per unit of width: an extended unit costs radius 6 and
   returns capacity 1 + excess 1 (worth 1/k), i.e. T_0 = (1/k − 2)/6 per unit of radius; a twisted unit costs radius 2 and returns capacity 1 + excess 0, i.e. 0. At the same H,
   z_r(twist) − z_r(no twist) = −(τ/3)(1/k − 2). That is negative exactly when k < 1/2. Along the chain thr ∈ [0.3626, 1/2), and exactly:
   z_∞(r) < 7/6 − 3r/2 ⇔ 3(3r−2) < 9r/7 ⇔ r < 7/9. So V14's hypothesis "the twist is paid" fails at the first stage below R = 7/9 and at every later
   stage, for every τ > 0, i.e. every selection of twisted positions.
4. **Barrier extends (exact).** On the comparison domain T(R, z_R) ≥ 0, so for k ≥ thr every twisted stage satisfies
   Y_tw ≤ Hz_R + T(R, z_R)(1 − H). This is B2Below §4 step 1, and its steps 2–4 then give: every finite chain built from PR 137's stage lemma with V14 twists
   has gain ≤ z_∞, and the half-plane Re s > θ it gives has θ ≥ θ_∞ = 0.8733458561067185768.
5. **What would pass θ_∞ (lead, unproved).** A new twisted-family estimate: Σ over the twist family ≤ Z^{(1−γ)τ}·(child bound). With such a saving the optimum
   is ρ′ = τ/u′ = 1 (no Gauss extension at all), where T = γ/(2k) > 0 even at k = 1/2. ANY γ > 0 then passes θ_∞:
   γ = 0.03 → 0.8733327 (float). V14's Thm tw:thm:moving (V l.10930) is the nearest proved analogue, but it covers V14's returned sums, not marked anchors.
   The natural tools give γ ≤ 0 in the anchor regime (rows ≫ columns).

## Chain (continuum limit of PR 137's stage, PR 141's flow; count window + return as B2Below/ode.py; code/v14_import.py)
| mode (bookkeeping) | T in r z′ = z − T | start r0 | z(0.7523) | θ (float) | θ − θ_∞ |
|---|---|---|---|---|---|
| untwisted (PR 137/141) | (1/thr − 2)/6 | 7/9 | 8.8719e−4 | 0.8733458561067 | −5e−15 |
| **V14 twist, honest** (count τ in excess), σ = τ/e = 0.028470 (ρ′ = τ/u′ = 8/289 = ρ_max(9/2)/3) | (1/thr − 2)/(6 + 2σ) | 7/9 | 8.7898e−4 | **0.8733462063** | **+3.50e−7** |
| same, σ = 0.085409 / 0.5 | 〃 | 7/9 | 8.6300e−4 / 7.6211e−4 | 0.8733468868 / 0.8733511556 | +1.03e−6 / +5.30e−6 |
| V14's own count R_loc + 2s (e = u′ − 2τ), σ = 0.028470 | (1/thr − 2 − 2σ)/(6 + 8σ) | 0.772744 | 5.9044e−4 | 0.8733584122 | +1.26e−5 |
| CeilAudit2 model (cost 4 − ρ, no count), ρ = 24/289 | (1/thr − 2 + ρ/2)/(6 − ρ) | 0.781760 | 1.1402e−3 | 0.8733349830 | −1.09e−5 |
| hypothetical count-free twist (γ = 1), ρ′ = 8/289 | (1/thr − 2 + 2ρ′)/(6 − 4ρ′) | 0.783152 | 1.2336e−3 | 0.8733309399 | −1.49e−5 |
| hypothetical saving γ = 0.03 / 0.05 / 0.1 / 0.2, ρ′ = 1 | max(T_0, γ/(2thr)) | 7/9 (pinned) | 1.196e−3 / 1.767e−3 / 3.522e−3 / 6.995e−3 | 0.8733327 / 0.8733083 / 0.8732278 / 0.8730627 | −1.3e−5 / −3.8e−5 / −1.2e−4 / −2.8e−4 |
Finite zero-pad chain (step 1/4000, k = thr, wb = R; code/v14_checks.py): z(0.7523) = 8.7971e−4 (untwisted), 8.7156e−4 (σ = 0.028470), 8.5571e−4 (σ = 0.085409);
z(0.7335) = 3.1302e−3, 3.1016e−3, 3.0461e−3. Break-even family saving γ* = (1 − 2k)/3 ∈ [0, 0.0916]; 0.0494 at the contact crossing (k = 0.425831).
No certified chain: the twist does not hold at any stage (Statement 3), so item (3) of the brief does not apply. The tracker-submission draft is moot.

## Evidence (logs; all single-threaded under bash .scratch/memcap.sh 290, footprints from /usr/bin/time -l)
- logs/v14_import.log: the table above, 8 runs (memcap peak 36 MB, max RSS 50 MB, 23 s). logs/v14_gamma.log: γ rows (20 MB, 47 MB, 15 s).
- logs/v14_checks.log: thr range [0.362614, 0.500000], max_r z_∞ − (7/6 − 3r/2) = −1.5e−9 (at r → 7/9), finite chains, γ* (44 MB, 15 s).
- logs/v14_exact.log (sympy): the closed form z_∞ = (33r−14)/48 − (35/144)(9r/7)^{4/5}(9r−6)^{1/5} solves r z′ = z − T exactly (residual 0), z(7/9) = 0, and the
  sign of the margin matches [3(3r−2) < 9r/7] on (0.6671, 0.7780) (60 MB, 8 s).

## Checks
- The pipeline reproduces θ_∞ (untwisted −5e−15; PR 143's root) and CeilAudit2's 0.8733349830 to every printed digit (ca2 mode, same code).
- Finite chains and the ODE agree in sign at every σ: the twist lowers the gain in both.
- The stage identity z_r = Hz + e(1/k − 2) was re-derived by hand from High/Algebra.lean l.43–59, and again by the reviewer (B). The margin is exact.
- Lean sources read: Select/Scale.lean (selection_exponents l.24–51, numerics l.318–445), Select/Separate.lean l.291–420, Select/AtScale.lean (plan, mass_le),
  Select/Assembly.lean RestoreChildAt l.61–82, Poisson/Interface.lean (IsChildRay, ChildBound, prefactor_identities), High/Algebra.lean l.31–61,
  High/Marked/Twist.lean and Restore/Twist.lean (fixed-modulus ray-class twists only).

## Review
Reviewer subagent (logs/reviewer.txt; brief with every formula inline: logs/reviewer_brief.txt). A CORRECT for PR 137's existing Cauchy/fibre/per-twist-max
argument: "the residual Z^τ gives xHc′ + τ ≤ u′; it also multiplies the diagonal, whose exponent becomes K′ + τ, improving by 3τ". B CORRECT: "strict negativity
requires τ > 0 and 0 < k < 1/2". C CORRECT "within the stated comparison domain (R > 2/3, T ≥ 0, e ≥ 0, 0 < H ≤ 1)". Overall: CORRECT as a conditional,
method-specific reduction. "Unavoidable" concerns the existing argument, not a hypothetical joint-twist estimate (our Statement 5). The reviewer did not certify
the numerical margin; we then proved it exactly (v14_exact).

## Novelty
- Tracker (GitHub API, 20 Nov ~04:45 UTC): the newest PR is #143 (merged 15 Nov 22:38). None came later, and no submission uses a unit twist.
  PR 137's Twist.lean files are fixed-modulus ray-class twists of slot coefficients (SlotTwist with N(rayM′) < D; childTwist = conj θ·δ·η′).
- arXiv math.NT 2026-10/11 (lst.py 'zero-free|quasi-Riemann|half-plane'): 6 hits, all already in our files (2610.12234, 2610.15291, 2610.18491, 2610.26835,
  2611.00924 Chang–Wang, 2611.03356). Zenodo '"zero-free" "Dirichlet"': nothing on N7.
- Our files (rg 'V14|tw:thm:norm|0.87333'): the only prior is CeilAudit2 B1's lead (0.8733350, now corrected) plus Hexagon/Hexagon2. No file has the
  count-versus-radius bookkeeping.

## Caveat
1. Scope: the falsification covers PR 137's selection as written, i.e. the central Cauchy over the whole selected tuple followed by a per-twist maximum,
   together with V14's stated price. It does not exclude a new joint estimate over the twist family with γ > 0 (Statement 5). Even that estimate would
   still need (H2), uniformity over moving radicals, which RestoreChildAt does not give (q is fixed before Z).
2. (H3) is not Lean-checked: PR 137 proves highPriceF ≤ u only for u = (1−H)/6. We use only that thr is unchanged (G k = H·F0, exact) and the lemma's
   capacity form s = H·S_R + e/k, the same convention CeilAudit2 and B2Below used for their levers.
3. The θ columns are floats (RK4 plus golden search, 1e−16 against PR 143 in the untwisted case). The barrier (Statement 4) rests only on the exact stage
   inequality, the exact z_∞ and B2Below §4's comparison lemma (two reviewer rounds there), with stage radii H ∈ [47/50, 1] as in PR 137.
4. The modes 'full', 'fixed' and γ are bookkeeping variants, not proved bounds.
5. Correction to CeilAudit2 B1 "Lead": "Gauss cost 4 − ρ with the payment unchanged" counts V14's net origin saving 3s as a radius saving. In PR 137 the
   saved 4τ goes to the radius (worth 2τ) and the count τ comes out of the excess (worth τ/k > 2τ).

## For the orchestrator (suggested NOGO line; not written by us)
- v14twist (Orca, 20 Nov 2026): V14's unit twist χ_n(a)², put into PR 137/141's selection, loses at every stage. The payment N′+V does not depend on K′,
  so the saved extension 4τ only enlarges the radius (worth 2τ), while the Cauchy count Z^τ over free twisted tuples comes out of the child excess
  (worth τ/k > 2τ, k = thr < 1/2). So Δz = −(τ/3)(1/k − 2) and θ ≥ θ_∞ for every twisted chain. CeilAudit2's 0.8733350 left out the count.
  Exit: a twisted-family bound saving Z^{γτ}, any γ > 0 (float: γ = 0.03 gives 0.8733327).

RH relevance: N7. The ceiling θ_∞ of the zero-free half-plane for ζ and every Dirichlet L-function survives V14's twist. The only exit we found is a power
saving over the twisted-tuple count.

## Files
Riemann/Deep/V14Twist/: REPORT.md; code/v14_import.py (all modes, ODE + window + return), code/v14_checks.py, code/v14_exact.py, copied code/ca2_b1_twist.py
(CeilAudit2), code/system.py and code/infimum.py (B2Below); logs/v14_{import,gamma,checks,exact}.{log,err}, logs/reviewer.txt, logs/reviewer_brief.txt.
Scratch: .scratch/v14twist/STATE.md (downloads deleted).

## Notes (written as I go)
- V14 setup (V l.8907–8980): p = 6κ ∈ [112/25, 9/2], ϱ = 1 − 1/p, α_max = (4p−6)/(p(p+24)), ρ_max = 6(p−6)(2p−3)/(p²−193p+198) ∈ [24/289, 5662/67351].
  Thm tw:thm:norm (l.8953): under Thm en:variable-rectangle-norm's hypotheses, origin A_max + 4W < 1 → A_max + (4−ρ)W < 1, A_max = 2r + (1−2r)/p + W;
  retained: A_max + 3W < 5/6 + 2t/3, 2 + 2t − 4m_− + pW < 1, W ≤ 1/25. Input: returned-sum bound (tw:eq:returnbound) z < (1/p + α)G, 0 < α < α_max(p)
  (Cor tw:cor:return-bound) replacing z < G/p of Lemma en:exhaustion.
- Margin (l.8982–9013): Δ = ϱ(1−2r) − (5−ρ)W > 0; for m ≤ r, w ≤ W, z ≤ (1−2m)/p: d_m = 1−2m, G_0 = d_m − 6w, A_t = 2m + z + w;
  A_t − 1 + 4w ≤ ρw − ν (second origin), G_0 ≥ k_p w + ν/ϱ, k_p = (6 − p − pρ)/(p − 1) (homogeneous return).
- Selection (l.9015–9086): origin loss 𝒮 = 4w_f + 5f_1 + 5v + 4c_cor + 4(T_1+T) + S_1 + S + e_1 + e + μ + B_C; second-origin excess ≤ A_t − 1 + 4w − 𝒮.
  No twist if ρw ≤ ν/2 or w_f + f_1 + v ≥ ρw/4 (then 𝒮 ≥ ρw already). Else select whole UNIT positions (fixed or free, before extension)
  of total width s = ρw/3 + ϑ_sel; needs unit width w − w_f − f_1 − v > (1 − ρ/4)w ≥ s; s/w < 3/100.
- Operation (V l.9088–9157): unit locus G(n, a⁴h) = χ̄_n(a)⁴G(n,h) = χ_n(a)²G(n,h) (cubic); χ_n(a)² joins the common character of the remaining factors.
  Free twisted tuple: Cauchy loses U^{s+ε}. Shifts K ↦ K − 4s, m_* ↦ m_* + 4s, q_* ↦ q_* + s, M_* ↦ M_* + 5s; R_loc ↦ R_loc + 2s, O_total ↦ O_total − 3s
  (so the origin saving is 3s = ρw); exceptional exponent −7s/3; terminal −3s, −13s. The twisted positions are NOT extended.
- PR 137 is the same Eisenstein-cubic setting (Select/Scale.lean: O = ActualEisensteinCubic.O, primes p_i ∉ (λ), λ² | p_i − 1): its "fourth-power
  selection" is V14's extension k ↦ k·T⁴, T = ∏ selected p_i, cost 4u′ in K′ = 2(N+V) − 1 + κ + 4u′ (selection_exponents l.24–51).
  NB PR 137's κ is a small slack (5·mesh ≤ κ, Pmax + Smax + 4u + 10κ < 1), NOT Liu's prime-length constant p = 6κ.
- PR 137 child (Inverse/Poisson/Interface.lean): P_y = Σ_{A ⊆ pool∖t} μ(A) ψ(A) mask_A(m) χ_A(y)^{±1} mark(A) Φ(NA/Xc)/√NA over Ny ≤ Z^{H+κ}, H = 2(N′+V) − K′;
  ChildBound quantifies over ψ with IsChildRay Ψ ψ :⇔ ψ = (ray character mod 4)·Ψ only. prefactor_identities: K′ − N′ − V = K − F + 5u, H = 1 − 6u,
  payment (K − F + 5u) + (H + u) = F.
- V14 return (l.9158–9287): returned sums pay the twist: B_ret = C_hom + 2s + L_ret − Δ′, gap G′ ≥ G_0 + … + 5s; active case coefficient of w:
  E_α(p,ρ) = (7 − 5pα)ρ/18 − pα(6 − p − pρ)/(6(p−1)) ≤ 0 needed. At α = 0 (plain M2, z < G/p) E_0 = 7ρ/18 > 0: V14's twist closes ONLY through
  the new returned-sum bound z < (1/p + α)G (α > 0, Cor tw:cor:return-bound); ρ_max is the root of E_{α_max} = 0.
- PR 137 child bound (Select/Assembly.lean RestoreChildAt l.61–82): ChildBoundInj … Hc′ (C·Z^{Hc′(1+x)+ε}), x = u/(1−6u−κ): the child carries EXCESS
  Hc′·x = u (= the u in "(H + u)" of prefactor_identities). Base characters: FactorsModulo (fixedBaseConductor q) with q FIXED before Z (C, Z₀ depend on q).
- PR 137 stage algebra (High/Algebra.lean l.31–61): u = (1−H)/6; s = H S_R + u/k; g ν = (r + s − H/2 − H a0 − u/ν)/D; G ν = H + u − H/6 − (5/6−ν)(3gν/2) − ν(r+s);
  F0 = 5/6 − (5/4 − 3k/2)wb − k(1+R)/2 − kz; Stg.highPriceF_le: highPriceF a0 a1 k r s H ≤ u. All three u's (in s, g, G) are the EXCESS (price budget):
  at ν = k they cancel (G k = H·F0), so thr does not depend on the excess. The u in S_r = s + u is the SELECTED width. In PR 137 both are (1−H)/6.
- Select/Scale.lean l.365–445 (selection_at_scale numerics): central Cauchy mass Σ_d|∏a(d_i)|² ≤ 128^K·T, T = ∏_R H = Z^{u′} (AtScale mass_le, |a| ≤ 1, so mass =
  the TUPLE COUNT); the inner energy at column scale X/T carries Z^{−u′} (hTinv: T·Z^{−u′} = 1); exponents hB1: (K′−N′−V+εP) + (Hc′(1+x)+εR) ≤ N+V+e with
  hxHc: x·Hc′ ≤ u + x(30mesh+κP). The d-sum is absorbed into the enlarged rows only through the injection (d,k) ↦ k·P_d⁴ (hmap, hmult: fibres ≤ K^K·Cd·Z^η).
- DECISIVE STEP (twisted selection). Twist a part τ of the selected width u′ (V14: s = ρw/3): not extended, removed from the column (V14: "enter only through
  χ_n(a)²"). Then (i) K′ = 2(N+V) − 1 + κ + 4(u′−τ), N′ = N − u′, Hc′ = 2(N′+V) − K′ + κP = 1 − 6u′ + 4τ − κ + κP; (ii) payment K′ − N′ − V + Hc′ = N′ + V + κP
  (K′ cancels: the saved extension goes ENTIRELY into the radius); (iii) the twisted rows Y_{d_tw} are distinct functions, so Σ_{d_tw} is not absorbed by any
  injection: E ≤ mass·mult·#(d_tw)·max A_tw, #(d_tw) ≈ Z^τ for full lists (V14 l.9103: "Cauchy's inequality loses at most U^{s+ε}"). Hence the child's
  excess budget is Hc′x ≤ u′ − τ (V14's own count R_loc + 2s would give u′ − 2τ).
- Twisted stage (High/Algebra with excess u → e = u′ − τ, selected width u′ in S_r = s + u′; G k = H·F0 so thr unchanged):
  z_r = Hz + e/k + u′ − (1−H)/2 = Hz + e(1/k − 2), 1 − H = 6e + 2τ. A twisted width is a pure radius loss: same additive gain as the untwisted stage with u = e,
  radius 2τ smaller. At the same H: z_r(twist) − z_r(untwist) = −(τ/3)(1/k − 2) < 0 iff k < 1/2. Exchange rates: 4τ radius saved is worth 2τ of capacity;
  τ of excess is worth τ/k > 2τ. Continuum: r z′ = z − T_tw, T_tw = (1/thr − 2)/(6 + 2σ), σ = τ/e (V14 fraction: σ = (ρ/3)/(1 − ρ/3) = 0.028470).
- On the chain thr ∈ [0.3626, 1/2] (k = 1/2 only at the start 7/9): checks (1) max_r z_∞(r) − (7/6 − 3r/2) = −1.5e−9 < 0 ⇔ thr < 1/2. So EVERY stage with
  positive gain loses under the twist; the first one is the first stage below 7/9.
- Barrier extends: on the comparison domain T(R, z_R) ≥ 0 (⇔ z_R ≤ 7/6 − 3R/2), so for k ≥ thr: Y_tw ≤ Hz_R + T(R,z_R)(1 − H) (if 1/k − 2 ≥ 0 use e ≤ (1−H)/6;
  else Y_tw ≤ Hz_R); this is B2Below §4 step 1, so steps 2–4 give every twisted chain gain ≤ z_∞ and θ ≥ θ_∞. Stages from R ≥ 7/9 (z = 0, k > 1/2) give Y_tw ≤ 0.
- Break-even: a NEW twisted-family estimate saving Z^{γτ} on the count (Σ over the twist family better than #family × max) would make the twist profitable iff
  γ > (1 − 2k)/3: 0.0494 at the contact crossing (k = 0.425831), up to 0.0916 at r = 0.7335. γ = 1 (count-free, "fixed" tuples) = mode 'fixed'.
  Deferred extension in the child (y ↦ P_d⁴y, injective) costs 4τ of excess instead of τ: worse.
- Diagonal: the twisted count also multiplies the paid diagonal Z^{K′}: N + V + 4u′ − 3τ + κ ≤ 1 = V14's relaxed origin A_max + (4−ρ)W < 1 with 3τ = ρu′
  (reviewer: "improving by 3τ"). Non-binding in the continuum limit (u′ → 0 per stage; Pmax + Smax ≈ 0.87 < 1).
- Hypothetical family saving Z^{γτ}: T(ρ′) = [(1 − (1−γ)ρ′)/k − 2 + 2ρ′]/(6 − 4ρ′) is linear-fractional in ρ′ = τ/u′, optimum at ρ′ ∈ {0, 1}; ρ′ = 1 (no extension,
  H = 1 − 2u′, e = γu′) gives z_r = Hz + γu′/k, T = γ/(2k) > 0 even at k = 1/2. So ANY γ > 0 passes θ_∞ near R = 7/9 (logs/v14_gamma.log: γ = 0.03 → 0.8733327).
  Natural tools give γ ≤ 0 in the anchor regime (rows Z^{Hc′} ≫ columns Z^{N′}): the count; deferred extension y ↦ P_d⁴y (cost 4τ); Heath-Brown's cubic large sieve
  on the full ball Z^{Hc′+4τ} (rows term Z^{Hc′+4τ} > Z^{Hc′+τ}). Coefficients with χ_A(P_d) = 1 for all d make the count sharp; full prime lists exclude them.
