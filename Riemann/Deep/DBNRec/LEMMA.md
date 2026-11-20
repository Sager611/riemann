# LEMMA (dbnrec): Raymond's row source with the signed head and the unspent main slope

Source: D. J. Raymond, Zenodo 23752103, `Proof_and_Replay/Release_Support/reconstruction/baseline_originals/Unconditional_0107965.tex` (cited as U.tex l.N).
Scope: the retained affine rows ("R" in U.tex Table l.3619ff; rows 0–11 checked there), scalar (no joint jet: rows ≠ 0, 16, 44).

## Notation (all from U.tex)
- X = 6000000185827, N₀ = 690988, c₀ = ½ log N₀, L₀ = 4c₀ (l.1633, l.3301); L = log(x/4π), c = L/4, Δ = c − c₀ (l.1044, l.1633); τ = x/2 + πt/8 (l.1045, l.3427, l.3543).
- Row i: closed time cell [t₋, t₊], θ = (t − t₋)/(t₊ − t₋); nodes (h_j, a_j, w_j) of M8_LATE_PLAN.json; M = Σ_j w_j > 0 (signed node mass, l.1495, l.3303); F(n) = Σ_j w_j n^{−h_j/2} cos(a_j log n/2) (l.1597, l.3304); s_n = sgn F(n).
- d_t(n): l.1272 (triangular recurrence), d₀ = Λ; d_t(n) ≥ 0 and nondecreasing in t (Lemma src:coeff, l.1266–1269; again l.3311–3320); for a prime p: d_t(p) = (log p) b_t(p), b_t(n) = exp(t log²n/4).
- a_n(t,c) = ½ d_t(n) n^{−1/2−tc}|F(n)|, P_{t,c}(χ) = Σ_{2≤n≤4096} s_n a_n(t,c) χ_n, χ the actual completely multiplicative unit phase, Re χ_n = cos(τ log n) (l.3323–3325, l.3427).
- Signed head: Q_i(t,x) := Re P_{t,c}(χ) = Σ_{2≤n≤4096} ½ d_t(n) n^{−1/2−tL/4} F(n) cos(τ log n).
- Operational fields of MAIN_CLOCK_SOURCE.json: M_i = main_lo, N_i = native_debit, A_i = first_hi, B_i = last_hi; source J̲_i(t) = M_i − N_i − (1−θ)A_i − θB_i (l.997–1004, l.3276–3278).
- β_i ∈ (0,1) rational, the row's certified affine-native slope (l.1495–1509, l.3209–3218; ALL_67_ORIGINAL_AFFINE_NATIVE.json field `beta`; `retained_main_mass` = (1 − β)M checked exactly). s_i := (1 − β_i)M/4.

## Lemma
For a retained scalar row i, every t ∈ [t₋, t₊] and every x ≥ X:

  (L1) J_i(t,x) ≥ M_i − N_i + (1 − β_i)M·Δ − Q_i(t,x),  and  Δ ≥ ¼ log(x/X).

Consequently, if rationals A′ ≥ B′ > 0 satisfy

  (L2) Q_i(t,x) − s_i log(x/X) ≤ (1 − θ)A′ + θB′  for all t ∈ [t₋, t₊], x ≥ X,

then Theorem src:source (l.1011–1023) holds for row i with (first_hi, last_hi) replaced by (A′, B′), all other fields and rows unchanged.

## Proof
Step 1 (source identification, Raymond's own theorem; no new mathematics). U.tex l.3329–3332: "The previously established physical source theorem reduces the source to Mc − Re P_{t,c} minus its fully paid analytic and physical debit. The present section proves improvements of this finite polynomial bound; it does not replace that source-identification theorem." The fully paid debit is the native debit (activation labels n > 4096, reflected labels, physical and analytic errors; l.1099–1230, l.1340–1480, l.3036–3290), charged with absolute values (l.953–957, l.3204–3206). For a retained row its value at L is N_actual(L) ≤ N(0) + βM(L − L₀)/4 (eq. src:affine-native l.1505; l.3215–3217), after the refund, which "changes neither β nor its paid slope" (l.1523–1525, l.3241); the final debit N_i = D_i^final is that refunded anchor value (l.3266–3267: "For rows other than 0,16,44 the final native debit is exactly D_i^sc"), and "the later main-clock gains change only the first and last finite-head costs. They do not change this native debit" (l.3272–3274). Hence, with the main term Mc = Mc₀ + MΔ and M_i ≤ Mc₀ (checked in exact arithmetic, check.py),

  J_i(t,x) ≥ Mc − N_actual(L) − Q_i(t,x) ≥ M_i − N_i + MΔ − βMΔ − Q_i(t,x),

which is (L1) since (L − L₀)/4 = Δ. "The remaining main slope is precisely (1 − β)M; only that slope is available to the later arithmetic credit" (l.3218–3219). Raymond's own first_hi/last_hi are bounds for −(1 − β)MΔ + Q_i (raw chord l.3331–3345 minus credits C_comp + d paid by prime deficits on the near band and by AΔ on the far one, l.3494–3584); the lemma discards those credits and uses (L1) directly.

Step 2 (Δ ≥ ¼ log(x/X)). Δ = ¼[log(x/4π) − 2 log N₀] = ¼ log(x/X) + ¼ log(X/(4πN₀²)), and X/(4πN₀²) > 1 (anchor inequality log(X_*/4π) > 2 log 690988, X_* = X − 4, l.1090–1093; rechecked in Arb, check.py). As (1 − β)M > 0, (1 − β)MΔ ≥ s_i log(x/X).

Step 3 (conclusion). By Steps 1–2, J_i(t,x) ≥ M_i − N_i − [Q_i(t,x) − s_i log(x/X)] ≥ M_i − N_i − (1 − θ)A′ − θB′ under (L2). The flow consumer (make_flow, verify_arcs) uses only J̲_i and asserts first ≥ last > 0 (make_flow_b2.py l.73–76). ∎

## Part B: certifying (L2) (dbnrec's construction; new, to be reviewed)
Fix row i, a shell 𝒮 = [X e^{u₀}, X e^{u₁}] (or u₁ = ∞), L_lo = L(Xe^{u₀}), L_hi = L(Xe^{u₁}). Write m_n(t,x) = ½ d_t(n) n^{−1/2−tL/4}|F(n)| ≥ 0, so Q_i = Σ s_n m_n cos φ_n, φ_n = τ log n.
- Upper chord. For t ∈ [t₋,t₊], x ∈ 𝒮: m_n(t,x) ≤ (1 − θ)H⁻_n + θH⁺_n, where for primes H^±_p = ½ log p·p^{−1/2} e^{t_± log p (log p − L_lo)/4}|F(p)| (exact d_t(p); the exponent is linear in t, hence convex, and decreasing in L since log p ≤ log 4096 < L_lo), and for other n H^±_n = ½ d_{t₊}(n) n^{−1/2 − t_± L_lo/4}|F(n)| (d_t ≤ d_{t₊}; L ≥ L_lo; e^{−t(L_lo/4) log n} convex in t). Raymond's chord l.3331–3345 is the same with c₀ in place of L_lo/4 ≥ c₀.
- Lower bound, constant in θ. m_n(t,x) ≥ ℓ_n: primes ℓ_p = ½ log p·p^{−1/2} e^{t₊ log p (log p − L_hi)/4}|F(p)|; others ℓ_n = ½ d_{t₋}(n) n^{−1/2 − t₊ L_hi/4}|F(n)| (as l.3353, with L_hi); ℓ_n = 0 if u₁ = ∞.
- Termwise. With g_I(φ) = max(lo·cos φ, hi·cos φ) for an interval I = [lo, hi] and I^±_n = s_n·[ℓ_n, H^±_n] (ordered): s_n m_n cos φ ≤ (1 − θ) g_{I⁻_n}(φ) + θ g_{I⁺_n}(φ) (case s_n cos φ ≥ 0: upper chord; case < 0: m_n ≥ ℓ_n, the same in both). And g_{I⁻_n} − g_{I⁺_n} = (H⁻_n − H⁺_n)|cos φ|·[s_n cos φ ≥ 0] ≤ δ_n := H⁻_n − H⁺_n.
- Hence, with Q⁺_𝒮(x′) := Σ_n g_{I⁺_n}(φ_n) (a function of x′ = x + πt/4 alone, since φ_n = (x′/2) log n) and D_𝒮 := Σ_n δ_n:
  Q_i(t,x) ≤ (1 − θ)Q⁻ + θQ⁺ ≤ Q⁺_𝒮(x′) + (1 − θ)D_𝒮.
- bb_rec (dbnhead bb.c + pruning margin) bounds sup over x′ ∈ [Xe^{u₀} + πt₋/4, Xe^{u₁} + πt₊/4] of Q⁺_𝒮(x′), cell by cell; with x ≥ x′ − πt₊/4, each cell C gives Q⁺_𝒮 − s_i log(x/X) ≤ UB_C − s_i log((x′_C,min − πt₊/4)/X), pruned cells ≤ thr − s_i u₀. Let S_𝒮 be the max of these. Then on 𝒮: Q_i − s_i log(x/X) ≤ (1 − θ)(S_𝒮 + D_𝒮) + θ S_𝒮.
- Far range 𝒮_∞ = [Xe^{u_F}, ∞): ℓ_n = 0 and Σ_n m_n ≤ (1 − θ)ΣH⁻ + θΣH⁺, so Q_i − s_i log(x/X) ≤ (1 − θ)(tri⁻_F − s_i u_F) + θ(tri⁺_F − s_i u_F). As used (certify.py): tri^±_F are the last shell's sums, i.e. H^± with L_lo = L(Xe^{u₀,last}) ≤ L(Xe^{u_F}); valid (and conservative) because H^± decrease in L_lo.
- Choice. A′ = max over shells (incl. 𝒮_∞) of the θ = 0 values, B′ = max of the θ = 1 values (rounded up to the producer's 2^−96 grid). Then (L2) holds because (1 − θ)a + θb ≤ (1 − θ)max a + θ max b; A′ ≥ B′ since every D_𝒮 ≥ 0. The shells cover x ≥ X.
- Not used: Raymond's credits (composite, clock, prime deficits) on these rows; they are replaced, not added.

## Review
- reviewer (2026-11-20 09:43 UTC, inequality chain pasted inline): **GAP**, one point: the sup defining S must be taken over the phase variable x′ = x + πt/4 on [Xe^{u₀} + πt₋/4, Xe^{u₁} + πt₊/4], not over the x-shell, with amplitudes from the x-shell; "If 'the shell' already means this enlarged phase interval, this particular gap disappears." (L1) judged correct given the refunded native bound N_actual ≤ N_i + βMΔ, i.e. N_i is the entire remaining debit (l.3266–3274, l.3329–3332).
- Repair: already so. prep_end.py writes xs = Xe^{u₀} + πt₋/4 and xe = Xe^{u₁} + πt₊/4 and bb_rec scans K = ⌈(xe − xs)/D⌉ + 1 anchors from xs; amplitudes use L_lo, L_hi of the x-shell; the slope uses log((x′ − πt₊/4)/X) ≤ log(x/X) (certify.py). With this, the reviewer's verdict is CORRECT.
- reviewer, second pass (10:33 UTC, implementation facts pasted inline: x′ scan range [Xe^{u₀} + πt₋/4, Xe^{u₁} + πt₊/4], x-shell amplitudes, location-aware cells, far chord, endpoint maxima, Arb re-evaluation): **CORRECT** — "the earlier phase-domain gap is closed … yielding (L2) globally." Side claim (a box sup over the whole cell does not justify the brief's uniform shift of first_hi and last_hi): RIGHT, "a uniform box supremum alone does not certify the smaller terminal value of a downward-sloping shifted chord"; another certificate could also do it, this construction is sufficient, not necessary.
