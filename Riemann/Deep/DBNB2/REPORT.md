RECORD: Lambda <= 0.107559925318 (Raymond's chain with Re s > 0.87338196661320752)

# Raymond's de Bruijn-Newman chain restarted from the Lean-checked half-plane Re s > 0.87338196661320752
Orca agent dbnb2, 2026-11-05, 13:50-16:55 UTC (the Mac slept about 1 h 50 min in the middle). Evidence: evidence/ (small files only; the 25 MB flows are regenerated in 100 s each by evidence/make_flow_b2.py).

## Short answer
- **New bound: Lambda <= 2130443810762016970728680497723/19807040628566084398385987584000 = 0.10755992531713451... < 0.107559925318.** It improves Raymond's 0.107964910079 (Zenodo 23752103) by **4.0498e-4**.
- **Only one input changed: the starting height of the zero strip.** Raymond starts the comparison flow at the squared height q = 9/16 + 10^-12, which is what the 7/8 theorem gives. The Lean-checked half-plane gives q = (2 b0 - 1)^2 + 10^-12 = 0.5576563719687857 (b0 = 0.87338196661320752), lower by 4.8436e-3.
- **Nothing else needed to change.** The kernel and own-pair certificates of all 67 probe rows hold for every q in [0, cap], the source certificates do not involve q at all, and the new path stays below every cap (smallest clearance 1.05e-3, against 2.1e-5 before). The new path is 220,181 exact arcs instead of 221,010; it reaches q <= 10^-4 at T = 220181/2048000 = 0.10751025390625.
- **Raymond's own arc verifier passes on the new path**, with only the hard-coded numbers of the old run replaced (diff in evidence/verifier_A.diff): 220,181 arcs, 5,777,345 exact rational checks; his compact-witness stage passes too (3,963,258 checks). **The author's full 15-stage driver, with this arc stage swapped in, passes all 15 stages on this Mac** (2,132 s, python 3.12.12, python-flint 0.8.0); its access logs show that only the arc stage reads the flow.
- **One retained check had to go, and it is vestigial.** On row 55 the author still checks that the path stays in the tube 0.0289 <= q <= 0.0361. That tube belonged to an old local credit that the final source ledger no longer uses (author: "conservative ... not credited a second time"; his source checkers mark row55_local_native_credit_used = False). The new path leaves the tube on the last 30 arcs of row 55 (minimum 0.028862). **Keeping even this check** (by raising the path to the tube floor on 31 arcs, which the comparison allows) gives **Lambda <= 0.10757446927**, still 3.9044e-4 below the old record, and the author's verifier then passes with every one of its checks (variant B).
- **Cross-checks agree.** (1) My producer run with b0 = 7/8 reproduces the author's published path exactly: all 221,010 heights, T, q_T and the endpoint. (2) An independent float model written from the TeX formulas gives both endpoints to 1e-15. (3) The same model at Re rho <= 0.8748867 gives 0.107936602216, matching the 0.107936602217 that Raymond reports for an earlier research run with that strip (TeX line 210): he made the same substitution with the same plan, for a half-plane that was not proved at the time.
- **Sensitivity:** dLambda/db0 = 0.2507 (dLambda/dq0 = 0.0839). The R306 estimate (0.41, from Hakimoglu-Brown's chain) overstated the gain: Raymond's higher rates (12.4-22.7) damp the initial difference by e^-1.8 before the end. Bebel's sharper B_2 + 2.2e-254 lies within 10^-17 of b0 (b0 is B_2 to 17 decimals), so it would move the bound by under 3e-18.

## 1. Where the half-plane enters
Every occurrence of 7/8, 1/8, 3/8 or 3/4 in the TeX and in the 15 stages:

| Place | Role of 7/8 | Changed? |
|---|---|---|
| TeX l.136-146 (intro), l.331-350 (fnd:seveneighths, fnd:initialstrip) | zeros have Re rho in [1/8, 7/8] (1/8 by the functional equation), so H_0 zeros have abs(Im z) <= 2(7/8)-1 = 3/4, y^2 <= 9/16 | yes: [1-b0, b0], abs(Im z) <= 2b0-1 = 0.74676393322641504, y^2 <= q0 |
| TeX l.365-367 (fnd:initialcontraction), l.2301-2330 (startup) | de Bruijn contraction y^2 <= 9/16 - 2t; startup gap 10^-12 - r0 q_init tau, tau = 10^-14 | yes: q0 - 2t; gap 9.3094e-13 (old 9.3036e-13), still > 1/(2*10^12) |
| TeX l.2288-2291 (flow:initial); verify_arcs_fresh.py l.166, 172-173; FLOW.json.gz initial_q | q_init = 9/16 + 10^-12 | yes: q0 + 10^-12 |
| verify_bundle.py l.61-76 (--endpoint) | byte and arithmetic re-check of the published old numbers | not applicable (checks the old run only); the arc verifier checks the new endpoint |
| M8_LATE_PLAN.json "initial_sigma": "7/8" | metadata of the 67-row plan | unchanged; the rows' certificates cover q in [0, cap] |
| 67 row caps (kernel stage: Bernstein cover of q = cap*v, v in [0,1]; own-pair monotonicity on [0, cap]) | caps were set just above the old path (clearance 2.1e-5 at row 0) | unchanged; new path lower, clearance >= 1.05e-3 |
| Row-55 tube check in verify_arcs_fresh.py l.210-212 | relic of a removed local credit (section 3) | variant A: off; variant B: kept |
| Finite head (TeX l.585), compactness (l.839) | use the classical strip abs(Im z) <= 1, not 7/8 | no |
| Sources: native, activation, spectral, heads28/39, wall, overlays, affine, direct-source | no strip input: probes sit at sigma_eff = (1+h)/2 + tL/8 > 5/4, the wall uses Platt-Trudgian; the "/8" in tL/8 and in tau = pi(2e^L + t/8) is the heat-flow normalisation, not the zero-free region | no |

So the half-plane enters only through the start of the flow. The source S_i(t) of each row does not depend on q; the own pair O_i(q) is nondecreasing in q (certified), so a lower path only raises every rate.

## 2. Changed constants (exact)
- b0 = 87338196661320752/10^17 = 5458637291332547/6250000000000000; 2b0 - 1 = 74676393322641504/10^17 (exact; nothing rounded).
- q0 = (2b0-1)^2 = 5445863007497906841530941507209/9765625000000000000000000000000 = 0.5576563719677856...
- q_init = q0 + 10^-12 = 5445863007507672466530941507209/9765625000000000000000000000000 (old 562500000001/10^12).
- First rate r0 = 12.383392719721941 (old 12.380985261471059); startup gap = 36014121503520992269607139677392802690878467780645847777029127/(38685626227668133590597632 x 10^48) = 9.309432214325519e-13.
- Terminal (A): 220,181 arcs, T = 220181/2048000, q_T = 7870749227741121606689559/2^96 = 9.934282176901506e-5, T + q_T/2 = 2130443810762016970728680497723/19807040628566084398385987584000 = 0.10755992531713451; headline 53779962659/500000000000 = 0.107559925318.
- Terminal (B): 220,211 arcs, T = 220211/2048000, q_T = 7854192859820348154813855/2^96, T + q_T/2 = 2130731883412734385049125516963/19807040628566084398385987584000 = 0.10757446926926342; headline 0.10757446927.
- Unchanged: step 1/2048000, cushion 10^-12, eta = 10^-12, tau = 10^-14, 2^-96 outward dyadic rounding, stop at the first arc end with q <= 10^-4, all 67 rows (plan sha256 ec346b7c..., source b59905d7...).

Hand check of the startup stage: r0 q_init tau = 12.38339272 x 0.55765637 x 10^-14 = 6.9057e-14, so the gap is 10^-12 - 6.9057e-14 = 9.3094e-13 > 5e-13. Hand check of the end: T + q_T/2 = 0.10751025390625 + 4.967141088e-5 = 0.10755992531713. Hand check of the gain (evidence/damping.txt, old minus new height at equal times): 4.84e-3 at t = 0, 2.43e-3 at 0.05, 9.6e-4 at 0.1, 8.15e-4 at 0.1075 (factor 0.168 = e^-1.78; the band rates of the audit integrate to about 1.82). Near the end |Q'| = 2 + r q is about 2, so the path arrives 8.15e-4/2 = 4.07e-4 earlier; observed 4.05e-4.

## 3. The row-55 tube
- The checker requires, on row 55 only, 289/10000 <= Q <= 361/10000 for the whole ideal arc and its reset.
- TeX l.2392-2394: "This retained tube check is conservative. The old native improvement associated with that tube has already been replaced in the final source ledger and is not credited a second time."
- verification/source_overlays/ANALYTIC_SOURCE_JOIN_REVIEW.txt section 4: the old CONTACT_SOURCE differed from M8_LATE_SOURCE only by a row-55 native improvement of 5e-6; the final MAIN39 construction removes it; the direct consumer does not read CONTACT_SOURCE or the local-native certificate. verify_source_overlays_independent.py rebuilds final row 55 as M8 base + its global pure-source gain and records row55_local_native_credit_used = False; so does the credit audit.
- So no premise uses the tube; the kernel and own certificates of row 55 cover q in [0, cap]. Variant A turns the check off (one line, evidence/verifier_A.diff). Variant B keeps it by outward resets (q_out = max(ideal, floor) on 31 arcs; any q_out >= the ideal arc is allowed) and costs 1.45e-5.

## 4. Stage results
| Stage | Reads the flow? | This run (variant A swapped in) | s |
|---|---|---|---|
| native | no | PASS | 711 |
| activation | no | PASS | 5 |
| spectral | no | PASS | 12 |
| heads28 | no | PASS | 75 |
| heads39 | no | PASS | 44 |
| wall | no | PASS | 37 |
| kernel | no | PASS (all 67 rows on q in [0, cap]) | 865 |
| arcs | yes: new flow, old-run constants replaced, tube check off | PASS: 220,181 arcs, 5,777,345 exact checks, endpoint as above | 291 |
| affine | no | PASS | 2 |
| overlays | no | PASS | 25 |
| compact-arcs | yes: new transcript (arc count changed) | PASS: 3,963,258 checks | 56 |
| interfaces | reads arcs/ACTUAL_67_ROW_ARC_INPUTS.json, built from plan and source only; byte-identical to the author's (sha256 29902c66...) | PASS | 1 |
| physical-debit-map | no | PASS | 1 |
| arithmetic-credit-map | no | PASS | 1 |
| direct-source | no | PASS | 5 |

Variant B (tube kept): arcs PASS with the author's complete check set (220,211 arcs, 5,780,535 checks), compact-arcs PASS (3,963,798); the other 13 stages are the same as above. Receipt: evidence/replayA_RECEIPT.json; log tails: evidence/replayA_logs.txt.

## 5. Self-checks
- **Producer = author's path at 7/8.** make_flow_b2.py --b0 7/8: all 221,010 qin/qout identical to FLOW.json.gz (sha256 a664e829...), same T, q_T, endpoint 427693072077942473661578098491/3961408125713216879677197516800 and cap clearance 21059999/10^12. Only the last bits of the stored rate/source/own bounds differ (about 1e-28; both satisfy every check). evidence/repro78_compare.txt, FLOW_REPRO78_SUMMARY.json.
- **Independent float model** (evidence/float_check.py; own-pair jets via (p-y)^-(m+1) + (p+y)^-(m+1), not Newton sums): 7/8 -> 0.10796491007877 (exact 0.10796491007877157); b0 -> 0.10755992531713 (exact 0.10755992531713451); 0.8748867 -> 0.10793660221626 (Raymond, TeX l.210: 0.107936602217).
- **What the exact checks show** (A): smallest slacks: source rounding 5.0e-35, own rounding 1.2e-35, rate after eta 5.1e-30, Pade 2.2e-35; ideal-arc floor 9.93e-5 (last arc); cap clearance 1.05e-3 (row 64); rates 12.38-22.67.

## 6. Not re-optimised, and why
- Stop threshold 10^-4 -> 2*10^-6: model gain 2.8e-8. Continuous rates instead of frozen arcs: 3.5e-8. Not worth a new certificate.
- Probe mix and caps: the rows were tuned to the old path. A redesign for the lower path needs new kernel covers and the full source pipeline (native, activation, spectral, heads, overlays, direct-source), which the bundle checks but does not produce. By the envelope argument the gain is second order in the path shift except through the now slack caps (row 0's kernel is tightest at q = cap); rough size <= 1e-5.

## 7. Trust that remains
- Raymond's deposit is unrefereed and AI-assisted; our audit (Riemann/Deep/DBNAudit/REPORT.md) found it sound. This bound inherits all of its inputs except the 7/8 theorem.
- The half-plane Re s > 0.87338196661320752 is a Lean-kernel-checked theorem (Joseph Bebel, 31 Oct 2026; TREE.md N7); its sharper successor B_2 + 2.2e-254 (1 Nov) changes nothing here.
- The producer is ours; it is trusted for nothing, since the author's verifier checks every arc. Its 7/8 run reproduces the published path exactly.
- Novelty (5 Nov 2026, 14:26 UTC; evidence/novelty_check.txt): arXiv math.NT Oct and Nov 2026 listings and Zenodo's newest "de Bruijn-Newman" records show no bound for zeta below 0.107964910079 and no use of the sharper half-plane in Raymond's chain (newest: Planat-Raymond 0.156821, arXiv 2611.06245; Dirichlet 9/32, Zenodo 23931764).

## 8. Replay
```sh
EV=$PWD/Riemann/Deep/DBNB2/evidence; mkdir -p /tmp/dbnb2 && cd /tmp/dbnb2 && sh $EV/replay_dbnb2.sh
```
It downloads Proof_and_Replay.zip (201 MB, sha256 bf5c7b3d...), runs the 7/8 reproduction, builds flows A and B, patches and runs the author's arc and compact-witness checkers on both, runs the float model, and then the author's full 15-stage driver with the new arc stage swapped in. One core, about an hour, peak about 0.7 GB.

## Files
- evidence/make_flow_b2.py: the producer (exact FLINT rationals, the author's outward 2^-96 rounding of the path).
- evidence/patch_verifier.py, verifier_A.diff, verifier_B.diff: the only changes to the author's checkers.
- evidence/FLOW_B2_SUMMARY.json, FLOW_B2T_SUMMARY.json, FLOW_REPRO78_SUMMARY.json, repro78_compare.txt.
- evidence/A_*.json.gz, B_*.json.gz: the author's arc-verifier reports (all 67 row statistics, minima, maxima, transcript hashes), compact receipts and worst-arc witnesses.
- evidence/float_check.py, float_check.out, damping.txt; evidence/replayA_RECEIPT.json, replayA_logs.txt; evidence/A_verify_run.log, B_verify_run.log; evidence/replay_dbnb2.sh; evidence/novelty_check.txt.
- Flow sha256: A f8996d5e7af7807fbb1b8362d192fdf3a24f7d28e2a0709805d660026235cd02, B bc35d48806c41f3b7ba7c96a7ffe27e092360c781c11a1d05192257d6297aada (producer sha256 c3e8acd8...).

RH relevance: N2
