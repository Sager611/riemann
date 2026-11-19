COVERED: λ² ∈ [2, 227/6 = 37.8333] ∪ {40} independently certified (QW_λ has a simple, even ground state; chain stopped by the time limit at 227/6, no step failed; Botha's lower bounds μ remain imported on (37.8333, e⁴] \ {40}).

# BothaChain — our certificate of CCM's simple even ground state continued from λ² = 385/12 to 227/6 (bothachain, 14–19 Nov 2026)

## 1. Short answer
| λ² | lower bounds (λ₁ᵒᵈᵈ, λ₂ᵉᵛᵉⁿ) | upper bound (λ₁ᵉᵛᵉⁿ) | status |
|---|---|---|---|
| [2, 30] | ../WeilSimple [2, 20], ../CCMFMM [20, 30] | same | independently certified (before) |
| [30, 385/12 = 32.0833] | ../BothaAudit (CCMFMM's chain continued) | same | independently certified (before) |
| [385/12, 227/6] = [32.0833, 37.8333] | ours: 8 computed steps + 1 splice into CCMFMM's window 35 | ours (Lemma G) | **independently certified, this round** |
| 40 | CCMFMM | CCMFMM | independently certified (single window) |
| (37.8333, e⁴] \ {40} | Botha's Theorem V μ (imported) | ours (../BothaAudit point.py) | conditional on his μ only |
- In Botha's L = log λ: the new stretch is L ∈ [1.7341, 1.8167]; it contains 8 of his 31 windows in (1.7006, 2] (L = 1.73505 … 1.81525).
- Every step passed at its first attempt; every process ≤ 154 MB (memcap peaks ≤ 120 MB); every step re-checked from its records (table_chain.py).
- No contradiction with Botha's records or our float scan (sec. 5): our certified levels ε stay ≥ 28× below the rigorous Galerkin upper bounds at
  every window below each step end; his μ exceed our U(c_i) by ≥ 21.5 at the next window.
- Stopped by the time limit (and a 5-day stall of this session), not by a failure; e⁴ needs ≈ 24 more steps, ≈ 40–50 CPU hours on this route (sec. 6).

## 2. Setup (what was copied, what changed)
- ccm/: ../BothaAudit/ccm (= ../CCMFMM's scripts) copied; certificate code unchanged: gstate.py (U = Q(g)/‖g‖², Lemma G with one vector),
  prep_ccm.py, slt.c, weil_matrix.py, o1fix.py, rscan.c (md5 bd77d55f0a0bff7e7955e6c734aab9a9 = WeilZeta's).
- Driver changes only (ccm/chain_ccm.py header): one ok step per call (CCM_MAXSTEPS=1, so every call fits memcap 290 MB / 1800 s); the scan
  premise by ../WeilZeta/scan_zeta.py + rscan2.c (proved early exit, ../WeilZeta/RIGOR.md; CCMFMM's own route at c = 35, 40) before prep, written
  to the same cache file prep reads; a failed (c_i, c_{i+1}) pair is not retried; Lemma G's non-degeneracy checked on every U (finite Arb ball,
  mid > 0, radius < 10⁻³·mid, so ‖g‖² > 0 and g ≠ 0; record field lemmaG); run_ccm.py caches a finished LDL run (resume after a kill).
- Step heuristics (logs/chain_int.jsonl, mkint.py): bothaaudit's rows c = 30..35 unchanged; rows 36..56: o1 := 0.6 × the float-scan odd value at
  Botha's windows (log-interpolated), bits := bitsrule of the scan e1. No rigour (they place c_{i+1} and the precision only).
- Validation before the run: (i) a chain step 20 → 165/8 from CCMFMM's U(20) with the new driver: ok, ν = 3.314883e-94 (defl), 3.314883e-96 (odd),
  n = 987/986, lemmaG true, 56 s, 54 MB; (ii) rscan2 against CCMFMM's rscan caches at c = 715/24 (K2 = 1481) and 30 (K2 = 1501): identical
  k-ranges and identical selected odd/even node sets (1129/537, 719/1029), 7 s each instead of minutes.
- Binaries: cc -O2 slt.c -lflint -lgmp -lmpfr -lm (Homebrew FLINT 3), cc -O2 -ffp-contract=off rscan.c, rscan2.c; python-flint 0.9.0.

## 3. Certificate logic (unchanged) and the two driver economies
- One step [c_i, c_{i+1}] (FULLPROOF sec. 5 of ../BothaAudit; ../CCMFMM/DERIVATION.md): U(c_i) = Q(g)/‖g‖² ≥ λ₁ᵉᵛᵉⁿ(c_i) (Lemma G, one exact vector, NG modes);
  at c_{i+1} the odd certificate (O) and the deflated even certificate (D) at the level ε_i = 2.00001·U(c_i): a structured LDL at the shift ν = ε_i/10
  with the expected inertia and ν > κη gives λ₁ᵒᵈᵈ(c_{i+1}) ≥ ε_i and λ₂ᵉᵛᵉⁿ(c_{i+1}) ≥ ε_i. Monotonicity (every λ_k non-increasing in c) then gives
  λ₁ᵉᵛᵉⁿ(c) ≤ U(c_i) < ε_i ≤ min(λ₁ᵒᵈᵈ, λ₂ᵉᵛᵉⁿ)(c) on [c_i, c_{i+1}]: simple, even ground state.
- Pairs (ccm/chain_pair.py, from 67/2 on): (D) is not the binding sector (e₂/e₁ ≈ 2–3·10⁹ against o₁/e₁ ≈ 6–7·10⁴), so one (D) at c_c serves two steps
  c_a → c_b → c_c: (O) at c_b at ε_a, (O) at c_c at ε_b, (D) at c_c at ε_a; λ₂ᵉᵛᵉⁿ(c) ≥ λ₂ᵉᵛᵉⁿ(c_c) ≥ ε_a > ε_b for all c ≤ c_c. Same certificates, same
  monotonicity; tried only if 0.6 × the scan's U₂ at c_c exceeds 7.1·ε_a, else the standard (D) at both points. Saves one (D) (≈ 40 % of a step) per pair.
  Tested at c ≈ 20–22 from CCMFMM's U(20): both paths (a failed second step closed alone with its own (D); a shared (D) at 67/3 at level 8.475e-102)
  certify; U(62/3) = 3.87093071196e-100 equals CCMFMM's record. table_chain.py re-checks every record: (O) level = ε_i at c_{i+1}; (D) level > U(c_i)
  at a point ≥ c_{i+1}; ν > κη in both; Lemma G check; U(c_i) < ε_i.
- Splice (splice.py): a chain point c_i < c_w with U(c_i) < CCMFMM's certified level at the window c_w ∈ {35, 40} closes [c_i, c_w] with CCMFMM's (O)
  and (D) at c_w (ε := its level); the chain resumes from c_w with CCMFMM's U(c_w).
- Incident (20:52): chain.sh was edited while a call was still reading it; that call's tail failed after its step had been recorded (record intact).
- Precision (from c = 35 on): bits := 64⌈(1.25·log2(1/e₁) + 64)/64⌉ instead of CCMFMM's factor 1.6 (rows ≥ 36 of logs/chain_int.jsonl: 832 bits at c ≤ 36,
  896 at 37–39, 960 at 40–41, instead of 1088/1152). Rigour does not depend on it (Arb balls; a certificate passes only if ν > κη). Validated at the
  overlap point c = 35 against CCMFMM's certified window (val35.sh, logs/val35.jsonl): rscan2's node sets = CCMFMM's cache (664 odd, 1531 even);
  U(35) at 832 bits = [3.35474138011541e-177 ± 1.10e-192], identical to CCMFMM's 1024-bit value in every printed digit; the odd certificate at CCMFMM's
  level 8.972773e-175 passes at 832 bits with ν = 8.972773e-176, κη = 2.2e-211 (κ = 3.4e20; at 1024 bits κη was 3.6e-269), n = 5499 as there;
  peak 92 MB; 446 s on a machine at load 15–19 (slt 330 s wall, 122 s CPU).

## 4. The chain (logs/chain.jsonl; table_chain.py → logs/chain_table.md)
Rows 1–3 bothaaudit's, 4–12 this round (8 computed steps, row 8 the splice). ε/U(c_i) = 2.000 by the rule (95.2 at the splice); ν = ε/10 is the LDL
margin and κη its proved error bound; n = particles of (D)/(O); for a shared (D) the "by" column names its point and level.
| # | c_i → c_{i+1} | U(c_i) ≥ λ₁ᵉᵛᵉⁿ on the step | ε: λ₁ᵒᵈᵈ, λ₂ᵉᵛᵉⁿ(c_{i+1}) ≥ ε | ε/U(c_i) | certified ν (defl, odd); max κη | U(c_{i+1}) | NG; n (defl/odd) | bits | peak RSS MB (g; prep, LDL defl; prep, LDL odd; scan) | s | by |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 30 → 92/3 | 3.27353749418e-150 | 6.547108e-150 | 2.000 | 6.547108e-149, 6.547108e-151; 6.7e-234 | 6.75090730136e-154 | 320; 3332/3331 | 896 | 111; [83, 16]; [72, 15]; None | 361.6 | bothaaudit |
| 2 | 92/3 → 251/8 | 6.75090730136e-154 | 1.350188e-153 | 2.000 | 1.350188e-154, 1.350188e-154; 1.2e-251 | 1.00300264109e-157 | 320; 4450/4449 | 960 | 115; [None, 34]; [100, 32]; None | 611.1 | bothaaudit |
| 3 | 251/8 → 385/12 | 1.00300264109e-157 | 2.006015e-157 | 2.000 | 2.006015e-158, 2.006015e-158; 5.1e-254 | 1.72477599362e-161 | 320; 4763/4762 | 960 | 106; [99, 31]; [95, 28]; None | 975.8 | bothaaudit |
| 4 | 385/12 → 787/24 | 1.72477599362e-161 | 3.449569e-161 | 2.000 | 3.449569e-162, 3.449569e-162; 1.1e-250 | 2.40770250625e-165 | 320; 4815/4814 | 960 | 115; [105, 36]; [102, 34]; 69 | 836.4 | bothachain |
| 5 | 787/24 → 67/2 | 2.40770250625e-165 | 4.815429e-165 | 2.000 | 4.815429e-166, 4.815429e-166; 2.7e-271 | 3.83864482386e-169 | 320; 4884/4883 | 1024 | 123; [105, 36]; [102, 33]; 71 | 919.7 | bothachain |
| 6 | 67/2 → 821/24 | 3.83864482386e-169 | 7.677328e-169 | 2.000 | 7.677328e-170, 7.677328e-170; 3.5e-272 | 6.11584390297e-173 | 320; 5073/4940 | 1024 | 124; [107, 37]; [95, 33]; 71 | 1587.2 | bothachain (defl shared: at 419/12, level 7.677328e-169) |
| 7 | 821/24 → 419/12 | 6.11584390297e-173 | 1.223175e-172 | 2.000 | 7.677328e-170, 1.223175e-173; 7.4e-271 | 9.42783071633e-177 | 320; 5073/5072 | 1024 | 115; [107, 37]; [103, 35]; 71 | 1131.2 | bothachain |
| 8 | 419/12 → 35 | 9.42783071633e-177 | 8.972773e-175 | 95.173 | 8.972773e-174, 8.972773e-174; 3.6e-269 | 3.35474138012e-177 | 320; 5500/5499 | 1024 | 120; [84, 21]; [81, 19]; None | 0 | splice: CCMFMM window |
| 9 | 35 → 857/24 | 3.35474138012e-177 | 6.709516e-177 | 2.000 | 6.709516e-178, 6.709516e-178; 4.7e-212 | 4.54260306834e-181 | 400; 5097/5036 | 832 | 150; [97, 28]; [100, 29]; 71 | 1329.2 | bothachain (defl shared: at 437/12, level 6.709516e-177) |
| 10 | 857/24 → 437/12 | 4.54260306834e-181 | 9.085252e-181 | 2.000 | 6.709516e-178, 9.085252e-182; 5.2e-232 | 6.47832192481e-185 | 400; 5097/5096 | 896 | 154; [97, 28]; [94, 26]; 71 | 731.4 | bothachain |
| 11 | 437/12 → 297/8 | 6.47832192481e-185 | 1.295671e-184 | 2.000 | 1.295671e-185, 1.295671e-185; 1.6e-231 | 9.75519723648e-189 | 400; 7046/6854 | 896 | 154; [None, 44]; [None, 34]; None | 586.8 | bothachain (defl shared: at 227/6, level 1.295671e-184) |
| 12 | 297/8 → 227/6 | 9.75519723648e-189 | 1.951049e-188 | 2.000 | 1.295671e-185, 1.951049e-189; 4.8e-231 | 1.51634909623e-192 | 400; 7046/7045 | 896 | 154; [None, 44]; [None, 42]; None | 586.8 | bothachain |

covered by the chain: [30, 227/6] = [30, 37.833333]; ok steps 12 (bothaaudit 3), failed attempts 0; max RSS (our runs) 154 MB; re-check failures: none

## 5. Against Botha's μ and the float scan, every 5 steps (compare.py → logs/compare.out; run after the splice and at the end)
Uo, U2 = our rigorous Galerkin upper bounds (../BothaAudit point.py) at his windows c_B; his μ from his records (../BothaAudit/logs/compare2.out).
Required by monotonicity: ε_i ≤ Uo(c_B), U2(c_B) at every window c_B ≤ c_{i+1} (a ratio < 1 would refute a certificate): all hold, min 6.3 (bothaaudit's
step 2), 28.2 in this round (step 297/8 → 227/6), U2 ratios ≥ 1.25e6. μ(c_B′)/U(c_i) ≥ 21.5 in this round: his levels and our upper bounds keep the
same order. Float scan at the windows: o₁/e₁ = 5.9–8.3·10⁴, e₂/e₁ = 1.9–3.7·10⁹ along the new stretch, no crossing; all steps passed (no lead).
| step c_i → c_{i+1} | eps_i (our lower bound at c_{i+1}) | U(c_i) | Botha window c_B ≤ c_{i+1} (nearest) | Uo(c_B)/eps, U2(c_B)/eps (≥ 1 required) | next window c_B′ ≥ c_i | mu_odd(c_B′), mu_l2(c_B′) | min mu(c_B′)/U(c_i) | scan Uo/U, U2/U at c_B |
|---|---|---|---|---|---|---|---|---|
| 30 → 92/3 | 6.547108e-150 | 3.27353749418e-150 | 30.570758 (L = 1.71002197265625) | 19, 5.43e+05 | 30.570758 | 5.649e-149, 7.393e-145 | 17.3 | 5.52e+04, 1.58e+09 |
| 92/3 → 251/8 | 1.350188e-153 | 6.75090730136e-154 | 31.345426 (L = 1.7225341796875) | 6.3, 2.08e+05 | 31.345426 | 3.692e-153, 5.61e-149 | 5.47 | 5.95e+04, 1.96e+09 |
| 251/8 → 385/12 | 2.006015e-157 | 1.00300264109e-157 | 31.345426 (L = 1.7225341796875) | 4.24e+04, 1.4e+09 | 32.139724 | 2.316e-157, 3.363e-153 | 2.31 | 5.95e+04, 1.96e+09 |
| 385/12 → 787/24 | 3.449569e-161 | 1.72477599362e-161 | 32.139724 (L = 1.73504638671875) | 1.49e+04, 4.72e+08 | 32.139724 | 2.316e-157, 3.363e-153 | 1.34e+04 | 6.03e+04, 1.91e+09 |
| 787/24 → 67/2 | 4.815429e-165 | 2.40770250625e-165 | 32.926002 (L = 1.74713134765625) | 6.07e+03, 2e+08 | 32.926002 | 1.319e-161, 1.883e-157 | 5.48e+03 | 5.91e+04, 1.95e+09 |
| 67/2 → 821/24 | 7.677328e-169 | 3.83864482386e-169 | 33.719166 (L = 1.759033203125) | 2.3e+03, 7.53e+07 | 33.719166 | 7.709e-166, 1.164e-161 | 2.01e+03 | 6.58e+04, 2.15e+09 |
| 821/24 → 419/12 | 1.223175e-172 | 6.11584390297e-173 | 34.518793 (L = 1.770751953125) | 783, 2.81e+07 | 34.518793 | 3.943e-170, 6.768e-166 | 645 | 7.24e+04, 2.59e+09 |
| 419/12 → 35 | 8.972773e-175 | 9.42783071633e-177 | 34.518793 (L = 1.770751953125) | 1.07e+05, 3.82e+09 | 35.315821 | 1.915e-174, 3.517e-170 | 203 | 7.24e+04, 2.59e+09 |
| 35 → 857/24 | 6.709516e-177 | 3.35474138012e-177 | 35.315821 (L = 1.78216552734375) | 706, 2.84e+07 | 35.315821 | 1.915e-174, 3.517e-170 | 571 | 7.36e+04, 2.96e+09 |
| 857/24 → 437/12 | 9.085252e-181 | 4.54260306834e-181 | 36.113614 (L = 1.7933349609375) | 269, 1.18e+07 | 36.113614 | 9.545e-179, 1.943e-174 | 210 | 7.7e+04, 3.38e+09 |
| 437/12 → 297/8 | 1.295671e-184 | 6.47832192481e-185 | 36.920414 (L = 1.80438232421875) | 88.6, 3.86e+06 | 36.920414 | 4.545e-183, 8.955e-179 | 70.2 | 7.93e+04, 3.46e+09 |
| 297/8 → 227/6 | 1.951049e-188 | 9.75519723648e-189 | 37.731419 (L = 1.81524658203125) | 28.2, 1.25e+06 | 37.731419 | 2.097e-187, 4.231e-183 | 21.5 | 8.31e+04, 3.67e+09 |

violations: 0 (windows checked: all c_B <= c_{i+1} for each of 12 steps)

## 6. Cost, memory and what e⁴ needs (cost.py; one process at a time, background QoS, machine load 4–19 on 10 cores)
- Per step (CPU s, wall s; peak RSS): 787/24: 671, 836; 67/2: 672, 920 (960–1024 bits, standard steps, ≈ 0.7 in c each); pair 67/2 → 821/24 → 419/12:
  1101 CPU s, 1587 wall s for 1.42 in c; pair 35 → 857/24 → 437/12 at 832/896 bits: 939 CPU s, 1329 wall s. Pairs + the precision rule cut the CPU
  per unit of c by ≈ 30 % although n grows (4814 → 5096 particles, G ≈ 0.9–1.3·10⁵ sources). Max RSS 154 MB (gstate at NG = 400); prep ≤ 112 MB,
  LDL ≤ 37 MB, rscan2 ≤ 75 MB; memcap peaks ≤ 120 MB.
- The next pair 437/12 → 297/8 → 227/6 needed > 1800 s wall at load ≈ 14 (killed by memcap at 22:27 UTC, 14 Nov, during its shared (D) LDL; every
  other piece cached); it was resumed on 19 Nov after a 5-day stall of this session (sec. 4).
- Step length is fixed by the method, not by cost: ε_i = 2.00001·U(c_i) must stay below λ₁ᵒᵈᵈ(c_{i+1}), and o₁/e₁ ≈ 6–7·10⁴ against a decay
  of ≈ 10^5.4 per unit of c allows ≈ 0.71 per step (all 8 computed steps passed at their first attempt: the 0.6 × scan heuristic is conservative).
- Outlook to e⁴ = 54.59815 from 227/6: ≈ 24 more steps (12 pairs). The cost is the Cauchy sum n·G (≈ 5·10⁸ at c ≈ 36; ≈ 7·10⁸ at 40, CCMFMM;
  ≈ 5·10⁹ at 47 and ≈ 10¹⁰ at 53, ExploreFMM's counts) times (bits/64)²: ≈ 0.3 CPU h per pair at c ≈ 36, ≈ 3 at 47, ≈ 6–7 at 53, in total
  ≈ 40–50 CPU h on one core (wall ≈ 1.2–2.7× that here). Memory stays under the cap on this route as far as we can predict: gstate ≈ 157 MB at
  832 bits grows ∝ bits to ≈ 240 MB at 1280 bits (c ≈ 54.6; with CCMFMM's factor 1.6, 1664 bits, it would exceed 290 MB, so the precision rule
  matters there), the scan at c ≈ 53 is ≈ 10¹¹ nodes (≈ 30 min with rscan2, chunked); prep's growth with G (3–4× at c ≥ 47) is unmeasured: check first.
  The wall is time, not memory; nothing failed.

## 7. Reviewer (one reviewer subagent, 20:40 UTC; values pasted inline, no code run)
Step 385/12 → 787/24 (U(c_i) = [1.72477599361988e-161 ± 1.83e-176], ε = 3.449569e-161, (D) and (O) at 787/24 with ν = 3.449569e-162, κη = 4.8e-261 /
1.1e-250, inertia n_U + 1 / n_U, U(787/24) = [2.40770250625287e-165 ± 6.59e-182]) and the splice rule at c = 35: **CORRECT**. Its notes: monotonicity
directions right at both ends; the certified level is ε, ν only the margin; 2.00001·U_up = 3.44956923…e-161 rounds down to ε, harmless since ε > U_up;
finite Arb division excludes ‖g‖² = 0 (exactness of g and the 321-mode span hold by construction: gstate rounds g to exact dyadics); the next level
2.00001·U_up(c_{i+1}) = 4.815e-165 < ε; the splice is valid under U_up(c_i) < L, the lower bounds transported backward from 35, U(35) < L starting the next step.
Second reviewer (19 Nov; the pair logic, which closes 6 of the 8 computed steps): pair 437/12 → 297/8 → 227/6 with (O) at 297/8 (level ε_a =
1.295671e-184), (O) at 227/6 (ε_b = 1.951049e-188), one (D) at 227/6 at ε_a: **CORRECT**; "the shared higher-level (D) suffices for both intervals;
no (D) at c_b is needed"; note: printed upper ends are 12-digit roundings (9.75519723648e-189 vs 9.755197236484…e-189), still far below ε_b.

## 8. Conclusion, trust, novelty
- Independently certified (our lower and upper bounds only): QW_λ has a simple, even ground state for every λ² ∈ [2, 227/6 = 37.8333] and at λ² = 40.
- Still imported trust on (37.8333, e⁴] \ {40}: Botha's lower bounds μ (certify.py, sr.c; PROOFS P1–P23), Elston's stage-S input s(L) and integer
  Cholesky, the operator setup (../BothaAudit sec. 8); our rigorous upper bounds already replace his fr_sweep there.
- Trust base of our part, as in ../CCMFMM sec. 9: Arb (python-flint 0.9.0; FLINT 3 in slt.c), IEEE binary64 in rscan2.c (../WeilZeta/RIGOR.md),
  our Python and C (no proof assistant). This round changed only drivers (step bookkeeping, pairs, splice, scan route, precision choice); no certificate
  code. Printed U_up is the 12-digit rounding of the ball's upper end (relative error ≤ 5·10⁻¹³; second reviewer's note), far inside the factor 2.00001.
- Novelty: none claimed. Botha (Zenodo 23852617) claims all of λ² ∈ (1, e⁴]; ACLG (arXiv 2611.13129) report the same ordering numerically up to
  λ² = 300. An audit result: our independent coverage of his imported-trust Piece IV is now L ∈ (1.4035, 1.8167] and L = 1.8444. No DISCOVERIES entry.
- Suggested TREE N13 clause: "… and for λ² ∈ [385/12, 227/6 = 37.8333] (bothachain, 19 Nov 2026: the chain continued in 9 steps, one deflated
  certificate per pair of steps, spliced at λ² = 35, ≤ 154 MB; independent coverage λ² ∈ [2, 37.8333] ∪ {40}; Botha's μ imported only on
  (37.8333, e⁴] \ {40}; e⁴ ≈ 24 more steps, ≈ 40–50 CPU h; [BothaChain](Riemann/Deep/BothaChain/REPORT.md))".

## 9. Files (Riemann/Deep/BothaChain/)
- ccm/: certificate scripts (gstate.py, prep_ccm.py, slt.c, weil_matrix.py, rscan.c, o1fix.py: unchanged copies), run_ccm.py (+ LDL-run cache),
  chain_ccm.py (driver changes in its header), chain_pair.py (pairs), scan_zeta.py and rscan2.c (copies of ../WeilZeta's).
- chain.sh (one memcap'd call), loop.sh, loop2.sh (call loops with stop time), splice.py, mkint.py (step heuristics), val35.py / val35.sh (precision
  validation), table_chain.py, compare.py, cost.py.
- logs/: chain.jsonl (all records, bothaaudit's 3 first), chain_table.md, compare.out, cost.md, val35.jsonl, chain_int.jsonl, chain.out, chain.err;
  runs/ (every prep and LDL JSON line and /usr/bin/time record, gstate outputs and vectors), scan/ (rigorous scan caches of the chain points).
- Replay (binaries as in sec. 2, into .scratch/bothachain/ccm): PAIR=1 bash chain.sh 55 1800 continues from the last record (U(227/6) = 1.51634909623e-192).

RH relevance: N13
