23752103: SOUND; 23272725: SOUND - no false or unjustified step found; both author replays pass on new hardware, and every independent check (kernels, all 221,010 arcs, finite heads, clock credits, a direct physical test at t=0) agrees, so Raymond's Lambda <= 0.107964910079 is the N2 record (unrefereed, computer-assisted, resting on the 7/8 theorem).

# Audit of D. J. Raymond's Zenodo bounds Lambda <= 0.107964910079 (23752103) and Lambda <= 0.116397726408291 (23272725)
Orca agent dbnaudit, 2026-11-05, 11:14-13:45 UTC. Scripts and outputs: evidence/ (144 KB).

## Short answer
- **23752103 is sound as written.** It is a 63-page analytic argument plus 15 finite certificate stages. I re-derived every analytic step on the path from the 7/8 theorem to the terminal strip (table in section 2). The fifteen-stage replay passes on macOS arm64 with a separately installed arithmetic library. My own exact recheck of all 221,010 comparison arcs passes. Kernel, head and physical spot checks agree.
- **The gain from 0.1164 to 0.1080 comes almost entirely (98.7%) from a new kind of probe functional.** It is a signed combination of up to 26 translated probes at heights 1.36 to 6. Its kernel inequality is proved only on the boundary b = y and then extended to every other-zero height by a positive harmonic (Poisson) extension. Arithmetic credits give 1.3%. Derivative ("jet") rows, the finer time grid, the seed and the stopping rule together give less than 0.02%.
- **The most delicate ingredient is load-bearing.** After the common-clock and Gram relation credits, the charged finite-head cost is 0.0009-0.0108 below the head's unconstrained phase supremum (rows tested). The head estimate is therefore valid only because the relation constraints exclude those phase vectors. These credits use exact prime relations such as 3^4*5*11*19^2*23*47*71 / (2^4*7^5*13*17*31*67) = 1 + 8.1e-12. They hold because near X these relations pin combinations of prime phases near odd multiples of pi; far from X the growth of the main term pays instead. I re-derived the argument and could not falsify it numerically (section 4). Even if every arithmetic credit were dropped, the model bound would be 0.108077.
- **23272725 is sound under 7/8.** Its 44 new kernels, self caps, seed 9/16 + 1e-6 and the 8 stored histories all check. It is superseded by 23752103.

## 0. Versions audited
- **23752103** v1.0 (29 Oct 2026), the only version of concept 23752102.
  - The GitHub companion Dale-Dug/de-bruijn-newman-0107965 has head 9d1a340 (the cited commit) and mirrors the deposit byte for byte.
  - The other newer Raymond de Bruijn-Newman items claim no smaller bound: 23781981 (30 Oct) is v1.1 of the 0.149396 paper, and Planat-Raymond 0.156821 (4 Nov; arXiv 2611.06245) uses the classical strip.
  - sha256 of what I audited:
    - Unconditional_0107965.tex adbd147a63e2c556a3c5c332e78c011bf3420071a32678262d61aa3aac7ee4f4
    - Unconditional_0107965.pdf df40d365c71e09973342ec30eb20cdebdf17faf2d89de32bc7a04772cd1e1e2e
    - Proof_and_Replay.zip bf5c7b3dca9e17e84e146ecc5ca0c6f671b409613088212e3f0b8a0b15f6aaf6
    - Verify_Deposit.py 11a129ecffad544397eb8d576f94cb7f97fce06495ad8faa2ec485d831568f60
  - All 16 files match the Zenodo md5 values and SHA256SUMS. Verify_Deposit.py reports PASS.
  - Input identities: plan ec346b7c...; source b59905d7...; flow a664e829... (FLOW.json.gz).
- **23272725** (10 Oct 2026), the only version of concept 23272724.
  - Conditional_0116398.tex 20226820c128f11e677820b83995c3609e024b207ce26e0c07ee22818e7ccd00
  - Conditional_0116398.pdf 6eec50fdb28f31211da2cfb9d5d3640b64723ec6c044ecc0db9774a1fe2e1b18
  - Proof_and_Replay.zip 045bd7703019a78167e746c8260c0135803b80aa072d5ec0313d5a2b4cc85b89 (md5 matches Zenodo)

## 1. Gain decomposition (question 1)
**Their own tables.** 23752103 has no table that splits 0.1164 into 0.1080; it only states the difference 0.008432816329291. It does give per-row data:
- the post-M8 "pure-source gain" column in tables/ALL67_SOURCE_DECIMAL_GUIDE.csv;
- the credit split in DIRECT_SOURCE_RECONSTRUCTION.json.

23272725 does have its own table (0.118068908899880 -> 0.117794794642445 new kernels -> 0.116787767525899 retuned heights -> 0.116397726408291 matched-log source).

**My decomposition** (evidence/decompose.py) uses a float model of their exact scheme: arcs of 1/2048000, rate 2(S_i(t0) - O_i(q0)) - 1e-12, stop at q = 1e-4, bound t + q/2. The model reproduces their endpoint to 1.5e-10 (0.107964910225). Each counterfactual removes one ingredient, using their own fields:

| Ingredient removed | Model bound | Share of the 0.0084328 gain |
|---|---|---|
| none (model of 23752103) | 0.1079649102 | - |
| all arithmetic credits (composite, common clock, fixed triads, Gram orbit) | 0.1080770442 | 1.33% (1.121e-4) |
| ... only the 39 post-M8 clock gains | 0.1079755391 | 0.13% (1.06e-5) |
| derivative (jet) components of rows 0, 16, 44 | 0.1079652722 | 0.004% (3.6e-7) |
| fine arcs (1/2048000 -> the old 1/128000) | 0.1079654341 | 0.006% (5.2e-7) |
| seed 9/16+1e-12 -> old 9/16+1e-6 | 0.1079649935 | 0.001% (8e-8) |
| stop at 1e-4 + cleanup -> old stop at 1e-8 | 0.1079648819 | -0.0003% (-2.8e-8) |
| all of the above together | 0.1080779866 | 1.34% |
| **remainder: the new functional class and its source** | 0.1164 -> 0.10808 | **98.66% (0.008320)** |

**Where the time is saved** (old = 23272725 main history, new = 23752103 model):

| q band | old rate | new rate | old dt | new dt | saved |
|---|---|---|---|---|---|
| 0.5625-0.5 | 11.30 | 12.56 | 0.007808 | 0.007212 | 0.000596 |
| 0.5-0.4 | 11.79 | 13.22 | 0.013717 | 0.012614 | 0.001103 |
| 0.4-0.3 | 12.47 | 14.16 | 0.015752 | 0.014393 | 0.001359 |
| 0.3-0.2 | 13.36 | 15.42 | 0.018801 | 0.017123 | 0.001678 |
| 0.2-0.1 | 14.74 | 17.38 | 0.023988 | 0.021965 | 0.002023 |
| 0.1-0.05 | 16.27 | 19.40 | 0.015611 | 0.014577 | 0.001034 |
| 0.05-0.01 | 17.45-18.15 | 20.86-21.81 | 0.015938 | 0.015341 | 0.000597 |
| 0.01-1e-4 | 18.61-18.85 | 22.39-22.68 | 0.004732 | 0.004690 | 0.000041 |

The new rates are 11-20% higher at every height.

**Why this beats DBN6's derivative-lever models (0.112-0.115).** DBN6 kept one probe height per cell and a scalar kernel constant (lambda <= 2) with a fixed shift. It added a derivative term with a uniform bound G.

23752103 instead solves, row by row, a design problem over signed weights of shifted probes at several heights. Row 0, for example, has 9 nodes at heights 2.5-6 with weights from -14.3 to +16.9. The only kernel constraint is the boundary inequality R_i(d, sqrt q) <= 4/(d^2+4q) (section 3.8). Lemma fnd:poisson then extends it to every 0 <= b <= y.
- Early on, the signed main mass is close to the ceiling 1/y that this constraint allows: M = 1.261 against 1/y = 1.333 at q = 9/16.
- With signed weights, the effective height sum(w h)/sum(w) = 1.59 can lie below the lowest actual height, 2.5.
- The probes sit at sigma_eff = (1+h)/2 + tL/8 > 5/4 (guard in Prop src:physical). The Dirichlet series therefore converge absolutely, and the Riemann-Siegel errors are about e^{-15} or smaller.
- The derivative lever itself is almost worthless here (3.6e-7).

## 2. Proof audit of 23752103 (question 2)
Every step below was re-derived from the TeX. **First false or unjustified step: none found.**

| Step | Verdict | Note |
|---|---|---|
| 7/8 input -> q <= 9/16 at t=0; seed 9/16+1e-12 [2.1, 4.4] | OK | (1+iz)/2 map; functional equation gives Re rho in [1/8,7/8]. 7/8 is a Lean theorem (our QRH port). |
| Startup: tau=1e-14, gap >= 1e-12 - r0 q0 tau = 9.30e-13 [4.4] | OK | exact; de Bruijn contraction on [0,tau]. |
| Wall at X=6000000185827, t<=41/250 [2.2] | OK | same wall as 22867547/23272725 (REPORT3/4); replay stage "wall" recomputes all 2,763,952 prefixes at 512 bits. |
| Finite head: zeros with abs(x) <= X real, simple for t>0 [Prop fnd:finite-head] | OK | Platt-Trudgian ordinate 3000175332800 > (X+1)/2 by 175239886; Rouche + forward Hermite splitting. |
| Kernel: boundary inequality R(d,sqrt q) <= 4/(d^2+4q) => K_y(d,b) >= (y/2)R(d,b) for 0<=b<=y [Lemma fnd:poisson] | OK | re-derived: R and K_y harmonic in abs(b) < y, even in b; strip Poisson multiplier cosh(b xi)/cosh(y xi) is positive; FT of K_y/y = (2pi/y)exp(-y abs(xi))cosh(b xi) = extension of boundary data 2/(d^2+4y^2) + (pi/y)delta_0. Signed nodes are fine because the boundary inequality is for the whole signed sum. |
| Current inequality q' <= -2 - 2q(J - O_pair(q)) [Prop fnd:current] | OK | Polymath Prop 3.1(i); own conjugate gives -2; J = O_pair + sum_other R (absolute convergence, Lemma fnd:convergence). |
| Compactness for t in [tau,T] [Lemma fnd:compact] | OK | Polymath Thm 1.5(i) (x >= e^{C/t} => real, one absolute C for 0<t<=1/2); only existence of C is used, so ineffectivity is harmless. |
| No repeated nonreal first contact [Lemma fnd:multiple] | OK | backward heat scaling -> Hermite polynomial He_m(w/(i sqrt2)) has a root i sqrt2 lambda, lambda>0: height excess ~ sqrt(eps) vs barrier change O(eps). Works at complex zeros. |
| Global first contact [Thm fnd:comparison, Prop flow:global-terminal] | OK | attained first contact (compactness + Rouche), at x>X by the finite head, all other zeros below; upward resets harmless; contact at arc endpoints handled from the left. Zeros entering from abs(x) <= X impossible (wall); from large x impossible for t>=tau (compactness). |
| Arc lemma: S_i nondecreasing in t (A_i>=B_i), O_pair nondecreasing in q => q' <= -2-(r+eta)q < Q' [Lemma flow:strict-arc] | OK | positivity of the ideal arc from e^{-ru} >= 1-ru. |
| [6/6] Pade: e^{-x} <= P(-x)/P(x) [4.3] | OK | coefficients are the diagonal Pade numerator; remainder -(x^13/12!)int e^{-tx}t^6(1-t)^6 <= 0. |
| Jet own pair (1/m!)d^m/dp^m 2p/(p^2-q) = (-1)^m S_{m+1}, Newton recurrence [4.3] | OK | 2p/(p^2-q) = 1/(p-y)+1/(p+y). |
| Terminal step: all zeros of H_T in abs(Im z) <= sqrt(q_T), then de Bruijn => Lambda <= T + q_T/2 [Lemma fnd:cleanup, Cor flow:endpoint] | OK | the comparison bounds every zero (head zeros are real); de Bruijn Thm 13 / Polymath Thm 3.2. |
| Source: common model, Riemann-Siegel errors on Re s in [-3,4], Gaussian 11/12 [Prop src:physical] | OK | Polymath Prop 6.1 holds for every real sigma; remainder rederived as in HB (audited in HBDBN); 4(u^2+11/12)-(u^2+6 abs(u)+2/3) = 3(abs(u)-1)^2. Errors ~e^{-S}, S >= 15 at every native height, so tiny. |
| Matched log residual abs(V-QT_K) <= B-(1-A)D_K-2Z_K [3.3] | OK | coefficient identity (log n)c_t(n) = sum over k dividing n of d~_t(k) c_t(n/k), all terms >= 0. |
| One positive spectral measure; envelope abs(F_i(lambda)) <= sum c_l e^{-h_l lambda/2} [3.5] | OK | one Laplace transform for all translations (frozen common model); Laplace transform of convolution = product. |
| Future cutoffs: activation integrals, three-exponential closed forms, unique maximum [App A] | OK | re-derived the closed forms A(s), D(s) and the derivative sign arguments; sums <= integrals (first branch decreasing, reflected branch shift paid by rho). |
| Affine native bound N(L) <= N(0)+beta M(L-L0)/4 [3.4] | OK | convexity near + global cap far. |
| Composite credit: 1-Re prod z_p^{e_p} <= Omega sum e_p(1-Re z_p) [App D, 39 rows] | OK | = abs(1-prod)^2/2 <= (sum abs(1-z))^2/2 <= Omega sum abs(1-z)^2/2 (Cauchy-Schwarz). Global in L by radial contraction + polydisc maximum principle. |
| Common clock: k.theta/pi = 2m+1+delta, delta = a(1+v+q)-h, sign duals with conjugate max{Psi(u), pi u-2} [3.6, App D] | OK | tau = pi(2e^L+t/8), e^L = N0^2(1+v); weak duality; F concave (affine + (A/4)log(1+v)) so corners suffice; far side paid by main slope A Delta >= d. |
| 28 shared rows: triad deficits, coordinate reserve mu kappa0, Gram LDL, cosine envelope [App D] | OK | triad identity (C/2) abs(sqrt(A/B)u+sqrt(B/A)v^bar - sqrt(AB)/C)^2; U-UWU = LDL^T => V^T W V <= I => sum r_p theta_p^2 >= sum w_j (s_j rho_j)^2; (1-cos v)/v^2 decreasing; c-c0 >= (1/4)(1-4 pi N0^2/x) = 2.24e-7 > kappa0. |
| Jet rows: Cauchy bound on one holomorphic log, analytic envelope sup abs(P_b(-p lambda/2)) e^{-a lambda} [3.7] | OK | abs(g^(r+1)(p)) <= (r+1)! eta / r^(r+1). |


**Known failure points tested**
- *Floors uniform in x >= X, every future cutoff, unbounded frequencies:* covered.
  - Activation integrals over all s >= 0 (App. A).
  - The affine native bound.
  - Residual terms decrease in L term by term.
  - The composite credit, globalised by radial contraction and the polydisc maximum principle.
  - The clock credit covers the near interval L0 <= L <= L0 + omega (omega <= 0.1); the main slope A*Delta >= d covers the far side.
  - The Gram credit covers the band X0 <= x <= 6.6e12; coordinate growth mu*(c - c0) covers beyond it.
- *First contact with multiple zeros, zeros entering from |x| <= X or from large x:* excluded.
  - Wall and finite head, Thm 1.5(i) compactness and backward Hermite splitting.
  - Effectivity of Polymath's C is not needed.
- *Terminal step:* the comparison bounds all zeros of H_T, including the head zeros, which are real. The cleanup then uses only de Bruijn's theorem.
- *Imported theorems within range:*
  - Platt-Trudgian: ordinate 3000175332800 >= (X+1)/2 + 175239886.
  - Wall: t <= 41/250 = 0.164 > T + q_T/2 = 0.10797.
  - Source horizon: 19/125 > T.
  - Polymath Prop 6.1: any real sigma (used on [-3,4]).
  - Thm 1.5(i): 0 < t <= 1/2.
  - The 7/8 statement: the Lean declaration OAI.riemannZeta_ne_zero_of_seven_eighths_lt_re at openai/math fd4aeeb; it is a theorem in our QRH port.

**Tightest points** (all hold, none is a gap)
1. The relation credits are load-bearing (section 4, C3/C3b).
2. Row 0's kernel defect is only 7.1e-8 relative, at d ~ 9.1, q = cap.
3. The startup gap is 9.30e-13 against the 1e-12 cushion.
4. kappa0 = 2e-7 against c_X - c0 = 2.24e-7.
5. The self-cap margin in 23272725 is 1.8e-9.

**Trust that remains**
- The deposit is unrefereed and AI-assisted.
- The transcendental enclosures come from python-flint/Arb.
- The Riemann-Siegel constants rest on Polymath's propositions as published.
- The 2,924-module Lean build of the 7/8 theorem is external; our own port proves the same statement.

## 3. 23272725 (question 3: what REPORT4 did not cover)
- **44 new kernels** 2K_y(d,b) >= lambda (y/h) J_{h,r}(d,b), for 0 < y <= rho h and 0 <= b <= y.
  - Adversarial minimisation in the original variables, independent of the author's polynomial: minimum relative defect 1.0e-5, which equals the far-field limit 1 - lambda/2. No negative value.
  - The reduction to P = (D+v-B)D4 - (lambda/2)QG was checked: J_{1,r} = 2G/D4 and 2K_y = 4y(D+v-B)/Q at 8,800 random points (error 1.2e-14), with D4 > 0. So P >= 0 is equivalent to the inequality.
  - The author's stage INDEPENDENT_GEOMETRY regenerates all 44 Bernstein covers: PASS.
- **Self caps:**
  - r = 23/40: max s_r = 1.50306713000 <= M = 1.503067131; the exact quadratic has minimum 1.8e-9 >= 0.
  - r = 3/5: equality 25/17 at w = 0.
- **7/8 seed** q0 = 562501/10^6:
  - q0 - 9/16 = 1e-6.
  - Startup q0 - (2 + c0 q0)1e-8 - 9/16 = 9.17e-7 > 0.
  - Thm 1.5(i) compactness from t = 1e-8.
  - Terminal Lambda <= t + tau_e + q_e/2 with q_e = 1e-8. Correct.
- **Steps that differ from 22867547:**
  - The matched-log residual |QT - V| <= B - (1-A)D_K - 2Z and the quotient bound Re(V/Q) <= (B-2Z)/(1-A) - (D_K - Re T): re-derived from (log n)c_t(n) = sum d~_t(d)c_t(n/d) >= 0.
  - The paired log deficit uses (c+(1-eps Re zu) + c-(1-eps Re z u^bar))/2 >= c0(1-|Re u|): correct.
  - The frozen per-cell provider selection copies whole certified rows: correct.
  - The Gaussian 11/12 fix is the same as in 23752103.
- **Replay:** RUN.py --audit passes all stages in 683 s.
  - The baseline 22867547 audit, the new geometry and the source rebuild from stored enclosures all pass.
  - All 8 histories pass. The main primary history (14,899 steps) gives exactly 116397726408291/10^15; the alternate gives 0.116397726408312.
  - Not done: the optional native regeneration of the 608 retuned cells. This is the same boundary as REPORT4.

## 4. Replay and independent checks (question 4)
**Author replay, 23752103.** Run with python 3.12.12 and uv python-flint 0.8.0 / FLINT 3.3.1 on macOS arm64 (not the bundled Linux runtime): **all 15 stages PASS**, 2,109 s, one job.

| Stage | s | Re-derives or re-hashes |
|---|---|---|
| native | 716 | re-derives the 469 native-height masses and moments from n = 1..690988, with fresh bins and a heat-log recurrence |
| activation | 5 | re-derives the closed-form maxima by rational bisection on s >= 0 |
| spectral | 11 | re-derives all 67 envelope inequalities (Taylor covers + tails); the saved partition is an untrusted hint |
| heads28 | 68 | re-derives capacities, triads, Gram LDL and near/far credits (1.5M checks) |
| heads39 | 43 | re-derives heads, composite allocations, relation duals (Arb arcsin) and all sign patterns |
| wall | 29 | re-derives all 2,763,952 prefixes at 512 bits plus the transfer fees |
| kernel | 793 | re-derives all 67 boundary polynomials and Bernstein positivity; the partitions are untrusted |
| arcs | 334 | re-derives all 221,010 arc inequalities exactly, with Pade, Newton jets, startup and cleanup |
| affine | 2 | re-derives the near/far affine splits from raw moments |
| overlays | 27 | re-derives the final source fields, refunds, jet circle caps and clock four-corner bounds |
| compact-arcs | 73 | re-derives 128-bit margin witnesses from the regenerated transcript |
| interfaces | 1 | identity joins across layers plus a rebuilt shared analytic debit |
| physical-debit-map | 1 | exact ledger recomputed from primitives |
| arithmetic-credit-map | 1 | composition audit of the credit proposals (feasibility comes from heads28/39) |
| direct-source | 4 | re-derives all 268 fields from primitives |

Integrity: Verify_Deposit.py PASS; verify_bundle.py --endpoint PASS (3,388 files plus the endpoint arithmetic). Both are byte checks only.

**My own checks** (written from the TeX, not from the author's code):

| Check | Scope | Result |
|---|---|---|
| C2 arcs_check.py (exact fractions) | all 221,010 arcs: joins, cells, caps, s <= S_i(t0), 0 < r <= 2(s-o) - 1e-12, ell > 0, q1 < q0, dyadic q1; exact own pair and [6/6] Pade on 6,048 arcs; float ideal endpoint on every arc | PASS. Minimum margins: rate 5e-30 beyond eta, source 1.5e-30, cap 2.106e-5, ell 9.97e-5. Re-integration with my own rates ends at 0.1079649100787704 (exact 0.10796491007877157). T + q_T/2 = B* < 0.107964910079 exactly; startup 9.30e-13 |
| C1 kernel_check.py | 67 rows x 64k (d,q) samples + Nelder-Mead + 40-digit recheck | all positive; minimum relative defect 7.1e-8 (row 0), 8.2e-7 (row 16), about 1e-4 otherwise; own pair nondecreasing; far field <= 1 |
| C3 head_check.py | rows 0,1,8,16,20,26,27,30,44,62: own heat recurrence and filter | raw heads equal the author's to 12+ digits. The L-BFGS torus supremum stays below the phase-free charge (minimum slack 1.8e-4, row 62 t+; a 16-start search moved this supremum by only 5e-6) but above the fully credited charge by 0.0009-0.0108, so the relation credits are needed |
| C3b clock_check.py | rows 0, 8, 26, 60; maximise the true head on {abs(k_j.theta) >= rho_j}, all sign patterns (SLSQP) | the maximum stays below H - cc - d by 0.0048-0.0115: the dual certificates are not falsified, and the constraints cost more than the credit claimed |
| C4 phys_t0_check.py | row 0 at t = 0, straight from xi'/xi (digamma + von Mangoldt sum, 34k prime powers, tail <= 1.6e-4); 5.4M centres in [X0, X0+2e6] | minimum J = 7.7303 > floor 7.6401; J(X0) = 7.8164; main term at X0 8.478318463 >= main_lo 8.47831818 |
| C5 kernels116_check.py | 23272725 kernels, self caps, seed | see section 3 |

Resources: one compute process at a time, nice 19, single thread. Peak memory 645 MB (arcs). About 2 CPU-hours, including an over-long first head run that I killed.

## 5. Consequences
- **N2 record:**
  - Lambda <= 427693072077942473661578098491/3961408125713216879677197516800 < 0.107964910079 (Raymond, Zenodo 23752103).
  - It replaces our 0.1383621 (DBNOpt) and Hakimoglu-Brown's 0.13963.
  - Its inputs: the 7/8 theorem, Platt-Trudgian, Polymath (Thm 1.5(i), Props 3.1/6.1), de Bruijn and Rodgers-Tao.
  - The DBN6 claim "bar = 0.116398" is superseded.
- **Robustness:** with every arithmetic credit removed, the same functionals give about 0.10808. Their most delicate part therefore moves the record by at most 1.1e-4.
- RH relevance: N2

## Files
- evidence/: arcs_check.py, kernel_check.py, head_check.py, clock_check.py, phys_t0_check.py, decompose.py and kernels116_check.py, each with its .out file.
- replay_23752103_RECEIPT.json and replay_23752103_log_tails.txt.
- replay_23272725_audit.log and integrity_23752103.log.
- Downloads and replay work directories were in .scratch/dbnaudit (deleted after the run).
