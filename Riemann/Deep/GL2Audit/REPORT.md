AUDIT PASSED — every reduction from Bernard's theorem to the certified functional now has a complete written proof, both GL(2) certificates (and the ζ one) replay with 0 failures under my rerun, and no admissible choice I tried beats 0.07005633 (0.08971026 under Ramanujan).

# REPORT — gl2audit: referee of discoveries/gl2-levinson-ceiling.md (one-piece Levinson–Conrey ceiling for GL(2))

Agent gl2audit, 12 Nov 2026, 17:40–18:16 UTC (all runs single-threaded under memcap 290 MB; peaks ≤ 62 MB). Full proofs: FULLPROOF.md (§0 input, §1 reductions, §2 Lemma D, §3 cells, §4 tails,
§5 code audit, §6 falsification). No step failed. Two clarifications of the claim's scope, one error found in the input papers
(it does not change their numbers), and two minor code remarks (below).

## Checked, step by step

| Step | Status | Where |
|---|---|---|
| Bernard Thm th:mm2 → c(P,Q,2R,ν/2); zero count → κ_f ≥ 1 − (1/r) log c(P,Q,r,ν/2), r = 2R; ξ = ν/2 | proved from TeX (B l.159–161, 770–795; AT l.191–231) | §0 |
| (E): c = (1+|w(1)|²)/2 + (C_P/ξ)W₀ + ξB_PW₁ for real P, complex Q | proved | §1.1 |
| complex Q of Bernard's shape never helps with real P (c(P,Q) = c(P,Q_r) + ∫∫A_i²) | proved | §1.2 |
| unique P-minimiser sinh(λx)/sinh λ; barrier needs only c ≥ M(R,λ(Q)) (no density) | proved | §1.3 |
| B_λ non-increasing, C_λ non-decreasing; B_PC_P ≥ 1/4; κ* non-decreasing in ξ ⇒ open ν-range | proved | §1.4–1.5 |
| Lemma D (weak duality) | proved, reviewer CORRECT (cell form) | §2, §3 |
| cell bound: λ-subcells, F₀, G₁, K₁″, concave-vs-convex endpoint test | derivatives recomputed; reviewer CORRECT | §3 |
| tails R ≤ 0.5367, λ ≥ 3 (R ≤ 0.6045 + strip), R ≥ 28 with C = 1, λ = 0 | proved, numbers recomputed by hand | §4 |
| replay.py, covcheck.py vs DERIVATION | read line by line, no error | §5 |
| reruns on /tmp copies (memcap 290 MB) | 25/256: REPLAY OK 1257+19 cells, 250 s, peak 59 MB; 1/8: REPLAY OK 1014+3, 214 s, 58 MB; ζ 4/7: REPLAY OK 429+97, 79 s, 56 MB; covcheck COVER OK ×3 | §6.6 |

## Independent numbers (own code, not gl2ceil's)

| Quantity | Page | gl2audit | Code |
|---|---|---|---|
| float ceiling ξ = 25/256 | 0.0700563167 at R = 9.867773, λ = 0.770589 | 0.07005631669 (K = 42, degree 83; limit ≈ 0.0700563167), λ = 0.770589037, max at R = 9.867773 | myceil.py |
| float ceiling ξ = 1/8 | 0.0897102510 at R = 7.688073 | 0.08971025099 (K = 36) | myceil.py |
| ζ, ξ = 4/7 | 0.4089593260 at R = 1.291586 | 0.40895932599 | myscan.py |
| κ(R) off the optimum | — | 0.0630 (R = 15), 0.0377 (28), 0.0218 (35), < 0 for R ≤ 3 | myscan.py |
| certificate vs own minimiser | — | target ≤ D ≤ M at 25 points of 5 cells incl. the worst-margin cell | spot.py |

## Findings (scope and inputs; none breaks the claim)

1. Class = the stated hypotheses, verbatim: P real (B l.161 "P is a real polynomial"; AT l.194 "P ∈ R[x]"), P(0)=0, P(1)=1;
   Q complex of Bernard's shape (B l.770, l.787), any degree; R > 0; ν < 25/128 (< 1/4). P(1) = 1 is no restriction (scaling).
   Not admitted by the input and not covered: complex P, shifts, other mollifiers, other derivative scales.
2. Input error: B l.796–800 and AT l.993–996 display inf_P c with √(AC − (Im B)²). That is the infimum over complex P (Euler–
   Lagrange with the 2i(Im B)P′ term); with real P, as both papers state, the infimum is √(AC)coth(ξ√(C/A)), strictly larger when
   Im B ≠ 0. Direct quadrature (cplxP.py): real P 4722.9017460952 = √(AC) formula; complex P 4720.4201831735 = papers' formula.
   Harmless for both papers (they restrict to real Q, B l.801, AT l.1002). DERIVATION 1.1 notes it for AT; B l.796 is the same.
3. The one natural extension outside the barrier — jointly complex (P, Q), which needs B's theorem for complex P — does not beat
   the ceiling numerically: at the real optimum the imaginary directions decouple and the second variation is positive definite
   (coupling/cost ratio t = 0.2328 at ξ = 25/256, 0.2330 at ξ = 1/8, converged in degree; cplx.py), and BFGS from 12 random complex
   starts per R returns to the real optimum (cplxglobal.py). Not certified. Complex P with real Q provably cannot help (§1.6(i)).
4. Minor, no effect: (a) B l.779 writes [·]², not |·|²; the audit (and AT l.209, LEAN.md) read |·|², the only real-valued reading.
   (b) replay.py/ceil.py do not assert gd(R) ≤ 4Rξ (i.e. K ≤ e^R) in the R ≥ RMAX tail; it holds with margin (4·28·25/256 = 10.94,
   4·22/8 = 11 > π/2). (c) AT21's own values 0.0696 and 0.0896 (AT l.138) lie below the ceilings, as they must.

## What survives, what is not proved

Survives in full: κ*(ξ) ≤ 0.07005633 for ξ ≤ 25/256 and ≤ 0.08971026 for ξ ≤ 1/8, over real P, Bernard's complex Q of any degree
and every R > 0; the ζ side result ≤ 0.40896 at ξ = 4/7 (replayed). Conditional only on the Arb arithmetic of python-flint (both
ceil.py and replay.py use it; the spot check of §6.5 is an independent float confirmation, not a proof). Not proved by anyone:
the ceiling for jointly complex (P, Q) (finding 3), which would anyway need an extension of Bernard's theorem to complex P.
Suggested page edit: name B l.796 next to AT l.996 in the side remark (finding 2).

## Files (Riemann/Deep/GL2Audit/)

FULLPROOF.md; myceil.py (own Ritz minimiser), myscan.py, cplx.py (second variation, complex directions), cplxglobal.py, cplxP.py
(which P-class the papers' formula needs), spot.py (certificate vs own minimiser), flintsem.py (Arb comparison semantics); logs/: flintsem.log, myceil_conv.log, myceil_25_256.log,
myscan.log, cplx.log, cplxglobal.log, cplxP.log, spot_25_256.log, replay_25_256.log, replay_1_8.log, replay_4_7.log, covcheck.log.
Reviewer (verdict-only, two items inline): CORRECT, CORRECT.

RH relevance: N8, N12
