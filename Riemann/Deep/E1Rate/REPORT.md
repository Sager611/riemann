RATE (moment only): Conrey's Möbius structure (Vaughan + Deshouillers–Iwaniec 1984 Lemma 1 + Weil, with the window's Bessel phase separated at cost √(T/H)) proves E1-M for ν < ν_K(θ) = min{1/2, max{(3θ−1)/4, (4/7)(2θ−1)}}, beating (3θ−1)/4 exactly for θ > 9/11 (and SL's own Bettin–Chandee range for θ ∈ (9/11, 67/68)), but h(θ) is unchanged on the whole table because Wang's c(θ) dominates wherever ν improves, and for θ ≤ 9/11 (θ ≤ 3/5 for any conceivable θ = 1 input) the Bessel phase is a barrier.

# E1Rate — Kloosterman input for the short-window mollifier (N8; e1rate, 14 Nov 2026)

## Answer
1. **Mollifier length: yes, on (9/11, 1).** ν_C(θ) = (4/7)(2θ − 1) > (3θ − 1)/4 ⟺ θ > 9/11 = 0.8182; at θ = 1 it is Conrey's 4/7 (consistency ✓).
   With SL's cap ν < 1/2 (its diagonal lemma, l. 285–316, and K-branch): new range (4/7)(2θ−1) on (9/11, 15/16], 1/2 on [15/16, 67/68);
   for θ ≥ 67/68 SL's ν_B = min{1/2, (17/33)(2θ−1)} already reaches 1/2.
2. **h(θ): no.** Where ν improves, c(θ) ≥ 0.5103 > h(θ;κ) and > (c+2κ)/3 for every κ ≤ κ(1/2) = 0.3635; the lower bound max{0, c, (c+2κ)/3, h}
   stays c(θ) (margin ≥ 0.0648 on [9/11, 1] even for κ = 0.409 > Conrey's 4/7 one-piece ceiling 0.40896: dominance.py/.log, so lifting
   the cap cannot help either). The detector matters only for θ ≲ 0.6, where no Kloosterman gain exists.
3. **Barrier below 9/11.** Separating the Bessel weight turns the short-window dual into the θ = 1 off-diagonal at height T' = H²/T (transfer
   principle, DERIVATION §8): every separated method with θ = 1 range ν1 gives (2θ−1)ν1 — BC 17/33 → SL's ν_B, DI 4/7 → ν_C, ν1 = 1 → 2θ−1 —
   which beats E1's unseparated (3θ−1)/4 only for θ > (4ν1−1)/(8ν1−3): 9/11 (DI), 35/37 (BC), 3/5 (ν1 = 1). As θ ↓ 1/2 all → 0; E1 keeps 1/8.

## Table (κ: optimal Levinson–Conrey detector in SL's class, kappa_opt.py, re-evaluated in mpmath; h = h(θ;κ); LB = max{0, c, (c+2κ)/3, h})
| θ | old ν (max E1, SL ν_B) | new ν_K | κ old → new | h(θ;κ) old → new | c(θ) (Wang) | LB old = new |
|---|---|---|---|---|---|---|
| 0.51 | 0.1325 (E1) | 0.1325 | 0.094636 | 0.023474 | −0.1300 | 0.023474 |
| 0.60 | 0.2000 (E1) | 0.2000 | 0.144010 | 0.140866 | 0.1346 | 0.140866 |
| 0.75 | 0.3125 (E1) | 0.3125 | 0.228095 | 0.295702 | 0.4191 | 0.419075 |
| 0.82 | 0.3650 (E1) | 0.365714 | 0.267653 → 0.268187 | 0.355579 → 0.355908 | 0.5103 | 0.510319 |
| 0.90 | 0.4250 (E1) | 0.457143 | 0.311746 → 0.334478 | 0.416074 → 0.429785 | 0.5931 | 0.593102 |
| 0.95 | 0.463636 (SL ν_B) | 0.5 (cap) | 0.338975 → 0.363501 | 0.450397 → 0.465022 | 0.6357 | 0.635679 |
| 0.99 | 0.5 (SL ν_B, cap) | 0.5 (cap) | 0.363501 | 0.477890 | 0.6656 | 0.665554 |
ν values are suprema (open ranges). Uncapped (conditional on Conrey's ν ∈ [1/2, 4/7) main term in the window, not done) θ = 0.99: ν = 0.56,
κ = 0.401271, LB unchanged 0.665554. Side note (routine, not Kloosterman): E1's own range with an optimised detector already gives
h(0.51) → 0.023474 (E1 §4 printed 0.016352 with the fixed ν = 3/25 detector); float optimum, mpmath re-evaluation, not Arb-certified.

## The bilinear form and the input (DERIVATION §1–§3)
SL (14) after E1 Step 1 + reciprocity: O = Σ_g (4π/g) Σ_{(p,q)=1} Σ_n β_p γ_q α_n e(n q̄/p) Φ(n/(pq)), p, q ≤ y/g, n ≲ N0 = pqT/H², with
β, γ the mollifier coefficients μ·(smooth)/p, /q and Φ = e(−x)𝒲(x), |Φ| ≍ H X0^{−1/2}, Bessel phase X0 = 2√(2π)√(nT/(pq)) ≤ 2√(2π)T/H.
Modulus p (or q), inverted variable q (or p), dual variable n. Mellin separation costs Λ = √X0 (SL (15); F1 numerics), so a bound with gain G
gives relative error N·G. Vaughan on q (cut U = y^{1/4}): type II via DI Lemma 1 (Baluyot L3.1.6: A, B ∈ [U, y/U], V = P, length N);
type I via Lemma I (Poisson + Weil + second-derivative test: conductor p·X0). Binding: DI's spectral term at the balanced split A = B = y^{1/2},
N0^{1/2}y^{−1/8} = T^{(1−2θ)/2+7ν/8}, and type I at a = U, T^{1−2θ+7ν/4}. exponents.py (all blocks, monomials, moduli, orders, U-grid):
ν_K = max{(3θ−1)/4, (4/7)(2θ−1)} at θ = 0.51, 0.6, 0.75, 0.8, 0.8182, 0.85, 0.9, 0.9375, 0.99, 1.

## Falsification (DERIVATION §6)
- Smallest gain case θ = 0.82, ν = 0.3657: E1's bound fails (+0.0014), DI spectral term −0.000013, type I −0.000025, all other DI monomials
  ≤ −0.045, K-branch/residues need ν < 1/2 ✓; the margins are exactly 7/8·(ν_C − ν) and 7/4·(ν_C − ν).
- θ just above 1/2 (0.51, ν = 0.1325): DI exponent +0.106, type I +0.212 — far worse than E1; ν_C(0.51) = 0.0114. No gain, no error.
- F1: Mellin L¹ cost / sup|Φ| grows like √X0 (ratio/√(1+X0) → 0.09 at X0 = 640–2560): the √(T/H) charge is real.
- F2: Lemma I holds in 24/24 cases (max|S|/bound ≤ 0.19); single sums are ≈ √D (r-sums of Kloosterman sums cancel) — usable only with
  an average over the modulus p (open).
- Optimizer sanity: κ(0.12) = 0.085244 ≥ E1's frozen detector 0.084478; κ(0.56) = 0.40127 just below Conrey's 4/7 one-piece 0.4088–0.40896.

## Review
Reviewer (Claims A–D inline with exact values): GAP — c2 range a ∈ (U, U²] said to be uncovered; derivative hypotheses for Lemma I / Mellin
cost; C(a,ℓ) must not depend on (b, v). Resolved (DERIVATION §7: (a, k) is a bilinear pair under DI with a inside, max exponent
T^{(1−2θ)/2+7ν/8}; hypotheses are E1's audited Steps 2–4; coefficients separate). Re-review of the c2 resolution: CLOSED ("the strict claimed
condition controls every term"; ν < 4δ/7 implies N0 < U). Exponent arithmetic of the balanced block and the 9/11 crossover: confirmed.
Not re-derived here: DI 1984 Lemma 1 itself (taken as stated in Baluyot's audited thesis; its θ = 1 output 4/7 is reproduced exactly).

## Prior art / pre-emption
- BC route PRE-EMPTED by SL v0.03 itself (§3.4, l. 494–545: ν_B = min{1/2, (17/33)(2θ−1)}; reproduced here exactly).
- Conrey-structure route in windows: Conrey 1989's Gaussian-window proposition ν < 4θ − 24/7 (SL l. 190–192) and Das–Pujahari arXiv 2104.10243
  Thm 1 (window version of Conrey's Thm 2; their bound gives ν < 4a − 24/7, their statement 2a − 10/7 takes the larger of two conditions —
  misstatement, l. 487–494); both < (4/7)(2θ−1) for every θ < 1. Tang 2608.14852: single-pair reciprocity, no mollifier.
- Scans (arXiv listings math.NT 2026-10/11, arXiv API, Zenodo, scout's 8 arXiv queries): (4/7)(2θ−1) in windows T^θ not found.

## For the orchestrator (gate)
- Not a discovery for N8's zero counts: h(θ) and every lower bound for simple/critical zeros in windows are unchanged.
- RECORDS.md candidate: short-window mollified moment (E1-M) for ν < min{1/2, (4/7)(2θ−1)} on θ > 9/11 (before (3θ−1)/4, SL ν_B).
- NOGO.md candidate (barrier): separating the Bessel weight of the short-window dual reduces it to the θ = 1 off-diagonal at height H²/T,
  so Kloosterman/large-sieve inputs give ν1(2θ−1): nothing beyond E1 for θ ≤ 9/11 (DI), θ ≤ 3/5 (any ν1 ≤ 1).
- REFUTATIONS/minor: Das–Pujahari 2104.10243 Thm 1 range misstated (2(7a−5)/7 should be 4(7a−6)/7 for a < 1, from their own display).

## Files
DERIVATION.md (§1 form, §2 separation, §3 DI/Weil, §4 rest of SL + cap, §5 obstruction, §6 falsification, §7 review, §8 transfer barrier,
§9 prior art); exponents.py/.log (grid minimax over blocks; peak 14 MB); kappa_opt.py/.log, kappa_opt.json (detector optimum; peak 47 MB);
h_table.py/.log, h_cols.py/.log (mpmath κ and h tables; peak 15 MB); falsify.py/.log (F1, F2; < 50 MB).

RH relevance: N8

