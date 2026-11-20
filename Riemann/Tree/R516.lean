/-!
# Round 516 — node N9 (moments, Lindelöf, subconvexity), round 18, falsification of our barrier explore49 (the determinant method has no grip on 13/84): the fibre bound (≤ 2 counted pairs per generic conic) holds for Q = 2R…8R and α = 0.43–0.49 and Bourgain's chain has no other lossy step, but the exact counts and the conic identity cover odd q only; Markov's inequality, quasiconformal maps and Horn's problem transfer no tool

* Node: N9 (18th round; R503 before). Round type: falsification of our own unaudited barrier (explore49, 14 Nov 2026, NOGO; in no ceilaudit), plus three exploration scouts. explore49's spacing2.c recompiled (cc -O2, .scratch/r516), every run < 6 s and < 20 MB.
* (1) Coverage (scout evenq; orchestrator read spacing2.c l.8 and Bourgain–Watt 1505.04161 TeX l.2225–2275): the counts keep only odd q ∈ [Q, 2Q) (t = 10¹²: 393 arcs, 255 kept), and the conic identity uses A = 4a with [[A, c], [q, d]] ∈ SL₂(ℤ), which needs gcd(4a, q) = 1. For even q the Gauss sum G(a, k; q) vanishes on one parity of k (q ≡ 0 mod 4: k = 2k′ with phase −ā k′²/q; q ≡ 2 mod 4: modulus q/2), so x₂ changes form and mixed-parity pairs are not covered; Watt's Weyl-differenced normalization has coordinates for every q. A gap of scope, not of the conclusion: the barrier's measurements hold for odd q; even q (about a third of the arcs) are untested.
* (2) Q ≠ R (scout qmult; orchestrator reran Q = 8R at t = 10¹⁸): max points per generic conic 1–2 for Q = R, 2R, 4R, 8R at t = 10¹⁶, 10¹⁸ (few pairs at 4R and 8R: 25 at 8R, t = 10¹⁸): HOLDS.
* (3) α ∈ {0.43, 0.47, 0.49} (scout alpha; orchestrator reran α = 0.49 at t = 10¹⁸: 277,328 pairs, max 2 per generic conic, line-pair conics 1.0 %): HOLDS; saving(all4/x2x3)·H² = 0.061–0.069 throughout, the σ = 2 shape.
* (4) Inputs (scout inputs): Bourgain 1408.5794 l.123: 'the only input of this paper is to provide an optimal result for the first spacing problem'; (3.8) is the double large sieve, B₁^{1/2} its Cauchy–Schwarz, R ≤ N ≤ R² checked; no other step loses a power: HOLDS ('attained' is explore49's word, 'optimal' Bourgain's).
* (5) Literature (scout lit320): no statement of 3/20 or θ(σ) in 1408.5794, 1505.04161, 2308.14859; arXiv API queries 'second spacing problem' and 'Bombieri–Iwaniec' with zeta: nothing unseen (Huxley's books not checked).
* (6) Exploration, A. A. and V. A. Markov's inequalities → N8: on P = (t − γ₁)(t − γ₂) the bound is tautological, no gap bound: no transfer (scout markov).
* (7) Exploration, quasiconformal maps (Teichmüller's extremal theorem) → N13: ζ has no canonical Beltrami coefficient or extremal family; |ρ(1 − ρ)| ≤ γ² + 1/4 is algebra: no transfer (scout quasiconf).
* (8) Exploration, Horn's problem (Weyl, Lidskii, Horn inequalities) → N13: they need the self-adjoint operator itself; diagonal encodings of the ordinates are tautological: no transfer (scout horn).
* (9) NOGO correction (orchestrator): explore42's 'an RH rate needs Pavlov A₂ control of |ζ|^{±2} near the line' is superseded by the audited nyman-beurling-continuous-rate (RH alone; interpolation cost bounds instead of Riesz/A₂ bounds).
* Next: an Orca falsification of the even-q classes together with four more unaudited barriers (explore29, explore39, explore42's remaining claims, explore45): ceilaudit8.
-/

namespace Fayed.R516
end Fayed.R516
