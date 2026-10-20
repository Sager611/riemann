DISCOVERY: Ryazanov (arXiv 1808.10774v11) Theorem 1 fails at the sentence before (4.3): φ*_n is NOT n!-periodic and does NOT equal χ off the intervals [n!k, n!k+1); proved false in Lean for every n ≥ 3 (elementary error; no published critique found).

# Audit: Ryazanov, "A proof of the Riemann hypothesis on zeros of ζ-function" (arXiv 1808.10774v11, 31 Aug 2026)

**Verdict: the proof fails at one sentence.** "Note also that functions φ*_n are n!-periodic because each function a_k, k=1,…,n! is n!-periodic. Thus, φ*_n ≡ 0 on each interval [n!k,n!k+1), k∈ℕ, and coincide with the function χ outside of all these intervals." The premise is false: a_k(t) = {t/k} has period k, so it is n!-periodic only when k | n!. Most k < n! do not divide n!. Everything after that sentence, including the error formula (4.3), has no basis.

## Chain (§2–§4)
1. Theorem D (Báez-Duarte 2003, after t ↦ 1/t): if some sequence in the class 𝔑 = span{a_k} converges to χ = 1_{[1,∞)} in L²((0,∞), t⁻²dt), then RH holds. **Correct** (cited).
2. (4.1)–(4.2): φ*_n = −Σ_{k<n!} μ(k)a_k + δ(n)n!·a_{n!}, with δ(n) = Σ_{k<n!} μ(k)/k. The slopes cancel, so φ*_n is a step function. φ*_n ≡ 0 on (0,1) and ≡ 1 on [1,n!) by Σμ(k)[t/k] = 1. **Correct.**
3. **Periodicity, vanishing on [n!k, n!k+1), φ*_n = χ elsewhere: FALSE.**
4. (4.3): ‖χ−φ*_n‖²_* = Σ_k 1/(n!k(n!k+1)) < ζ(2)/n!² → 0, so RH follows. **Unsupported**, because it rests on step 3.

## Why step 3 is false (for every n ≥ 3)
- On whole numbers, φ*_n(m) = Σ_{k<n!} μ(k)⌊m/k⌋ − δ(n)n!⌊m/n!⌋. So the jump of φ*_n at m is Σ_{d|m, d<n!} μ(d) − δ(n)n!·[n! | m].
- For 2 ≤ j < n! the jump is Σ_{d|j} μ(d) = 0. For n! < m < 2n! the jump is −μ(m), because the only divisor of m left out is m itself.
- Bertrand's postulate gives a prime p with n!/2 < p < n!. Then 2p lies in (n!, 2n!), μ(2p) = 1, and the jump at 2p is −1. But the jump at 2p − n! is 0. So φ*_n is not n!-periodic.
- Neither [2p−1, 2p) nor [2p, 2p+1) is one of the exceptional intervals, yet φ*_n changes by −1 between them. So φ*_n ≠ χ = 1 on at least one of them.
- Numbers (`.scratch/ryazanov/check.py`, exact rational arithmetic): φ*_n ≠ 1 at **every** whole number in (n!, 2n!) for n = 3,4,5,6. Also φ*_n(n!) = 1/5, 1.027, 2.083, 2.203 for n = 3,4,5,6, where the paper claims 0.
- The true error stays large. ∫ over [1, 20·n!] alone gives ‖χ−φ*_n‖²_* ≥ 0.0398, 0.0391, 0.0386, 0.0319 for n = 3..6, against the paper's bound ζ(2)/n!² = 0.0457, 0.0029, 1.1e-4, 3.2e-6. So (4.3) is already false at n = 4, 5, 6.
- The sequence's convergence is not refuted in Lean; only its proof is. The numbers give no sign that it converges.

## Lean (`Riemann/Deep/Ryazanov/Paper.lean`: builds, sorry-free, axioms propext/choice/Quot.sound)
- Defined: `a`, `delta`, `phiStar` (exactly (4.1)–(4.2)), `chi`, `InR` (class 𝔑), `distSq` (L²_*), `TheoremD`, `Converges`.
- Proved: `phiStar_mem` (φ*_n ∈ 𝔑) and `ryazanov_RH : TheoremD → Converges → RH` (the paper's logic). Also `periodic_of_ak`: the inference "every a_k periodic ⇒ φ*_n periodic" is valid; it is the premise that is false.
- Refutations for all n ≥ 3: `akPeriodic_false` (a_{n!−1}(1) = 1/(n!−1) but a_{n!−1}(1+n!) = 2/(n!−1)), `not_periodic`, `coincide_false`. Helpers: `phi_nat`, `jump`, `jump_small`, `jump_two_p`, `bertrand_p`. All of it is bundled in `ryazanov_gap`.
- Not yet imported in `Riemann.lean` (orchestrator). Axiom check: `.scratch/ryazanov/Axioms.lean`.

## Novelty (web_search, Semantic Scholar; no contact with the author)
- No critique of 1808.10774 was found. Semantic Scholar lists 0 citing papers. The search turned up only the author's own earlier "A disproof of the RH" on academia.edu (403, unread).
- The error is elementary (a period-k function is not n!-periodic). It is significant only because it is the paper's sole substantive step. **Pending reviewer judgement.**

## Files
`.scratch/ryazanov/src/RIEMANN_2026_08_31.tex` (v11 source), `check.py` (exact numbers), `cites.json`.
