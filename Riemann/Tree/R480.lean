/-!
# Round 480 — node N4 (Robin/Lagarias), round 13, the unseen Zenodo N4 records (first N4 sweep without a date filter: 93 records for Robin, Nicolas, superabundant and colossally abundant phrases): Chua's two claimed proofs via the Nicolas and Robin criteria break (Lemma 4's unsigned error term cannot give the needed gap; the 2-adic range log log log n ≤ α ≤ log n is never treated), Vega's conditional superabundant route is vacuous (its hypothesis P_n ≥ Q fails for every large superabundant n by its own Proposition 8): three entries in REFUTATIONS; Arneth's signed-Möbius route is conditional as stated; Polak's prime harmonic envelope states no bound; N4 records stand; Garside theory, the polynomial method and condensed mathematics transfer no tool

* Node: N4 (13th round; R427 before). Round type: unseen Zenodo records without a date filter (lesson 227; N4 had none):
  bash .scratch/zsw.sh n4 with 8 phrases ('"Robin inequality"', '"colossally abundant"', '"superabundant"', '"Nicolas inequality"', …):
  93 records, most off-topic or Vega's 2021–2022 series; six on-topic texts pre-extracted (pdftotext).
* [falsify] Chua, 'Nicolas Criterion and Riemann Hypothesis' (20790544, revised 22 Jun 2026): Lemma 4 needs the factor
  (1 + O(e^{−c√log P_k})) to lift 1 + (2+β)/(√log N_k · log log N_k) above the same expression at ⌊k/2⌋, a gap of about
  (√2 − 1)(2+β)/(√log N_k log log N_k); an unsigned, far larger O-term cannot give it; the base 'known for k < k₀' is
  unsupported: FLAW → REFUTATIONS.
* [falsify] Chua, 'On the RH via Robin Criterion' (18323203, 19 Jan 2026, his 'correct proof'): n = 2^α·odd treated only for
  α < log log log n and α > log n, then 'this means' every α; colossally abundant numbers have α ≈ log log n/log 2 in the gap;
  the chain e^γ log log 2n + 2.5/log log 2n < e^γ log log 2^{21}n fails at n = 3 (5.3254 > 4.8993): FLAW → REFUTATIONS.
* [falsify] Vega, 'Robin's criterion on superabundant numbers' (7294629, 2022): Theorem 1 is conditional on P_n ≥ Q =
  1 + 1.5809·10⁻⁶ for all large superabundant n; his Proposition 8 gives log P_n ≤ (2 + o(1))/q_k → 0; colossally abundant
  numbers with largest prime 10⁶, 10⁷, 10⁸ give P_n − 1 = 4.5·10⁻⁷, 5.2·10⁻⁸, 7.4·10⁻⁹ < Q − 1 (orchestrator; the scout's
  'sign error' misread the proof by contradiction): vacuous → REFUTATIONS.
* [check] Arneth, 'Signed Möbius Contraction and the RH' (24201597, 13 Nov 2026; RH-58 24201481): no unconditional claim that
  implies RH; Open Hypothesis B is sufficient, not shown equivalent; the finite inversion (2) holds for x ≤ 100: HOLDS.
* [check] Polak, 'Derived Prime Harmonic Envelope on CA Support' (20548637): Σ_{j≤k} 1/p_j = log log p_k + B +
  C_osc/(√p_k log p_k) with no stated bound on C_osc (its Analytic Target 1 is explicitly unproved); for p_k ≤ 10⁷,
  C_osc ∈ [0.2531352682 (p = 302831), 1.7551596165 (p = 113)]: NOT TESTABLE.
* [explore] Garside theory (stable translation lengths) → N13: zero ordinates carry no braid-monoid order: no transfer.
* [explore] polynomial method (Croot–Lev–Pach / Ellenberg–Gijswijt slice rank) → N7: ordinates mod q are not 3AP-free
  (first six zeros give S₃ = {0, 1, 2}): no transfer.
* [explore] condensed mathematics (solidification) → N1: no trace map or positive polarization for Weil's form: no transfer.
-/

namespace Fayed.R480
end Fayed.R480
