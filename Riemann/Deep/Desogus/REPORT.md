DISCOVERY: Desogus (arXiv 2609.20367v3) Theorem 8.6, ineq. (201), is false: on window log 8 an explicit odd test v has Q(v)/⟨e_7,v⟩² ≤ 1.1e-10 (Arb-rigorous), against the claimed ≥ g_7 ≥ 1.878.

# Audit: Desogus, "The Three Gates" (arXiv 2609.20367v3, 6 Oct 2026; Zenodo v4 the same)

**Verdict: the proof fails at Theorem 8.6, the induction step's coercivity inequality (Gate II).** The inequality is false by 10 orders of magnitude; it is not merely unproved. All the other links in the chain are correct or plausible.

## Chain (paper; window L = 2a = log Y, endpoints a_k = ½ log k)
1. Thm 1.2, restricted odd Weil criterion: RH ⟺ W[f] ≥ 0 for every real odd C_c^∞ f. **Correct** (checked the proof).
2. Lemmas 1.3 and 1.4: zero extension, so positivity on cofinal endpoints gives every a. **Correct.**
3. Certificate 6.64: A_{7,−} ≻ 0 (window log 7), smallest LDL pivot 9.2e-25. **Plausible**: our Galerkin λ_min(odd, log 7) ≈ 1.6e-24, and the window is below Mori's log 15. A referee found its production chain incomplete (gap F5).
4. Gate II: Mellin cells, Feshbach/Young branch costs, a "one-defect" lower model with uniform gap log 2 (Lemma 6.48, Cor. 6.49), MASTER-P2/P3 debit ledger, target diagonal B_k = ½ + log(1/w_k) > ½ + log k (Lemma 4.1, eq. 45).
5. Cert. 8.4: g_k^{M+} = ½ + log k − (439/(250 log 2))Σ_m V_{k,m} − (5/2)W_k > 0 for k ≥ 7 (interval sweep, min 1.878 at k = 33).
6. **Thm 8.6 (201)**: if A_k ≻ 0, then for every f in the collar K_k, ⟨H_k f, A_{k+1} H_k f⟩ ≥ g_k^{M+}|⟨e_k,f⟩|² + R⊥ with R⊥ ≥ 0. Hence T_k ≻ 0, so A_k ≻ 0 ⇒ A_{k+1} ≻ 0.
7. Induction from k = 7, then 2, then 1, gives RH.

## Failing step: (201)
- H_k f = (−A_k^{-1}B f, f) minimises the form over every old-core completion. So (201) implies Q(v) ≥ g_k⟨e_k,v⟩² for every odd v on window log(k+1). Here e_k is the unit collar "true ground": q_Y(t) = t^{-1/2} (Lemmas 6.4 and 6.8), which in physical coordinates is sgn(y)e^{-|y|/2} on a_k < |y| < a_{k+1}.
- **Rigorous refutation** (`.scratch/desogus/rigorous.py`; Arb balls, 3000 bits, `acb.integral`, exact rational G(s) = ∫P(x)P(x−s)): v is an explicit odd polynomial (Galerkin solution of A v = e). The explicit-formula pieces (pole, archimedean, primes) use the convention verified against zeta zeros to 1e-28 (Riemann/Deep/Weil).

| k (window log(k+1)) | deg v | Q(v)/⟨e_k,v⟩² (enclosure) | paper's bound |
|---|---|---|---|
| 7 (log 8) | 79 | 1.1016027e-10 ± 3e-18 | ≥ 1.878 |
| 7 (log 8) | 39 | 5.6298056e-7 | ≥ 1.878 |
| 8 / 12 / 33 / 100 | 39 | 1.6e-6 / 2.5e-5 / 3.8e-3 / 0.26 | ≥ 1.878 |

  The ratio keeps falling as deg v grows. The true infimum 1/⟨e,A_{k+1}^{-1}e⟩ is driven by A_8's near-null odd modes (λ_min ≈ 1e-29).
- The result does not depend on which ground vector is meant. The sinh-on-collar polar vector t_k and the plain sign-indicator give the same order (`refute.py`).
- **No debit accounting can save it**: the bare diagonal is already too small. Q(e_7) = 1.2263 (mpmath, `diag.py`; Galerkin projections 1.142 → 1.164 → ↑) is below g_7 ≥ 1.878, and ⟨He,AHe⟩ ≤ Q(e). The paper's target diagonal B_7 = 2.513 does not match the operator (Q(e_7) = 1.23; Q(e_12) = 1.73 vs B = 3.03; Q(e_33) = 2.73 vs B = 4.01).
- Smooth-test version: dilate and mollify v. C_c^∞ is dense in the log-weighted form domain, and Q is continuous there, so the strict inequality survives.
- Root cause: every Gate II debit bound uses local lower models with O(1) gaps (log 2, β = 0.048, Carleman d*). But the real Schur debits invert compressions of A_k, whose spectrum reaches 1e-24 to 1e-29 and decays doubly exponentially in L (Zhu). Any lower bound on the Schur complement that is uniform in k with an O(log k) coefficient is incompatible with that spectrum.
- The step's *conclusion* (A_k ≻ 0 ⇒ A_{k+1} ≻ 0) is not refuted: it holds if RH holds. Its proof is refuted.

## Lean (`Riemann/Deep/Desogus/Paper.lean`, builds, no sorry, no axioms)
`OddReal`, `WeilPosOdd`/`WeilPosOddStrict` (over `Weil.Defs.weilQ`). Proved: `WeilPosOdd.mono` (L1.3), `cofinal` (L1.4), `all_endpoints` (induction), `desogus_RH : Base7 → (∀k≥7, Step k) → OddCriterion → RH`.
Failing step: `Thm86Quant k g` (consequence of (201)); `Witness7` (the Arb fact, cited Prop); `thm86_false : Witness7 → 1e-6 ≤ g → ¬Thm86Quant 7 g`; `desogus_gap` bundles both halves.

## Novelty (web_search, Zenodo; no contact with anyone)
- G. Chen, "Referee Report on The Three Gates" (Zenodo 10.5281/zenodo.22892645, 22 Sep 2026, on v2). It lists 5 *gaps*: F1 no global harmonic-graph decomposition; F2 target diagonal B_k not derived from the operator; F3 fold model identified circularly; F4 Thm 8.5 (now 8.6) debit accounting; F5 Y=7 chain. It says the checks "do not reveal a simple numerical mistake". It gives no counterexample.
- v3/v4 claim to close all of these ("seven-point ledger"). Not found elsewhere: a proof that (201) is *false*, the quantitative witness, or the fact that B_k is inconsistent with Q(e_k). **New, pending reviewer judgement.**

## Files
- `.scratch/desogus/refute.py` (Galerkin, float-free Arb midpoints)
- `rigorous.py` and `rigorous.log` (the certificate; about 40-75 s each, 1 core, < 60 MB)
- `diag.py` (Q(e_k)), `diag_check.py` (Galerkin cross-check, λ_min of A_7 and A_8)

Paper text: /tmp/riemann/desogus.txt. Referee report: /tmp/riemann/referee.txt.

## Reviewer note (orchestrator)
- Dichotomy: (201) needs A_7 ≻ 0 to define H_7. If A_7 is not ≻ 0, the base certificate 6.64 fails; otherwise (201) fails by the witness. Either way the proof breaks.
