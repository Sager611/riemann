AUDIT PASSED: every step of Robin2Adic's Theorem A (Robin for every n > 5040 with v_2(n) ≤ 72) checks out, K = 72 included; only a rounding slip in the side clause of Theorem A′ is wrong.

# R2AAudit — referee of the Robin record v_2(n) ≤ 72 (Riemann/Deep/Robin2Adic, 12 Nov 2026), node N4

Agent r2aaudit, 12 Nov 2026, 21:10–21:50 UTC. Details: FULLPROOF.md. Code and logs: falsify.py, b3check.py, kln_table1.py, logs/ (falsify_id/ex, b3check, kln_table1, cert_rerun).

## Verdict
- **Theorem A holds exactly as written (v_2(n) ≤ 72)**, on the inputs it lists. K = 72 needs KLN's Lemma 4.14 with δ ≈ 0.30. The lemma itself assumes only δ > 0, and its proof
  also needs δ < log H_0(log log H_0 − 1)/2 = 33.9 and V > 1, both met by the rows used (δ = 0.303–0.307). Only Theorem 1.1's summary line says δ ≥ 1, and no step uses it (§7, §9).
  Under the stricter reading δ = 1, K = 71 (E* ≤ 1.7707·10⁻²²).
- **Erratum (Theorem A′, not load-bearing):** the clause "a_j ≥ 2 whenever p_j ≤ 1.07·10¹¹" is rounded the wrong way. √(1/E*) = √(1.14336·10²²) = 1.06928·10¹¹,
  so it should read **p_j ≤ 1.069·10¹¹**. Everything else in A′ checks: a_1 ≥ 73, a_2 ≥ 46, a_3 ≥ 31, a_4 ≥ 26, a_5 ≥ 21, a_6 ≥ 19, p_j^{a_j+1} > 1.14·10²². Theorem B's 2⁷⁶ also checks.
- **Erratum (brief, not R2A):** KLN is arXiv 2101.12263. The id 1708.06672 in the brief belongs to an unrelated physics paper.

## What was checked (line by line, FULLPROOF.md §1–§8)
| step | result | how |
|---|---|---|
| Prop 2: exact identity, definitions of B, D, H, E, m ≤ k, q_j ≥ p_j, signs, T > 1 | ✓ | re-derived; GL's PDF definitions match DERIVATION; Arb tests (below) |
| Lemma 3: every slot, a_j = 1, unused slots, small v_2, even/odd | ✓ | re-derived; for small x, c' = min(1/((T+U)log(T+U)), F(x)/L) makes it hold for every n ≥ 6 |
| threshold n > 10^(10^46.79) ⇒ x ≥ e^{108.57} | ✓ | log Λ = 108.5719899; RS62, Bertrand, BKLNW two-sided (TeX l. 5752: M_0 = m_0 = 2.4530·10⁻¹²) |
| (5), G(d,L), (11), (12), (14), (16) | ✓ | re-derived; termwise integration via Ingham's ψ_1; BKLNW l. 5734: M_0 = 1.1199·10⁻⁹ |
| Abel summation (running min, ½ and B_UP, c_j ≥ 0), τ, closed form of W(s) | ✓ | re-derived: telescoping identity, tail = 2W(s), 2W(½) ≤ τ |
| monotonicity in x, term by term | ✓ | derivatives |
| KLN at δ < 1 | ✓ given [Ram, Lemma 5.5] as KLN quote it (τ > 1) | all δ occurrences in the TeX; fallback via Ramaré–Zuniga Alterman (margin ≈ 25); the KLN code reproduces Table 1 to 4.4·10⁻⁴ |
| K = ⌊−log₂E*⌋ − 1 | ✓ | −log₂ 8.7461707·10⁻²³ = 73.276 ⇒ 72; 2⁻⁷³ > E* |

## Falsification (no counterexample to any lemma)
- falsify.py (Arb, 320 bits): 1.08·10⁵ cases. Every n ∈ [6, 40000]; all 58778 n ≤ 10²⁰ with non-increasing exponents (incl. all 120 superabundant numbers);
  684 colossally abundant numbers (largest prime ≤ 10723); 2^v·(odd part of CA), v ≤ 80; 3996 random even n with v_2 ≤ 80; 3000 primorial multiples.
  Results: identity residual ≤ 3·10⁻⁹²; B, D, H ≥ 0; the slot inequality never fails (equality exactly at B = H = 0); GL's c-chain never fails; Robin ⇔ B + D + H > E always;
  Robin fails exactly at the 23 known n ∈ [6, 5040], cross-checked against σ(n) as an integer; whenever E(x) < 2^{−v_2−1}, Robin holds (3253 cases).
- E(x) at all 664411 primes in [10³, 10⁷]: E > 0 throughout; E ≤ 0.36 × the majorant with small-x constants; E√x/f(L) ∈ [1.872, 2.021];
  the normalised zero part lies in [−0.0115, 0.0310], inside ±0.0468 = the bound S_low(1+ε). The x-dependence is the predicted 2f(L)x^{−1/2}.
- b3check.py: X·Σ_{X<n<5X} λ_X(n)²/n² = 0.4534, stable for X = 10⁴ … 4·10⁶, below Ramaré's b_3 = 0.605 (inside KLN, never checked here before).

## Not checked here (inputs, outside this brief)
- Ramaré's Lemmas 5.5 and 5.6 in their own words (Math. Comp. 2016; paywalled, not on arXiv). KLN quote them for τ > 1. Partial cover: the RZA bound and b3check above.
- KLN's proof at H_0 = 3·10¹² beyond the δ scan. Lemma 4.14 is printed with H_0 = 3.0610046·10¹⁰; we rely on Ingham4's reading, which follows the precedent of Cully-Hugill–Johnston's refereed Lemma 7.1.
- RobinCH's range proof, Bebel's half-plane (audited as hpaudit, TREE N7), CLMS, HSW, Platt–Trudgian: used as stated.
- cert.py: rerun here on an unchanged copy, identical (8.7461707·10⁻²³, K = 72; 39 MB). check.py (independent float code) was not rerun.

## Reviewer
Verdict-only reviewer on the δ step (pasted inline): "GAP in 'every δ > 0', not in allowing δ < 1. No quoted step requires δ ≥ 1, and δ = 0.303 satisfies the displayed parameter
restrictions." The gap is in my wording only: the lemma also needs δ < 33.9 and V > 1, and the rows used meet both. Corrected in §7 and §9.

## Process note
The first b3check.py run reached 1036 MB RSS for about 2 s, over the 300 MB cap: it allocated whole-range numpy arrays, and memcap's 3 s sampling missed it.
Rewritten in chunks and rerun: 119 MB RSS, 19 MB footprint. The first full falsify.py run was stopped by memcap at 299 MB (its test sets were held as lists);
rewritten with generators, it peaked at 136 MB. Every other run stayed ≤ 136 MB RSS.

RH relevance: N4
