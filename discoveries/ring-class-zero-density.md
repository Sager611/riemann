# Zero-density estimate for ring class L-functions

[Back to DISCOVERIES](../DISCOVERIES.md)

**Claim.** For a negative discriminant D and 0.6 ≤ σ ≤ 1, Σ_χ N(σ; t, χ) ≪_ε |D|^{2σ(1−σ)+ε} (|t| + 3)³, where χ runs over the ring class characters of D and N(σ; t, χ) counts the zeros β + iγ of L(s, χ) with β ≥ σ and |γ − t| ≤ 1.

**Before.** Every published density bound for this family is weaker in |D| for each σ in (1/2, 1) (Pascadi–Thorner 2025, Kowalski–Michel 2002).

**Method.** Classical: Carlson's zero detector and a weighted large sieve over primitive ideals. The method is known; the statement is new.

**Status.** Paper proof; audited by rcaudit (12 Nov 2026): AUDIT PASSED — Lemmas 1–5 and Theorem 1 (0.6 ≤ σ ≤ 1) re-derived line by line, Theorem 2 and Corollary 3 correct under (L); numerical tests passed (Lemma 1 and the exact sharp constant of Lemma 2 for all 5000 D ≥ −10⁴: at most 1.5·10⁻³ of the claimed constant; Lemma 4 for 208 D by PARI); the σ-range was missing here (for σ < 1/2 the range-free bound is false by ρ ↦ 1 − ρ̄; 1/2 < σ < 0.6 is a routine extension); still new against Pascadi–Thorner 2508.14888 Thm 1.3 and the Oct–Nov 2026 listings; RCApps2 Prop 1 / Thm 1′(a) give the smaller exponent: Σ_χ N(σ; t, χ*) ≪ (|D|T)^{2(1−σ)/(3−2σ)+ε} for 0.6 ≤ σ ≤ 1, from the same refereed Lemmas 2–5 (refereed by rc5audit, 12 Nov 2026: correct, [REPORT2](../Riemann/Deep/RC5Audit/REPORT2.md)). Not in Lean. [RCAudit](../Riemann/Deep/RCAudit/REPORT.md)

**Review.** Core argument correct (weighted large sieve; per unit window).

**Novelty.**
- An arXiv sweep of about 105 queries found nothing implying it ([novelty](../Riemann/Deep/Combos/NOVELTY.md)).
- An adversarial re-check (Zenodo, Hexagon series, OpenAI, ATLAS, GitHub, classical papers; exponents compared at σ = 3/4 and 0.9) still found it new ([novelty 2](../Riemann/Deep/Combos/NOVELTY2.md)).

**Files.** [note](../Riemann/Deep/Combos/NOTE.md)
