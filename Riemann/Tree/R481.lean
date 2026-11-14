/-!
# Round 481 — node N8 (critical-line proportions and simple zeros), round 30, audit of a load-bearing input of our Theorem E1: Wang's local pair theorem (arXiv 2609.07918 v2) is unconditional as stated (0 < λ < θ < 1, g ∈ C_c^∞ real even, supp g ⊂ [−λ, λ]), its cited Lamzouri Proposition 2.1 (2609.02882) matches the source, its error bound E_x ≪ 1/x holds for 1 ≤ x ≤ T^λ (a scout's FLAW read x^{1/2}/T² as x^{1/2}/T), c(1/2), c(0.6), θ₀ = 0.5501939647 and c(1) = 0.6725007 reproduce, and the withdrawn 2609.24167 (same author) is not used; no dispute found; the Kubilius model, Ehrhart h*-positivity and Bernstein–Sato polynomials transfer no tool

* Node: N8 (30th round; R479 before). Round type: audit of a load-bearing input (AGENTS.md priority; E1 and the e1arb record
  rest on Wang's unrefereed preprint, so far checked only at statement level); Wang's TeX pre-extracted.
* [audit] statement (TeX l.120–260): thm:local-pair (l.229–240) assumes only fixed 0 < λ < θ < 1 and a fixed real even
  g ∈ C_c^∞ with supp g ⊂ [−λ, λ]; error O_g(H + T^λ L²); no RH, GRH, simplicity or pair-correlation hypothesis anywhere;
  built on the unconditional explicit formula of Baluyot–Goldston–Suriajaya–Turnage-Butterbaugh (Lemma 1): HOLDS.
* [audit] citation: Wang's Proposition 2.1 = Lamzouri 2609.02882 Proposition 2.1 (source l.304–322), same hypotheses, two of
  Lamzouri's conclusions omitted, nothing strengthened: HOLDS.
* [falsify] proof: (eq:E) l.263–266, E_x(t) ≪ 1/x + x^{1/2}/T² + x^{−5/2}/T ≪ 1/x for 1 ≤ x ≤ T^λ: the middle term needs only
  x ≤ T^{4/3} (a scout's FLAW 'λ ≤ 2/3' read T² as T; rejected after reading the line): HOLDS; the lemmas lem_norm,
  lem:localization, lem:coeff and prop:pointwise were not refereed line by line (open).
* [check] c(θ) = 2 − θ/2 − cot(θ/√2)/√2 (prop:optimal l.790–819, the cosine test-function optimum proved): c(1/2) =
  −0.1659638503650, c(0.6) = 0.1345543, θ₀ = 0.5501939647 reproduce; c(1) = 0.6725007037 matches Cheer–Goldston's rounded
  0.6725 (a scout's FLAW on the 7·10⁻⁷ rounding difference rejected): HOLDS.
* [check] status: 2609.07918 v1 7 Sep 2026, v2 6 Oct 2026, no comments or journal-ref; 2609.24167 (Biao Wang, withdrawn
  6 Oct 2026) is not cited by it; no Zenodo or arXiv paper disputes it: HOLDS.
* [explore] Kubilius model → N3: no probability space for Li indices (an additive-function variance is negative at n = 2):
  no transfer.
* [explore] Ehrhart theory (h*-vector nonnegativity) → N3: no map to Li coefficients beyond the RH-equivalent positivity
  (R301 tried Ehrhart): no transfer.
* [explore] Bernstein–Sato polynomials (Kashiwara's rationality) → N7: ζ has no b-function (R299 tried Igusa zeta): no
  transfer.
-/

namespace Fayed.R481
end Fayed.R481
