NO DISCOVERY: Zhang's step (19) confuses a conjugate-pair factor with a zero quadruplet; an off-line quadruplet is two distinct factors (α and 1−α), so "l ≠ i" is no contradiction. Elementary and circular; formal refutation new but not significant.

# Audit: Weicun Zhang, "A Proof of the Riemann Hypothesis via a New Expression of ξ(s)", Preprints.org 202108.0146 v60 (24 July 2026)
PDF: `.scratch/zhang/p.pdf`, text `.scratch/zhang/p.txt` (curl is blocked by the CDN; Node `fetch` works).

**Verdict: Lemma 6 (⇒ of eq. 13) is false; it fails at eq. (19).** §1, Lemmas 3, 4, 5, 7 and the convergence work (22)–(24) are fine.

## Chain
1. (22)–(25): pair ρᵢ with ρ̄ᵢ, so ξ(s) = c ∏ gᵢ(s)^{mᵢ}, where gᵢ(s) = 1 + (s−αᵢ)²/βᵢ² and ∑1/βᵢ² < ∞. **Correct.**
2. (16)–(17): from ξ(s) = ξ(1−s), divisibility and Lemma 3, every gᵢ(s) equals some gₗ(1−s). **Correct**: it only restates that the zero set is symmetric under s ↦ 1−s.
3. **(19) fails.** Quote: "due to the uniqueness of the multiplicity mᵢ, the only solution to Equation (17) is" l = i, "since otherwise" (19c) gives "αᵢ + αₗ = 1, βᵢ² = βₗ², l ≠ i, which means the i-th and l-th quadruplets of zeros coincide …, violating the uniqueness of mᵢ."
4. (20)/(30): αᵢ = 1/2. This depends entirely on step 3.

## Why step 3 is false
- gᵢ has only the two zeros αᵢ ± iβᵢ. The paper also calls mᵢ "the multiplicity of the i-th quadruplet", but that label only makes sense when αᵢ = 1/2.
- If a zero lies off the line, its quadruplet {α ± iβ, 1−α ± iβ} contributes two different factors, gᵢ (with α) and gₗ (with 1−α), both with the same multiplicity. The map s ↦ 1−s swaps those two factors. Case (19c) describes exactly this situation; it is not a coincidence of two quadruplets.
- So "l = i for every i" just says every zero is on the line, which is RH. The argument assumes what it sets out to prove.

## Davenport–Heilbronn test: fails
The argument uses only the Hadamard product, the functional equation, 0 < α < 1 and ∑1/β² < ∞. The paper's own remark says DH escapes because its zeros are not all in the strip. That does not help:
(a) the counterexample below lies inside the strip;
(b) DH has a certified zero 0.8085 + 85.699i, which is also inside the strip;
(c) DH's zeros have bounded real part, so ∑1/β² < ∞ holds for DH anyway.

## Lean: `Riemann/Deep/Zhang/Paper.lean`
Compiled with `lake env lean`; no `sorry`; axioms are only propext, Classical.choice and Quot.sound. Model: g, F = ∏' gᵢ^{mᵢ}; `Hyp` is Lemma 6's full hypothesis list (strip, β ≠ 0, |β| monotone, ∑1/β² < ∞, mᵢ ≥ 1, distinct factors, convergence).
- `chain`: (17) plus `Step19` imply `Lemma6`. The deduction itself is valid (`alpha_half` does the coefficient comparison in (19a)).
- `hyp0`, `symm0`: α = (3/10, 7/10, 1/2, 1/2, …), β = (1, 1, 2, 3, …), m ≡ 1 satisfies every hypothesis, and F(s) = F(1−s) because s ↦ 1−s permutes the factors (`Equiv.swap 0 1`, `Equiv.tprod_eq`).
- `not_lemma6 : ¬ Lemma6` (α₀ = 3/10). `not_step19 : ¬ Step19` (i = 0, l = 1). `div17_0`: (17) holds, so the failure is at (19).
- `dh_zero_in_strip`: DH has a zero with 0 < Re < 1 and Re ≠ 1/2 (from `Fayed.DH.Zero.offLineZero_holds`).

## Prior critiques
A search of MathOverflow, Math.SE, PubPeer, Reddit, Zhihu and Google Scholar found no public rebuttal. The acknowledgements thank several reviewers, and PNAS rejected the paper. The flaw is elementary: this "pairing" fallacy is a known pattern among amateur RH proofs. Significance is low. The author was not contacted.
