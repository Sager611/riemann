# LeanGale plan (T2: Gale's conjecture fails for every n ≥ 3)

Status: BLOCKED-MEMCAP before the first file (STATUS.md). This plan is for the resume.

## Layout (all `module` files, private imports; Real.Sqrt as module = 246 MiB leaves ≈ 40 MiB per module)
1. Lemma: copy Gale7/Lean/Lemma.lean (H(p,q) not (p+q−1)-colourable; module files cannot import Gale7's plain files).
2. Fam: `famPt`, class formula (copy Gale7 Fam/Desc), unit simplices for every m ≥ 1 (Helmert: u_i(t) = c_it·√((m+1)/(m(t+1)(t+2))), c_it = 1 (i ≤ t), −(t+1) (i = t+1), 0; Gram by telescoping Σ 1/((t+1)(t+2))).
3. Rows: family_params.txt n = 5…1000 as integers (×10¹²); one Bool checker per row (36 class values ×m₁m₂·10²⁴, 17 edge classes ≥ t, all ≤ D², σ test: even (n+2)t > (n+1)D², odd (2t−D²)² (n+3) > (n−1)D⁴ with 2t > D²); `decide +kernel` in chunks of ≈ 100 rows per module.
4. P8: points8.txt (85 pts) = V₈ on the 10⁻⁸ grid; labels give H(5,5) ⊆ G(≥ dc²), min over H edges = dc² = 9023994457279939/10¹⁶ exactly (tools/hom.py) → lemma; plus an explicit 9-colouring of G(> dc²) → β₉ = dc exactly.
5. P7: points7.txt (84 pts): no H(4,5)-homomorphism into G(≥ dc²) by nearest labels (tools/hom.py: 68 nodes, none) → kernel-checked search certificate (verify.py forward checking: 1327 nodes, clique colours fixed w.l.o.g.; adjacency as integer bitmasks) + explicit 8-colouring of G(> dc²).
6. Lit: n = 3, 4: Voronov arXiv 2610.28828 (unrefereed; check TeX statement: set, ratio, σ or ball value); n > 1000: a source stating Borsuk counterexamples in every n ≥ N, N ≤ 1001 (candidates: Hinrichs–Richter 2003 n ≥ 298; Jenrich–Brouwer 2014; OpenAI Cor. 7.1 d ≥ 9, unrefereed; Ji 2608.12561) — exact wording to check; Borsuk ⇒ Gale proved in Lean (σₙ < 1).
7. Statement (needs PiL2): σ n explicit; `GaleConjecture n := ∀ V : Set (EuclideanSpace ℝ (Fin n)), Bornology.IsBounded V → ∃ c : EuclideanSpace ℝ (Fin n) → Fin (n+1), ∀ k, Metric.diam (V ∩ c ⁻¹' {k}) ≤ σ n * Metric.diam V`; β_k(V) := ⨅ over colourings of the max part diameter; bridge `EuclideanSpace.dist_eq`; main theorem ∀ n ≥ 3, ¬ GaleConjecture n; claim (1) exact ratios β₈(X₇)/diam = √(8900856812318928/10000000213244525), β₉(X₈)/diam = √(9023994457279939/10000000157507010); claim (2) ∀ n ∈ [5, 1000], ∀ unit simplices U, W, ¬ Gale split of Vₙ(U, W).
