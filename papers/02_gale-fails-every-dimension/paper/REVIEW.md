# Review log: gale.pdf

Reader reviewers: fresh `reviewer` subagents, role "average STEM reader (BSc engineering/physics/CS), no discrete geometry", full text, captions and page images inline. Math reviewer: fresh `reviewer`, numbers and proofs.

## Re-verification (before round 1)

| Check | Result |
|---|---|
| `verify.py` (points7, points8) | β₈/D = √(8900856812318928/10000000213244525) = 0.94344351301…; β₉/D = √(9023994457279939/10000000157507010) = 0.94994706774…; both > σₙ; G(≥s) 1765 / 2195 edges, search nodes 1327; 23097 and 169807; 801496; 1.0 s and 36.8 s |
| `verify.py` X₃, X₄ (Gale/, Gale45/) | 0.8906210556…, 0.9164673596… |
| DRAT, drat-trim built from source | cert7.drat and cert8.drat: `s VERIFIED` |
| `family.py` | 996 rows n = 5…1000 pass; gap 7.9·10⁻⁴ (n = 5) … 1.43·10⁻³ (n = 7) … 2.9·10⁻⁷ (n = 1000); gap·n² = 0.290 at n = 1000 |
| `family_opt.py 3`, `4` | best 0.875161 < σ₃, 0.910647 < σ₄ |
| Integer σ test (2) | true for n = 3, 4, 7, 8 |
| X₇ vs V₇(θ₇) | every point within 0.1757 of a family point; X₈ norms 0.6166…0.6325 |
| Citations | Crossref DOIs (Gale, Borsuk, Kahn–Kalai, Hinrichs, Jenrich–Brouwer, DRAT-trim, Lean 4, mathlib); Voronov TeX (Gale's conjecture, σₙ, β₂,₃ = √3/2, Problem) |

## Rounds

### Round 1 (v1, 7 pages): six fresh reviewers

Split: R1 §1–3, R2 §4–5 + Fig. 2, R3 §6, R4 §7–8 + references, R5 global order + slop detection, M mathematics. Flags: R1 31, R2 32, R3 31, R4 18, R5 49, M 7. Verdicts: none CLEAR.

Fixed in v2:
- Title without "hence" (R1.1). Abstract: equal-diameter comparison, n ≥ 9 marked unrefereed, how the family is proved (R1.3, R1.5, R1.6, R5.1, R5.3). Note names the two cited hypotheses (R1.7, R1.8).
- §1: "smaller diameter than the set"; "cannot be lowered"; Fig. 1 caption states the lower bound with citations; why Hinrichs is used only above 1000; hand step vs exact arithmetic; X₃, X₄ provenance; shorter; "only n = 7, 8 open" tied to [8] (R1.9–14, R4.1, R4.3, R5.4–6, R5.11, M.6).
- §2: partition before the formula, φ⁻¹(i) written out, infimum wording, scaling, regular simplex, facets, cones as rays, shared boundary points, σₙ is the value of one partition, "counterexample" defined as the stronger σ-condition (R1.15–23, R5.7–10).
- Table 1: proof pointer, "points" column, V_n pointer, why X₇, X₈ are listed beside the family, lower bound vs exact ratio, blank cells explained (R1.24–31, R5.12).
- §4: purpose first, s is a pairwise distance, odd-n condition explained (R2.1–2, R5.13–14).
- §5: Idea paragraph before the construction; orthogonal blocks 3+4 and 4+4; regular simplex; 69/85 counts; statements about positions restricted to θ₇, θ₈; "between the directions"; cell defined, shared poles; coefficient example; counting rule for the 36 classes; files define X₇, X₈; X₇ history and objective; "diameter about 1"; proof: D² as maximum, where the colouring lives, symmetry breaking, CNF encoding incl. fixed 6-clique, H(5,5) ⊆ G≥s(X₈) (2170 edges, re-checked exactly); Table 4 next to Theorem 2 (R2.3–24, R5.15–25, M.3, M.4).
- Fig. 2: panel headings, larger poles, dotted pole polygons, axes and labels in (b), duplicate highlighted points explained, block norms drop signs, X₈ norms 0.6166–0.6325 (R2.25–32, R5.26–28, M.2, M.5).
- §6: θₙ source and search; Table 5 and Fig. 3 before Lemma 2; "and" in Table 5; Fig. 3 bands instead of lines, A₀₀/B₀₀/C₀₀ labelled, "colours mark types"; proof: colours per index, plan sentence, K = used colours, pairwise-to-line argument, meaning of g, h, |I₀| ≤ g + |I₁|, order of exclusions, distinctness reason, |K₁₂| = h, identification of the colour of B, explicit h = 1 symmetry; Theorem 3 wording; hand vs computer; 17 classes explained; gaps "about" 2.898·10⁻⁷ / 1.433·10⁻³; 0.29/n² an observation; n = 3, 4 with ≈ (R3.1–31, R5.29–43, M.1).
- §7: exact ratios and certificates for X₃, X₄; Borsuk ⇒ Gale argued via "every partition has largest part diam S", no infimum step (R4.1–4, R5.44–45).
- §8: Lean and Mathlib explained; σ-form vs ball form; hypotheses; sorry, axioms and kernel explained; Mathlib revision; separate table rows; Lemma 2 lower bound; data paragraph with pinned commit (R4.5–15, R5.46–49).
- References: OpenAI preprint path, venue names expanded, LNCS 8561 (R4.16–18).
- Layout: 10 pt; 6 pages plus references.

Not changed:
- R1.2 author line: the brief prescribes "github.com/Sager611/riemann", no personal names.
- R2.11, R3.25, R5.19, R5.41 list of all 36 classes: a one-sentence counting rule instead; the list is in family.py.
- R2.14 which 15 points were copied: discovery history; X₇ is defined by its file.
- M.7 details the reviewer could not verify: verified here (max displacement 0.1757; 18887 certificate tokens, P7; timings from today's verify.py runs; lean-toolchain; PROOF.md).

### Round 2 (v2, 7 pages): six fresh reviewers

Flags: R1 12, R2 10, R3 4, R4 2, R5 9, M 1 (Table 4 "18887-step case split": the certificate has 18887 tokens). No CLEAR.

Fixed in v3: regular simplex glossed in the abstract; "colouring argument by hand with exact arithmetic checks"; "computer-checkable certificates and machine-checked proofs"; Voronov's ask as "without exhaustive computer search"; family scope "for 5 ≤ n ≤ 1000"; simplex = convex hull, facets = convex hulls of n vertices, purpose of the cone partition first, points in several cones; "this stronger inequality violates"; X₃, X₄, X₇, X₈ listed explicitly; V_n pointer in Table 1; "exhaustive search"; Idea paragraph gives dimensions, centre, pole distances 0.88–0.89, "at least about 0.94"; P–A … Q–C class count; proof of Theorem 2 before Table 4, 10¹⁶s² and 10¹⁶D², DRAT glossed and tied back to colourings, H(5,5) ⊆ G≥s(X₈) ⇒ χ ≥ 10; Fig. 2 caption: small markers, block lengths, purpose of the circle, α per n; §6 search sentence after H is defined; Lemma 2 grid sentence (blocks contain no C-vertices), plain plan sentence; 17-class rule spelled out; Lean: plain words for sorry, axioms, kernel; Lean scope stated exactly (Theorems 1–3, Lemma 2 lower bound, the parts of Lemma 1 used); "replays an 18887-token search certificate".

### Round 3 (v3): six fresh reviewers

Flags: R1 6, R2 6, R3 6, R4 2, R5 9; math reviewer: VERDICT: CORRECT.

Fixed in v4: diameter glossed in the abstract; Note without software names; Voronov's request as "without a computer search over all ways to split the set" and the family's hand step vs exact arithmetic; convex hull glossed; Fig. 1 now shows the inscribed triangle (dashed) and the rays that make the sectors; Table 1 "matching lower and upper bounds, shown truncated"; Idea paragraph names the poles and the added points and states the target "farther apart than σₙ times the diameter"; Cᵢⱼ "inside the right angle at 0 between Pᵢ and Qⱼ"; class enumeration moved into the proof of Theorem 3, 17-class rule as explicit splits; X₇: all 84 points moved, Table 2 approximate; Fig. 2 moved to the top of page 3 (end of §4), caption: none of the three maps increases distances, labels for X₇, block lengths; Theorem 2 proof: DRAT sentence condensed, X₈ hand proof moved after Lemma 2 with a pointer; Lemma 2: contradiction assumption first with the concrete target, facts (i)–(iii) as a list, per-colour choice of K₀ lines, "left over" defined; sₙ/Dₙ ≈ for n = 3, 4; §7 inline ratios, Borsuk remark inside the proof; §8 GaleConjecture vs GaleConjectureBall, hypotheses vs sorry wording, data sentence merged; Table 6 caption kept with its table; top and bottom margins 2.3/2.6 cm: main text 6 pages, references from page 6.

AI-slop detection on v4 (no-ai-slop detect mode, by the author): binary contrasts "an observation, not a proved rate" and "Colours mark types, not a colouring" rewritten ("a rate we have not proved"; "Marker colours show vertex types only"); no banned words, no em dashes in prose, no colon reveals, no summary ending; run-in heads (Idea., Construction., The sets.) kept as standard math-paper labels.

### Round 4 (v4): six fresh reviewers

R4 (§7, §8, references): VERDICT: CLEAR. Math reviewer: VERDICT: CORRECT. Flags: R1 3, R2 4, R3 1, R5 2.

Fixed in v5: diameter glossed as a supremum; infimum glossed, best partition for finite V; cones "cover B^n by congruent pieces", shared points assigned to one piece; Idea: Cᵢⱼ "in the angle between them at 0", how V_n(θ), X₇ and X₈ relate; classes: "the squared distance is a quadratic polynomial in θ with rational coefficients"; DRAT: "a Boolean formula whose solutions are the proper (n+1)-colourings"; Table 4: "checks an exhaustive-search certificate"; n = 3, 4: "ratio of shortest H-edge to diameter"; Fig. 1 label moved off the chord.

### Round 5 (v5): six fresh reviewers

R2, R4, R5: VERDICT: CLEAR. R1: 1 flag (cones are unbounded: intersect them with B^n). R3: 1 flag (say why left-over rows and columns matter when they are defined). Math: 1 flag, the history "counterexamples for n = 3,…,6" was not in the supplied evidence; checked in Voronov's TeX (Table 2: n = 3, 4, 5, 6 with 101, 61, 87, 64 points, each above σₙ; the n = 7 row is commented out), so the text stays and the evidence is added for the next math reviewer.

Fixed in v6: "Intersecting B^n with the n+1 cones … gives n+1 congruent pieces that cover B^n"; purpose clause after the definition of left-over rows and columns.
