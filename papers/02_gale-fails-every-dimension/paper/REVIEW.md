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

### Round 6 (v6): six fresh reviewers

R4, R5: VERDICT: CLEAR. Math: VERDICT: CORRECT. R1: 1 flag (convex hull glossed via the undefined "convex set"). R2: 1 flag (pole distance across is 0.8944 for n = 8, outside "0.88 to 0.89"). R3: 1 flag (the purpose clause for left-over lines skipped the final pole step).

Fixed in v7: convex hull = all weighted averages with nonnegative weights summing to 1; "0.880 (n = 7) or 0.894 (n = 8) across"; "Counting the colours that the Cᵢⱼ in left-over rows and columns need forces g = 1 or h = 1; then a pole shares its colour with an adjacent vertex."

### Round 7 (v7): six fresh reviewers

R1, R3, R4, R5: VERDICT: CLEAR. Math: VERDICT: CORRECT. R2: 4 flags (why the threshold 1 in the Idea; "about" for the rounded pole distances; backtracking and search nodes; DRAT certificate unexplained).

Fixed in v8: "A partition into parts of diameter at most σₙ times the diameter of the set, which is less than 1, puts the poles of each simplex into different parts"; "about 0.880 (n = 7) or 0.894 (n = 8)"; searches "colour one vertex at a time and step back when two adjacent vertices get the same colour (search nodes: partial colourings visited)"; "a DRAT certificate, a machine-checkable list of logical steps checked by the program drat-trim".

### Round 8 (v8): six fresh reviewers

R1, R3, R4: VERDICT: CLEAR. Math: VERDICT: CORRECT. R2 and R5: the same flag ("which is less than 1" seems to say the diameter is below 1).

Fixed in v9: "If D ≈ 1 is the diameter of the set, a partition into parts of diameter at most σₙD < 1 puts the poles of each simplex into different parts". make_figs.py: ruff clean (no unused import, context managers, literals, f-strings), deterministic PDFs (no CreationDate; two runs give identical checksums).

### Round 9 (v9): six fresh reviewers

R1, R3, R4, R5: VERDICT: CLEAR. R2: the Idea paragraph's "every such sharing forces some part to contain two points farther apart than σₙD" suggests each pair's three points act alone. Math: "Rows with n ≤ 1000 are also proved in Lean" overstated for X₃, X₄ (raised because the evidence we gave it listed only their Lean lower bounds).

Fixed in v10: Idea: "For every pair Pᵢ, Qⱼ we add three points … placed so that, together with the poles, they rule out every partition into n+1 parts of diameter at most σₙD (§6)". Lean actually proves both bounds for X₃, X₄ (P3.lean, P4.lean: lower, upper, diam_X; Main.lean beta_eq combines them); Table 1 caption now "For n ≤ 1000, Lean (§8) also proves the bounds in each row", and Table 6 gains the row "lower and upper bounds for X₃, X₄ | lower, upper | P3.lean, P4.lean" (links pinned to a85a25e, L42–L46, L48–L52). Evidence for the next math reviewer corrected.

### Round 10 (v10): six fresh reviewers

R1, R2, R3, R5: VERDICT: CLEAR. Math: VERDICT: CORRECT. R4: 2 flags (the long exact fractions for X₃, X₄ in §7 hide why they are counterexamples; "use only Lean's three standard logical axioms" leaves a non-Lean reader unsure whether geometric facts are assumed).

Fixed in v11: §7 now states "β₄(X₃)/diam X₃ > σ₃ and β₅(X₄)/diam X₄ > σ₄ (Table 1)" (the exact values stay in Table 1, truncated, and in P3.lean, P4.lean); §8 "assume nothing beyond the three standard axioms of Lean's logic"; "This applies to every counterexample to Borsuk's conjecture" (avoids a one-word last line).

### Round 11 (v11): six fresh reviewers

R1, R2, R3, R4, R5: VERDICT: CLEAR. Math: 1 scope flag: §8 "Theorems 1–3 … are formalized" and Table 6 "Theorem 3 → claim2", but claim2 states β_{n+1}(Vₙ)/Dₙ > σₙ, not the intermediate bound β_{n+1} ≥ sₙ (checked: Family.lean rowOK_facts uses a stored threshold tZ, not "the shortest H-edge").

Fixed in v12: §8 "Theorems 1 and 2, the inequality β_{n+1} > σₙDₙ of Theorem 3, the lower bound of Lemma 2 and the parts of Lemma 1 that we use are formalized"; Table 6 row "Theorem 3, β_{n+1} > σₙDₙ | claim2".

### Round 12 (v12): six fresh reviewers

R1–R5: VERDICT: CLEAR. Math: VERDICT: CORRECT.

### v13: author block (Adrian, 27 Nov)

Class article → amsart (reqno; caption margin 0): \author{Adrián Sager}, \address{Independent researcher}, \urladdr{https://www.linkedin.com/in/a-sager/} (printed after the references, clickable); abstract before \maketitle; date as amsart footnote. github.com/Sager611/riemann now only as availability link: Note adds "Code, certificates and Lean files are at github.com/Sager611/riemann; file names in this paper link there"; §8 "The file links point to commit a85a25e". Provenance sentence unchanged. Visual check of all 7 pages: title block, running heads, equation numbers right, Table 6 spacing; main text ends on p-6.

### Round 13 (v13): six fresh reviewers

R2, R3, R4, R5: VERDICT: CLEAR. Math: VERDICT: CORRECT. R1: 2 flags (Note: "certificates" and Lean unexplained at first use; §1 "the lowest dimension of a known counterexample is 64 [5], or 9 by an unrefereed preprint" reads as two competing minima).

Fixed in v14: Note "using computer-checkable certificates (§4) and proofs verified by the program Lean (§8) that assume two cited results" (still one provenance sentence); §1 "the lowest dimension with a refereed counterexample is 64 [5], and an unrefereed preprint gives counterexamples in every n ≥ 9 [8, Cor. 7.1]". Pages 2–7 unchanged.

### Round 14 (v14): six fresh reviewers

R3, R4, R5: VERDICT: CLEAR. Math: VERDICT: CORRECT. R1: "supremum" unexplained in the abstract; §1 "whose splitting step is proved by hand" vague. R2: the move from colourings to a Boolean formula is unexplained.

Fixed in v15: abstract "diameter (least upper bound of the distances between their points)" (matches §2's gloss of infimum); §1 "in which a proof by hand replaces the search over splits and exact arithmetic checks the distances"; Theorem 2's proof "a Boolean formula with one true/false variable for each vertex and colour has a solution exactly when G≥s(X) has a proper (n+1)-colouring in which one 6-clique has fixed colours; a DRAT certificate … proves that the formula has no solution" (checked against verify.py's CNF: variable k·i+c+1 = vertex i gets colour c); \looseness=-1 on that proof (no one-word last line).

### Round 15 (v15, final): six fresh reviewers

R1, R2, R3, R4, R5: VERDICT: CLEAR. Math: VERDICT: CORRECT. No open flags.

### no-ai-slop detection mode on the final text (v15)

Banned words 0; empty adverbs 0; empty phrases 0; binary contrasts 0; puffery 0; trailing -ing clauses 0; weasel attribution 0; recap endings 0; hedges 0. Em dashes: 1, the empty cell of Table 1 (not prose). Colons followed by lowercase: 23, all introducing lists, definitions or glosses (e.g. "a certificate in two parts: a proof that …"), none a colon reveal. Visual check of all 7 pages and of each figure at full size: clean.

### v16: Adrian's wording (27 Nov)

§1 opening replaced verbatim: "In 1933, Borsuk asked whether every bounded set in ℝⁿ can be split into n + 1 parts, each with a smaller diameter than the set. He proved that n parts are not always enough: a ball in ℝⁿ cannot be split into n such parts [1]." Citation kept. p-1 and p-2 re-rendered and checked (one line moves to p-2; pages 3–7 unchanged).

### Round 16 (v16): six fresh reviewers

R1, R3, R5: VERDICT: CLEAR. Math: VERDICT: CORRECT (incl. the Borsuk 1933 sentences). R2: Fig. 2 "Dashed: radius 0.62" names neither style nor panels; the Boolean variables' meaning is unstated. R4: "commit" is Git jargon.

Fixed in v17: Fig. 2 "Grey arcs in (b), (c): radius 0.62 about 0"; Theorem 2's proof "a Boolean formula whose true/false variables say which vertex has which colour has a solution exactly when …" (not "exactly one colour per vertex": verify.py's CNF has only at-least-one clauses); §8 "The file links point to version a85a25e of the repository".

### v17: Adrian's wording (27 Nov), second sentence of §1

Replaced verbatim: "The answer is no for large n: some sets in ℝⁿ cannot be split into n + 1 parts of smaller diameter [6]. This holds for every n > 322 [4], with a refereed counterexample in ℝ⁶⁴ [5]; an unrefereed OpenAI preprint gives counterexamples for every n ≥ 9 [8, Cor. 7.1]." Voronov [9] checked: all 8 citations are in Gale contexts (Fig. 1 disk, Gale's n = 3–6 counterexamples and Voronov's question in §1, σₙ bound and its optimality, Gale's conjecture, refined sets in §7, ball bound in §8); none in a Borsuk sentence, so none removed. \looseness=-1 on the last paragraph of §1 (avoids "(Theo-/rem 1)."). Pages 1, 3, 4, 6 re-rendered and checked.

### Round 17 (v17): six fresh reviewers

R4, R5: VERDICT: CLEAR. Math: VERDICT: CORRECT. R1: §1's last sentence nests qualifications and does not name Hinrichs's theorem. R2: the unit vectors u_i, v_j can be mistaken for the poles (edge length 1). R3: "so they are all of K₁₂" is ambiguous.

### v18: Adrian's own wording (27 Nov, supersedes v17's sentence) and round-17 fixes

§1, after the Borsuk sentences, Adrian's wording verbatim with one factual fix: "Even for large n some sets cannot be split into n + 1 smaller diameter parts [6], in fact it was proven for all n > 322 [4], and then for a smaller n = 64 [5]. Recently, there has appeared a preprint from OpenAI for all n ≥ 9 [8, Cor. 7.1]." Factual fix: "there have appeared preprints" → "there has appeared a preprint": openai/math has one Borsuk preprint (item 156, cited as [8]; Riemann/Deep/Gale7/NOVELTY.md). Not reworded otherwise; reviewer flags on this paragraph are reported, not applied. Voronov [9]: all citations already in Gale contexts.

Round-17 fixes: §1 "With our sets for n = 3, 4, the family for 5 ≤ n ≤ 1000 and Hinrichs's theorem [4] for n > 1000, the conjecture fails for every n ≥ 3 (Theorem 1)"; §5 Construction "…: the vertices of two regular simplices centred at 0; scaled by r_A and r_B (Table 2) they give the poles"; Lemma 2 "so they make up K₁₂". p-1, p-3, p-5 checked; no TeX warnings.
