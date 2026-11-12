FORMALIZED: Gale's conjecture fails in ℝ⁷ and ℝ⁸ — Lean proves ¬GaleSplit 7 X7 and ¬GaleSplit 8 X8 (no partition into n + 1 parts of diameter ≤ σₙ·diam) for the explicit family members V₇ (69 points, ℝ³ ⊕ ℝ⁴) and V₈ (85 points, ℝ⁴ ⊕ ℝ⁴) with the exact parameters of family_params.txt, plus β₈(V₇) > 0.9433983·D > σ₇·D and β₉(V₈) > 0.949947·D > σ₈·D; the colouring lemma of NOTE §3 is proved for all p, q and the reduction (NOTE §3 Theorem) for the whole family. No named hypotheses; sorry-free; standard axioms only; no native_decide. Not covered: the 84-point points7.txt (DRAT certificate, no hand proof), the 10⁻⁸-rounded files, and σₙ ≥ β_{n+1}(Bⁿ) (the regular-simplex partition of the ball, the cited source of σₙ).

# Gale7 in Lean 4 + Mathlib (galelean, 2026-11-12)

Lean v4.35.0-rc2, Mathlib v4.35.0-rc2 (repo lakefile). 13 files, 1,314 lines, [Lean/](Lean/), namespace `Gale7`. Source: [NOTE.md](NOTE.md) §1–§3, [discovery](../../../discoveries/gale-r7-r8.md).

| file | content (NOTE) | true peak |
|---|---|---|
| Lemma | `Vx p q` (vertices P, Q, A, B, C of H(p, q)), `Adj` (edges), colouring lemma `card_le`, `not_colourable`, `lemma_H` (§3, all p, q) | 276 MB |
| Desc | `ix1`, `ix2`, orbits `Orb`, rational class value `clsQ`, edge descriptors `adjD` (17 classes), `adj_adjD`, `neD`, `ne_neD` | 264 MB |
| Fam | `sqd`, `IsUnitSimplex`, `famPt` (Vₙ, §2), class formula `sqd_famPt`, rational form `sqd_famPt_cls` | 251 MB |
| Defs | `eucl`, `diam`, `sigma`, `GaleSplit`, `galeSplit_iff` (§1) | 286 MB |
| Crit | `long_pair`, `diam_image_le`, `not_galeSplit_of_H` (§3 Theorem) | 283 MB |
| Simplex3, Simplex4 | unit regular simplices in ℝ³, ℝ⁴: u_i(k) = c_ik √ρ_k; Gram by `decide +kernel` on integers | 286, 288 MB |
| Q7, Q8 | parameters, t, D²; all 100 descriptors: class ≤ D², edge classes ≥ t, distinct-vertex classes > 0 (`decide +kernel`) | 258, 261 MB |
| N7, N8 | `pt7`, `X7`, `hD7`, `hH7`, `hσ7`, `not_galeSplit_X7`, `X7_long`, `diam_X7`, `sigma7_lt`, `ratio7_gt` (same for 8) | 288, 286 MB |
| Main | `card_X7` = 69, `card_X8` = 85, `gale_false_7`, `gale_false_8`, `beta7_gt`, `beta8_gt` | 288 MB |
| Check | `#print axioms` | 280 MB |

## Main statements (Lean; constructor and edge lists compressed)

    noncomputable def sqd {n : ℕ} (x y : Fin n → ℝ) : ℝ := ∑ k, (x k - y k) ^ 2
    noncomputable def eucl {n : ℕ} (x y : Fin n → ℝ) : ℝ := Real.sqrt (sqd x y)
    noncomputable def diam {n : ℕ} (X : Finset (Fin n → ℝ)) : ℝ :=
      if h : X.Nonempty then (X ×ˢ X).sup' (h.product h) (fun z => eucl z.1 z.2) else 0
    noncomputable def sigma (n : ℕ) : ℝ :=
      if n % 2 = 0 then Real.sqrt (((n : ℝ) + 1) / ((n : ℝ) + 2))
      else Real.sqrt (1 / 2 + Real.sqrt (((n : ℝ) - 1) / ((n : ℝ) + 3)) / 2)
    def GaleSplit (n : ℕ) (X : Finset (Fin n → ℝ)) : Prop :=
      ∃ c : (Fin n → ℝ) → Fin (n + 1), ∀ k, diam (X.filter fun x => c x = k) ≤ sigma n * diam X

    inductive Vx (p q : ℕ) | P : Fin p → Vx p q | Q : Fin q → Vx p q
      | A : Fin p → Fin q → Vx p q | B : Fin p → Fin q → Vx p q | C : Fin p → Fin q → Vx p q
    def Vx.Adj : Vx p q → Vx p q → Prop    -- all 25 ordered pairs of constructors:
      -- P i,P k: i ≠ k | P k,A i _: k ≠ i | P k,B i _: k = i | P k,C i _: k ≠ i | Q j,Q l: j ≠ l | Q l,A _ j: l = j
      -- Q l,B _ j: l ≠ j | Q l,C _ j: l ≠ j | A i _,A k _: i ≠ k | A i _,C k _: i ≠ k | B _ j,B _ l: j ≠ l
      -- B _ j,C _ l: j ≠ l | C i j,C k l: i ≠ k ∧ j ≠ l | P-Q, A-B: False | the reversed orders alike
    theorem Vx.card_le (hcol : ∀ x y, Adj x y → col x ≠ col y) [Fintype K] : p + q ≤ Fintype.card K
    theorem Vx.lemma_H (hpq : 1 ≤ p + q) (col : Vx p q → Fin (p + q - 1)) : ∃ x y, Adj x y ∧ col x = col y

    def IsUnitSimplex {m : ℕ} (U : Fin (m + 1) → Fin m → ℝ) : Prop :=
      ∀ i k, ∑ t, U i t * U k t = if i = k then 1 else -1 / (m : ℝ)
    noncomputable def famPt (U : Fin (m₁+1) → Fin m₁ → ℝ) (W : Fin (m₂+1) → Fin m₂ → ℝ) (rA rB a b c d e f : ℝ) :
        Vx (m₁ + 1) (m₂ + 1) → Fin (m₁ + m₂) → ℝ
      | P i => Fin.append (rA • U i) 0          | Q j => Fin.append 0 (rB • W j)
      | A i j => Fin.append (a • U i) ((-b) • W j) | B i j => Fin.append ((-c) • U i) (d • W j)
      | C i j => Fin.append (e • U i) (f • W j)
    theorem not_galeSplit_of_H {n p q : ℕ} (hpq : n + 2 ≤ p + q) (pt : Vx p q → Fin n → ℝ) {t D2 : ℝ}
        (hD0 : 0 ≤ D2) (hD : ∀ v w, sqd (pt v) (pt w) ≤ D2) (hH : ∀ v w, Vx.Adj v w → t ≤ sqd (pt v) (pt w))
        (hσ : sigma n ^ 2 * D2 < t) : ¬ GaleSplit n (Finset.univ.image pt)

    noncomputable def U3 : Fin 4 → Fin 3 → ℝ := fun i t => (c3 i t : ℝ) * Real.sqrt ((w3 t : ℝ) / 9)
      -- c3 = ![![1,1,1], ![-1,1,1], ![0,-2,1], ![0,0,-3]], w3 = ![6,2,1]
    noncomputable def U4 : Fin 5 → Fin 4 → ℝ := fun i t => (c4 i t : ℝ) * Real.sqrt ((w4 t : ℝ) / 48)
      -- c4 = ![![1,1,1,1], ![-1,1,1,1], ![0,-2,1,1], ![0,0,-3,1], ![0,0,0,-4]], w4 = ![30,10,5,3]
    theorem U3_simplex : IsUnitSimplex U3          theorem U4_simplex : IsUnitSimplex U4
    def rA7 : ℚ := 612372435696 / 1000000000000   -- … f7: family_params.txt row 7 digits / 10¹²
    def t7 : ℚ := 356000182722431308912917 / 400000000000000000000000
    def D7 : ℚ := 100000000000102500044289 / 100000000000000000000000
    noncomputable def pt7 : Vx 4 5 → Fin 7 → ℝ := famPt U3 U4 rA7 rB7 a7 b7 c7 d7 e7 f7
    noncomputable def X7 : Finset (Fin 7 → ℝ) := Finset.univ.image pt7
    noncomputable def pt8 : Vx 5 5 → Fin 8 → ℝ := famPt U4 U4 rA8 rB8 a8 b8 c8 d8 e8 f8   -- row 8
    -- t8 = 90239946609080831680751/10²³, D8 = D7

    theorem not_galeSplit_X7 : ¬ GaleSplit 7 X7
    theorem not_galeSplit_X8 : ¬ GaleSplit 8 X8
    theorem gale_false_7 : ∃ X : Finset (Fin 7 → ℝ), X.card = 69 ∧ ¬ GaleSplit 7 X
    theorem gale_false_8 : ∃ X : Finset (Fin 8 → ℝ), X.card = 85 ∧ ¬ GaleSplit 8 X
    theorem beta7_gt (c : (Fin 7 → ℝ) → Fin 8) : ∃ x ∈ X7, ∃ y ∈ X7, c x = c y ∧ 0.9433983 * diam X7 < eucl x y
    theorem beta8_gt (c : (Fin 8 → ℝ) → Fin 9) : ∃ x ∈ X8, ∃ y ∈ X8, c x = c y ∧ 0.949947 * diam X8 < eucl x y
    theorem sigma7_lt : sigma 7 < 0.9419652        theorem sigma8_lt : sigma 8 < 0.9486834
    theorem X7_long (c : (Fin 7 → ℝ) → Fin 8) : ∃ x ∈ X7, ∃ y ∈ X7, c x = c y ∧ (t7 : ℝ) ≤ sqd x y
    theorem diam_X7 : diam X7 = Real.sqrt (D7 : ℝ)  theorem card_X7 : X7.card = 69   -- same for 8
    theorem ratio7_gt : (0.9433983 : ℝ) < Real.sqrt ((t7 : ℝ) / (D7 : ℝ))   -- √(t7/D7) = 0.9433983553…

## Hypotheses and sources

No named hypotheses; the colouring lemma, the class formula, the Gram matrices, all class values and all σ comparisons are proved. σₙ: NOTE §1 and verify.py `above_ball`; Voronov, arXiv 2610.28828, main_en.tex l. 96–105 (σₙ defined by exactly this closed form, with β_{n+1}(Bⁿ) ≤ σₙ from the regular-simplex partition), Gale's conjecture l. 90–92 ("for every n ≥ 2 the ball is extremal, β_{n,n+1} = β_{n+1}(Bⁿ)"), citing D. Gale, Proc. AMS 4 (1953) 222–225. The step β_{n+1}(Bⁿ) ≤ σₙ is not formalized: GaleSplit is the brief's σₙ-form of Gale's property, so the theorems refute Gale's conjecture through that classical bound (as NOTE §1 does).

## Build (one module at a time, each under the 290 MB watchdog)

    bash Riemann/Deep/Gale7/Lean/audit/buildall.sh   # clean rebuild; per module (audit/mk.sh):
    # LEAN_PATH=.scratch/galelean/olean:<lake env LEAN_PATH> bash .scratch/memcap.sh 290 1800 /usr/bin/time -l lean -j1 -o <olean> M.lean
    # order: Lemma Desc Fam Defs Crit Simplex3 Simplex4 Q7 Q8 N7 N8 Main Check

Clean rebuild 13:19–13:26 UTC (audit/buildall.sh), all rc = 0, 17–71 s per module; true peaks (column above, max 288 MB) from `/usr/bin/time -l` "peak memory footprint" (audit/footprints.txt; the watchdog's 3-s samples read 163–276 MB). Only deprecation and linter warnings (`if_pos`, `push_neg`, unused `done`). rg over Lean/*.lean: no `sorry`, `admit`, `axiom`, `native_decide`, `skipKernelTC`, `implemented_by`. Check.lean (audit/Check.log), identical for all 24:

    'Gale7.gale_false_7' depends on axioms: [propext, Classical.choice, Quot.sound]
    ... (gale_false_8, beta7_gt, beta8_gt, not_galeSplit_X7, not_galeSplit_X8, X7_long, X8_long, diam_X7, diam_X8,
    sigma7_lt, ratio7_gt, sigma8_lt, ratio8_gt, Vx.lemma_H, Vx.card_le, Vx.not_colourable, not_galeSplit_of_H,
    sqd_famPt, galeSplit_iff, card_X7, card_X8, U3_simplex, U4_simplex): [propext, Classical.choice, Quot.sound]

Watchdog kills along the way: one Simplex file with `fin_cases` + `norm_num` over both Gram matrices (297 MB), Simplex4 by `fin_cases` rows (291 MB) → integer Gram by `decide +kernel`; an EuclideanSpace bridge (import `Mathlib.Analysis.InnerProductSpace.PiL2`, 310 and 311 MB) → dropped. The first layout (Family/Classes/N7 with the `decide` checks inside) passed the watchdog but had true peaks up to 304 MB → split into Desc/Fam/Q7/Q8 (no `Real.sqrt` import: its import alone is 270 MB).

## Proof route

1. Lemma (NOTE §3), as in the NOTE: (a) `inX_unique`/`inY_unique`, (b) `inX_of_C`/`inY_of_C`, (c) `C_line`; colour types `Kx` (X(i)), `Ky` (Y(j)), `Kxy`, `Kl`; `R0`, `C0`; p ≤ |Kx| + r₀, q ≤ |Ky| + c₀, r₀, c₀ ≤ z; each L class gets one line (`rowsL`/`colsL`), |R₁| + |C₁| ≤ ℓ; the cells of (R₀∖R₁) × (C₀∖C₁) need distinct XY classes. With a = |R₀∖R₁|, b = |C₀∖C₁| this gives a, b ≥ 1 and a·b ≤ z ≤ a + b − 1, so a = 1 or b = 1 (the NOTE's tightness step, in this form); if a = 1 the b classes of C_{i*j} are all XY classes, so P_{i*} and B_{i*j₀} share one (edge P_iB_ij); b = 1 symmetric with Q_{j*}, A_{i₀j*}.
2. Metric: `sqd_famPt` expands |α u_i − α' u_k|² + |β v_j − β' v_l|² with the Gram matrix; `sqd_famPt_cls` casts it to the rational `clsQ` of the descriptor (orbits, index agreements). `Q7`/`Q8` evaluate `clsQ` on all 5·5·2·2 descriptors in the kernel; `adj_adjD` maps every edge of H to a descriptor with `adjD = true`.
3. `not_galeSplit_of_H`: a Gale partition is an (n+1)-colouring of H (p + q = n + 2) with same-colour pairs at distance ≤ σₙ√D² < √t, but `Vx.not_colourable` gives a same-colour edge of length ≥ √t.

## Fidelity

- ℝⁿ: coordinate vectors `Fin n → ℝ` with `eucl x y = √(Σₖ (x k − y k)²)`, the Euclidean distance (= `dist` in `EuclideanSpace ℝ (Fin n)` by `EuclideanSpace.dist_eq`, `Real.dist_eq`, `sq_abs`; that one-line bridge was not compiled, import over the cap). Mathlib's `dist` on `Fin n → ℝ` (sup metric) is not used.
- diam: the largest pairwise distance of the finite set (`le_diam`, `exists_eq_diam`; 0 only for ∅); `diam_X7 = √D7` is attained by Q₀, Q₁.
- σₙ: NOTE §1's closed forms verbatim (`sigma_even_sq`, `sigma_odd_sq`); σ₇ < 0.9419652, σ₈ < 0.9486834 (NOTE: 0.9419651451, 0.9486832981).
- Gale(n, X) (brief: "∃ partition into n + 1 parts each of diameter ≤ σₙ·diam X"): `GaleSplit` = a colouring `c : ℝⁿ → Fin (n+1)`, part k = `X.filter (c · = k)`. Empty parts are allowed, so ¬GaleSplit also excludes partitions into fewer, or exactly n + 1 non-empty, parts. `galeSplit_iff`: same-colour pairs at distance ≤ σₙ·diam X.
- H(p, q): `Vx.Adj` lists NOTE §3's edges in both orders (`adj_symm`, `adj_irrefl`); the class A_ijB_ij is not an edge (NOTE: "not used"). `lemma_H` = NOTE's Lemma (p + q − 1 colours); `card_le` is the same with any finite colour type. Uniform in p, q (NOTE needs p, q ≥ 1; Lean needs p + q ≥ 1).
- Vₙ (NOTE §2): `famPt` with ℝⁿ = ℝ^{m₁} ⊕ ℝ^{m₂} as `Fin.append` (first m₁ coordinates, then m₂), P_i = r_A u_i, Q_j = r_B v_j, A_ij = a u_i − b v_j, B_ij = −c u_i + d v_j, C_ij = e u_i + f v_j, i ≤ m₁, j ≤ m₂. `IsUnitSimplex` = "unit vectors of regular simplices (u_i·u_k = −1/m₁)". U3, U4: the standard unit-edge simplex coordinates (1/√(2k(k+1)), −k/√(2k(k+1)), 0) scaled to circumradius 1, i.e. coordinates in ℚ(√…); times the parameters they reproduce points7_sym69.txt and points8.txt (the discovery's 69- and 85-point files, rounded to 10⁻⁸) to 4.8·10⁻⁹ per coordinate and the sorted squared distances to 2.4·10⁻⁸ (.scratch float check). X7, X8 are the exact-parameter sets, not the rounded files.
- Parameters: family_params.txt rows n = 7 (split 3 + 4) and n = 8 (4 + 4), 12-decimal digits as exact rationals (reviewer re-checked all 16). t = the smallest of the 17 edge classes (AC, i ≠ k, j ≠ l), D² = the largest class (QQ; for n = 8 also PP): exactly family.py's values.
- Descriptors: `clsQ` = family.py `class_d2`; `adjD` = family.py `EDGES` (17 classes, 21 ordered constructor pairs); the checks run over all 100 descriptors, including unrealisable ones, whose values repeat realisable classes or vanish (same point); D² and t agree with family.py's 36 / 17 classes (Python, exact).
- Long-distance graph (brief: pairs at squared distance > σₙ²D²): `hH7` puts every edge of H(4, 5) at squared distance ≥ t₇ > σ₇²D₇, so `X7_long` says G(≥ t₇) ⊆ G(> σ₇²D₇) of V₇ has no proper 8-colouring (V₈: H(5, 5), 9 colours).
- Quantitative: β₈(V₇) > 0.9433983·D and β₉(V₈) > 0.949947·D are stated without a β definition (every 8-/9-partition has a part with two points farther apart); exact √(t/D²) = 0.9433983553…, 0.9499470859… (`ratio7_gt`: > 0.9433983). The discovery's headline 0.9434435130 is the 84-point points7.txt (not formalized).

## Reviewer (statement fidelity, verdict-only, statements and NOTE §1–§3 pasted inline)

One pass ([audit/review_prompt.txt](Lean/audit/review_prompt.txt), verdict [audit/review1.txt](Lean/audit/review1.txt)): **FAITHFUL** — "Blocking: none. sqd and eucl are squared Euclidean distance and Euclidean distance … diam is the actual finite maximum; its empty-set value creates no loophole. sigma matches both closed forms exactly … Allowing empty classes weakens GaleSplit, so its negation is at least as strong as excluding exactly n+1 nonempty parts … All sixteen parameters match … Adj matches §3 in both endpoint orders. adjD represents precisely the 17 canonical edge classes … not_galeSplit_X7 and not_galeSplit_X8 are unconditional statements about the explicit real-coordinate sets." Minor caveat: no standalone β_k or threshold-graph definition, no exact optimality of β/D. After the review only the file split (statements unchanged) and U3's form (integer c3 and weights w3/9 instead of real cM3 and ρ3 = (2/3, 2/9, 1/9): same numbers).

## Not formalized

- β_{n+1}(Bⁿ) ≤ σₙ (simplex partition of the ball): cited (Voronov l. 96–105, NOTE §1).
- points7.txt (84 points, 0.9434435130) and the rounded 10⁻⁸ files; β = dc exactly (the explicit colouring of G(>dc²)); n ≠ 7, 8 (the lemma and `not_galeSplit_of_H` are general; per n only the 100 rational checks and unit simplices in ℝ^{m₁}, ℝ^{m₂} are needed).

## .lake

Left as it is: the real directory (another agent's Mathlib, 7.4 GB) was only read (no `lake build`, no `lake exe cache get`, nothing written there); my oleans went to .scratch/galelean/olean (deleted after the final build; audit/buildall.sh rebuilds them in 2–7 min). The orchestrator restores the placeholder.

## FORMALIZATIONS.md entry

- Gale's conjecture fails in ℝ⁷ and ℝ⁸ (our discovery [gale-r7-r8](discoveries/gale-r7-r8.md)) machine-checked in Lean 4 + Mathlib: Gale7.gale_false_7/8 — the explicit 69- and 85-point family members V₇ ⊂ ℝ³ ⊕ ℝ⁴, V₈ ⊂ ℝ⁴ ⊕ ℝ⁴ (exact family_params.txt parameters) admit no partition into n + 1 parts of diameter ≤ σₙ·diam; β₈(V₇) > 0.9433983·D, β₉(V₈) > 0.949947·D; the NOTE §3 colouring lemma for all p, q; no named hypotheses, sorry-free, standard axioms, no native_decide (class checks by `decide +kernel`), 1,314 lines; σₙ ≥ β_{n+1}(Bⁿ) cited; points7.txt (84 points) not covered. [LEAN.md](Riemann/Deep/Gale7/LEAN.md)

Orchestrator (12 Nov 2026, 13:4x UTC): accepted from the build logs (copied from .scratch/galelean, now deleted, to [Lean/audit/logs/](Lean/audit/logs/): all 13 modules rc = 0, true peaks 251–288 MiB, watchdog peaks ≤ 276 MB) and a static scan (no sorry, admit, native_decide, skipKernelTC or axiom declarations; Check.log: standard axioms only).
