FORMALIZED: Theorem A, both directions (char 𝕜 = 0, G finite connected simple, n ≥ 2: A(G) ≅ A_RW(G) ⟺ G regular, or bipartite biregular with d₁/d₂ ∈ (𝕜^×)³), via Theorem N "⇒" (BCS22 Thm 2.4(i) "⇒" proved, not assumed), Theorem I (1)+(2), Lemma 2, Corollary 4.1 and §5; Corollary A: A(K_{1,2}), A(K_{1,3}), A(K_{2,3}) ≇ A_RW over ℚ (the survey's statement is false over ℚ), and over 𝕜 the statement holds for all G ⟺ every positive rational is a cube in 𝕜 (cube roots of all primes suffice). No named hypotheses; sorry-free; standard axioms only. Not compiled: the instances 𝕜 = ℝ, ℂ (their imports exceed the 290 MB cap).

# EvoAlg in Lean 4 + Mathlib (evolean, 2026-11-12)

Lean v4.35.0-rc2, Mathlib v4.35.0-rc2 (repo lakefile). 15 files, 1,725 lines, Riemann/Deep/EvoAlg/Lean/, namespace `EvoAlg`. Source: [NOTE.md](NOTE.md) §0–§5, [discovery](../../../discoveries/evolution-algebras.md).

| file | content (NOTE) |
|---|---|
| Defs | `evoMul`, `EvoIso`, `adjM`, `rwM`, `AIsoRW`, `IsRegular`, `IsBiregular`, `BiregularCube`; natural-basis and square lemmas (§0) |
| Basic | `ind`; coordinates of `φ`; `QP = PQ = I`; trace trick `card_eq_of_support`; `lemma1` (§2.2 Lemma 1) |
| Lemma2, Lemma2M | Lemma 2: identity, square form, Mathlib-matrix form, block form (§4) |
| Diag | diagonal isomorphisms; "if" of Theorem A (§5 "Conversely") |
| TwinI | twin classes; Lemma 1 for A(G)→A_RW(G); τ well defined, injective, `τ ∈ Aut H`, `k_{τS} = k_S`, `d_{τS} = d_S`, (a), (b) (§3) |
| Cor41 | `scaling_ii` (Lemma 2 applied), `theoremI`, `cor41` (§3, §4) |
| TheoremA | `sec5`, `theoremA` (§5) |
| TheoremN | Theorem N "⇒", "⇐" for any non-degenerate evolution algebra (§2.2) |
| TheoremI2 | Theorem I (2) (§3) |
| CorA, Examples, CorA2, CorA3 | Corollary A (§0): K_{1,2}, K_{1,3}, K_{2,3} over ℚ; K_{1,8}, K_{2,16} iso over ℚ; `corA_iff`, `corA_primes` |
| Check | `#print axioms` |

## Main statements

Context: `variable {V : Type*} [Fintype V] [DecidableEq V] {𝕜 : Type*} [Field 𝕜] (G : SimpleGraph V) [DecidableRel G.Adj]`.

    def evoMul (M : V → V → 𝕜) (x y : V → 𝕜) : V → 𝕜 := fun w => ∑ v, x v * y v * M v w
    def EvoIso (M M' : V → V → 𝕜) : Prop :=
      ∃ φ : (V → 𝕜) ≃ₗ[𝕜] (V → 𝕜), ∀ x y, φ (evoMul M x y) = evoMul M' (φ x) (φ y)
    def adjM : V → V → 𝕜 := fun v w => if G.Adj v w then 1 else 0
    def rwM : V → V → 𝕜 := fun v w => if G.Adj v w then ((G.degree v : 𝕜))⁻¹ else 0
    def AIsoRW : Prop := EvoIso (adjM G 𝕜) (rwM G 𝕜)
    def IsRegular : Prop := ∃ d : ℕ, G.IsRegularOfDegree d
    structure IsBiregular (s t : Set V) (d₁ d₂ : ℕ) : Prop where
      bip : G.IsBipartiteWith s t
      cover : ∀ v, v ∈ s ∨ v ∈ t
      deg₁ : ∀ v ∈ s, G.degree v = d₁
      deg₂ : ∀ v ∈ t, G.degree v = d₂
    def BiregularCube : Prop :=
      ∃ (s t : Set V) (d₁ d₂ : ℕ), IsBiregular G s t d₁ d₂ ∧ ∃ c : 𝕜, c ≠ 0 ∧ (d₁ : 𝕜) / d₂ = c ^ 3
    theorem evoMul_single_single (v u : V) :
        evoMul M (Pi.single v 1) (Pi.single u 1) = if v = u then M v else 0
    theorem adj_sq (v : V) : evoMul (adjM G 𝕜) (Pi.single v 1) (Pi.single v 1)
        = ∑ w ∈ G.neighborFinset v, (Pi.single w 1 : V → 𝕜)
    theorem rw_sq (v : V) : evoMul (rwM G 𝕜) (Pi.single v 1) (Pi.single v 1)
        = ((G.degree v : 𝕜))⁻¹ • ∑ w ∈ G.neighborFinset v, (Pi.single w 1 : V → 𝕜)

    -- Theorem A (NOTE §0, §5)
    theorem theoremA [CharZero 𝕜] (hG : G.Connected) (hn : 2 ≤ Fintype.card V) :
        AIsoRW G 𝕜 ↔ IsRegular G ∨ BiregularCube G 𝕜
    theorem theoremA_if [CharZero 𝕜] (hG : G.Connected) (hn : 2 ≤ Fintype.card V)
        (h : IsRegular G ∨ BiregularCube G 𝕜) : AIsoRW G 𝕜
    theorem aIsoRW_of_biregular {s t : Set V} {d₁ d₂ : ℕ} (hb : IsBiregular G s t d₁ d₂)
        {c : 𝕜} (hc : c ≠ 0) (hcube : (d₁ : 𝕜) / d₂ = c ^ 3) : AIsoRW G 𝕜   -- any field
    theorem sec5 [CharZero 𝕜] (hG : G.Connected) (hnb : ∀ v, ∃ w, G.Adj v w) (c : V → 𝕜)
        (hc0 : ∀ v, c v ≠ 0) (hi : ∀ u v v', G.Adj u v → G.Adj u v' → c v = c v')
        (hii : ∀ u v, G.Adj u v → c v ^ 2 = c u * G.degree u) : IsRegular G ∨ BiregularCube G 𝕜

    -- Corollary 4.1 (NOTE §4), any field, twin-class sizes non-zero in 𝕜
    theorem cor41 (hnb : ∀ v, ∃ w, G.Adj v w) (hdeg : ∀ v, (G.degree v : 𝕜) ≠ 0)
        (hk : ∀ u, (#(twin G u) : 𝕜) ≠ 0) :
        AIsoRW G 𝕜 ↔ ∃ c : V → 𝕜, (∀ v, c v ≠ 0) ∧
          (∀ u v v', G.Adj u v → G.Adj u v' → c v = c v') ∧
          (∀ u v, G.Adj u v → c v ^ 2 = c u * G.degree u)

    -- Theorem I (1) (NOTE §3), any field; P a w := (φ e_a)_w, twin G u = {u' | N u' = N u}
    theorem theoremI (hnb : ∀ v, ∃ w, G.Adj v w) (hdeg : ∀ v, (G.degree v : 𝕜) ≠ 0)
        (φ : (V → 𝕜) ≃ₗ[𝕜] (V → 𝕜))
        (hφ : ∀ x y, φ (evoMul (adjM G 𝕜) x y) = evoMul (rwM G 𝕜) (φ x) (φ y)) :
        ∃ c : V → 𝕜, (∀ a, c a ≠ 0) ∧
          (∀ a i, φ (Pi.single a 1) i ≠ 0 → evoMul (rwM G 𝕜) (φ (Pi.single a 1)) (φ (Pi.single a 1))
            = c a • ind (G.neighborFinset i)) ∧
          (∀ a, ∃ i, φ (Pi.single a 1) i ≠ 0) ∧ (∀ i, ∃ a, φ (Pi.single a 1) i ≠ 0) ∧
          (∀ a b i j, φ (Pi.single a 1) i ≠ 0 → φ (Pi.single b 1) j ≠ 0 →
            (G.neighborFinset a = G.neighborFinset b ↔ G.neighborFinset i = G.neighborFinset j) ∧
            (G.Adj a b ↔ G.Adj i j)) ∧
          (∀ a i, φ (Pi.single a 1) i ≠ 0 →
            (#(twin G a) : 𝕜) = #(twin G i) ∧ (G.degree a : 𝕜) = G.degree i) ∧
          (∀ a b, G.neighborFinset a = G.neighborFinset b →
            ∑ w, φ (Pi.single a 1) w * φ (Pi.single b 1) w = if a = b then c a * G.degree a else 0) ∧
          (∀ v u w₀ w, G.Adj v u → φ (Pi.single u 1) w₀ ≠ 0 →
            G.neighborFinset w = G.neighborFinset w₀ → ∑ u' ∈ twin G u, φ (Pi.single u' 1) w = c v)
    -- Theorem I (2), τ = id
    theorem theoremI_converse (hdeg : ∀ v, (G.degree v : 𝕜) ≠ 0) (c : V → 𝕜) (hc0 : ∀ v, c v ≠ 0)
        (P : V → V → 𝕜) (hblock : ∀ a w, P a w ≠ 0 → G.neighborFinset w = G.neighborFinset a)
        (ha : ∀ a b, G.neighborFinset a = G.neighborFinset b →
          ∑ w, P a w * P b w = if a = b then c a * G.degree a else 0)
        (hb : ∀ v u w, G.Adj v u → G.neighborFinset w = G.neighborFinset u →
          ∑ u' ∈ twin G u, P u' w = c v) :
        ∃ φ : (V → 𝕜) ≃ₗ[𝕜] (V → 𝕜), (∀ a, φ (Pi.single a 1) = P a) ∧
          ∀ x y, φ (evoMul (adjM G 𝕜) x y) = evoMul (rwM G 𝕜) (φ x) (φ y)

    -- Lemma 2 (NOTE §4)
    theorem lemma2_identity (P : ι → ι → 𝕜) (α γ : 𝕜) :
        (Fintype.card ι : 𝕜) * (γ ^ 2 - α)
          = ∑ a, ∑ b, ((∑ w, P w a * P w b) - if a = b then α else 0)
            - ∑ w, ((∑ a, P w a) - γ) * ((∑ a, P w a) + γ)
    theorem lemma2_matrix (P : Matrix ι ι 𝕜) (α γ : 𝕜) (h1 : Pᵀ * P = α • (1 : Matrix ι ι 𝕜))
        (h2 : P *ᵥ (fun _ => (1 : 𝕜)) = fun _ => γ) : (Fintype.card ι : 𝕜) * (γ ^ 2 - α) = 0
    theorem lemma2_blocks (S S' : Finset V) (P : V → V → 𝕜) (α γ : 𝕜)
        (ha : ∀ a ∈ S, ∀ b ∈ S, ∑ w ∈ S', P a w * P b w = if a = b then α else 0)
        (hb : ∀ w ∈ S', ∑ a ∈ S, P a w = γ) : (#S : 𝕜) * α = #S' * γ ^ 2

    -- Lemma 1 = BCS22 Thm 2.4(i) "⇒" (NOTE §2.2); Theorem N
    theorem lemma1 (φ : (V → 𝕜) ≃ₗ[𝕜] (V → 𝕜)) (hφ : ∀ x y, φ (evoMul M x y) = evoMul M' (φ x) (φ y))
        (a i : V) (h : φ (Pi.single a 1) i ≠ 0) :
        M' i = ((φ (Pi.single a 1) i)⁻¹ * φ.symm (Pi.single i 1) a) • φ (M a)
    theorem theoremN_mp (M : V → V → 𝕜) (hnd : ∀ i, M i ≠ 0) (ψ : (V → 𝕜) ≃ₗ[𝕜] (V → 𝕜))
        (horth : ∀ a b, a ≠ b → evoMul M (ψ (Pi.single a 1)) (ψ (Pi.single b 1)) = 0) :
        (∀ a i j, ψ (Pi.single a 1) i ≠ 0 → ψ (Pi.single a 1) j ≠ 0 →
          ∃ μ : 𝕜, μ ≠ 0 ∧ M i = μ • M j) ∧
        (∀ a, evoMul M (ψ (Pi.single a 1)) (ψ (Pi.single a 1)) ≠ 0) ∧
        (∀ (T : Finset V) (lam u : V → 𝕜), u ≠ 0 → (∀ i ∈ T, M i = lam i • u) →
          ∀ a b, a ≠ b → (∀ i, ψ (Pi.single a 1) i ≠ 0 → i ∈ T) →
            ∑ i ∈ T, lam i * ψ (Pi.single a 1) i * ψ (Pi.single b 1) i = 0) ∧
        (∀ T : Finset V, (∀ i j, i ∈ T → (∃ μ : 𝕜, M j = μ • M i) → j ∈ T) →
          ((univ.filter (fun a => ∃ i ∈ T, ψ (Pi.single a 1) i ≠ 0)).card : 𝕜) = #T)
    theorem theoremN_mpr (M : V → V → 𝕜) (ψ : (V → 𝕜) ≃ₗ[𝕜] (V → 𝕜)) (lam : V → 𝕜)
        (u : V → V → 𝕜) (hM : ∀ i, M i = lam i • u i)
        (hblock : ∀ a i j, ψ (Pi.single a 1) i ≠ 0 → ψ (Pi.single a 1) j ≠ 0 → u i = u j)
        (horth : ∀ a b, a ≠ b → ∀ i, ψ (Pi.single a 1) i ≠ 0 →
          ∑ j, lam j * ψ (Pi.single a 1) j * ψ (Pi.single b 1) j = 0) :
        ∀ a b, a ≠ b → evoMul M (ψ (Pi.single a 1)) (ψ (Pi.single b 1)) = 0

    -- Corollary A (NOTE §0)
    abbrev K12 : SimpleGraph (Fin 1 ⊕ Fin 2) := completeBipartiteGraph (Fin 1) (Fin 2)
    theorem corA_K12 : ¬ AIsoRW K12 ℚ
    theorem survey_false_over_Q : K12.Connected ∧ 2 ≤ Fintype.card (Fin 1 ⊕ Fin 2) ∧
        (∃ (s t : Set (Fin 1 ⊕ Fin 2)) (d₁ d₂ : ℕ), IsBiregular K12 s t d₁ d₂) ∧ ¬ AIsoRW K12 ℚ
    theorem corA_K13 : ¬ AIsoRW (completeBipartiteGraph (Fin 1) (Fin 3)) ℚ
    theorem corA_K23 : ¬ AIsoRW (completeBipartiteGraph (Fin 2) (Fin 3)) ℚ
    theorem K18_iso : AIsoRW (completeBipartiteGraph (Fin 1) (Fin 8)) ℚ
    theorem K2_16_iso : AIsoRW (completeBipartiteGraph (Fin 2) (Fin 16)) ℚ
    def SurveyHolds (𝕜 : Type*) [Field 𝕜] : Prop :=
      ∀ (V : Type) [Fintype V] [DecidableEq V] (G : SimpleGraph V) [DecidableRel G.Adj],
        G.Connected → 2 ≤ Fintype.card V →
          (AIsoRW G 𝕜 ↔ IsRegular G ∨ ∃ (s t : Set V) (d₁ d₂ : ℕ), IsBiregular G s t d₁ d₂)
    theorem not_surveyHolds_rat : ¬ SurveyHolds ℚ
    theorem corA_iff [CharZero 𝕜] :
        SurveyHolds 𝕜 ↔ ∀ a b : ℕ, 0 < a → 0 < b → ∃ t : 𝕜, t ^ 3 = (a : 𝕜) / b
    theorem corA_primes [CharZero 𝕜] (hp : ∀ p : ℕ, p.Prime → ∃ z : 𝕜, z ^ 3 = p) : SurveyHolds 𝕜
    theorem corA_allCubes [CharZero 𝕜] (h3 : ∀ x : 𝕜, ∃ z, z ^ 3 = x) : SurveyHolds 𝕜

## Hypotheses and sources

No named hypotheses. The one literature input of the paper proof, BCS22 Thm 2.4(i) "⇒" (Boudi–Cabrera Casado–Siles Molina, Publ. Mat. 66 (2022), arXiv 2006.14460, quoted in NOTE §2.1), is proved: `lemma1` (for an isomorphism) and `lemma1_nat` (any natural basis), by NOTE's one-line argument in coordinates (`f_a e'_i = (f_a)_i e'_i² = (φ⁻¹e'_i)_a φ(e_a²)`). Theorem N's counting `|B'_T| = |T|` is replaced by a trace identity (`card_eq_of_support`: `|A| = Σ_{a∈A}(QP)_{aa} = Σ_{w∈B}(PQ)_{ww} = |B|` when the support of `P` respects `A × B`), which gives `k_{τS} = k_S` and `d_{τS} = d_S` in `𝕜`, all Lemma 2 needs. Mathlib theorems used: SimpleGraph API (`reachable_iff_reflTransGen`, `Reachable.nonempty_neighborSet_left`, `degree_pos_iff_exists_adj`, `IsBipartiteWith.completeBipartiteGraph`), `mul_eq_one_comm` for square matrices over a commutative ring (`Matrix.SemiringInverse`; Theorem I (2) only), `Rat.num_pow` (core).

## Build

    bash Riemann/Deep/EvoAlg/Lean/audit/buildall.sh     # clean rebuild, one module at a time
    # each module: LEAN_PATH=a scratch olean directory (build output, deleted):<lake LEAN_PATH> \
    #   bash Lean/audit/memcap.sh 290 1800 lean -j1 -o a scratch olean (build output, deleted) M.lean
    order: Defs Basic Lemma2 Diag TwinI Cor41 TheoremA CorA Examples CorA2 CorA3 TheoremN TheoremI2 Lemma2M Check

Clean rebuild 11:29–11:36 UTC, all rc = 0, watchdog peaks (MB): Defs 207, Basic 172, Lemma2 249, Diag 208, TwinI 242, Cor41 256, TheoremA 212, CorA 214, Examples 254, CorA2 254, TheoremN 208, TheoremI2 226, Lemma2M 218, Check 261; CorA3 207 and Check (with CorA3) 205 at 11:40–11:41 (audit/peaks.log). Only linter warnings ("automatically included section variable(s) unused", 13). No `sorry`, `admit`, `axiom`, `native_decide`, `skipKernelTC`, `implemented_by` (rg); `decide` only on degrees and vertex counts of K_{1,2}, K_{1,3}, K_{2,3}, K_{1,8}. Watchdog kills: TheoremI2 with `Matrix.NonsingularInverse` (315, 303 MB) → replaced by `Matrix.SemiringInverse` (226 MB); import probe `Mathlib.FieldTheory.IsAlgClosed.Basic` alone 314 MB (hence no ℂ instance). Check.lean output (audit/Check.log), identical for all 31 printed theorems:

    'EvoAlg.theoremA' depends on axioms: [propext, Classical.choice, Quot.sound]
    ... (theoremA_if, sec5, cor41, theoremI, theoremI_converse, scaling_ii, lemma1, lemma1_nat, theoremN_mp,
    theoremN_mpr, lemma2, lemma2_identity, lemma2_blocks, lemma2_matrix, card_eq_of_support, aIsoRW_of_regular,
    aIsoRW_of_biregular, corA_K12, survey_false_over_Q, not_surveyHolds_rat, corA_cubes, corA_necessity, corA_iff,
    corA_allCubes, corA_primes, corA_K13, corA_K23, K18_iso, K2_16_iso, not_cube_two): [propext, Classical.choice, Quot.sound]

## Fidelity

- Graphs: Mathlib `SimpleGraph` on a `Fintype` (finite, no loops, no multiple edges), `G.Connected`, `2 ≤ Fintype.card V`: NOTE §0 "G is finite and connected with n ≥ 2 vertices", §1 "twins are non-adjacent, since v ∉ N(v)". Survey: "for every finite connected graph G" (2610.22248 TeX l. 1352). CRR's graph definition (1710.10516 l. 159, `E := {(i,j) ∈ V × V : i ≤ j}`) formally admits loops; loops are outside the NOTE's and this formalization's scope.
- Evolution algebra (survey l. 107–109: a `𝕜`-algebra with a basis `(e_i)`, `e_ie_j = 0` for `i ≠ j`; `𝕜` an arbitrary field): structure-constant model `(𝕜^V, evoMul M)`, `e_v = Pi.single v 1`; `evoMul_single_single` shows `e_ve_u = 0` (`v ≠ u`), `e_v² = M v`; `evoMul` is the bilinear commutative extension (`evoMul_comm`, `evoMul_add_left`, `evoMul_smul_left`). Every evolution algebra with natural basis indexed by `V` is isomorphic to this one by its coordinate map, so nothing is lost.
- A(G), A_RW(G) (NOTE §0; survey l. 617–619: structure matrices "the adjacency matrix" and "the transition matrix of the symmetric random walk"): `adjM`, `rwM` with `adj_sq`: `e_v² = Σ_{w∈N(v)} e_w` and `rw_sq`: `e'_v² = (deg v)⁻¹ Σ_{w∈N(v)} e'_w`. `(deg v)⁻¹` is Lean's `0⁻¹ = 0` convention when `deg v = 0` in `𝕜`; never used: Theorem A has char 0 and every vertex has a neighbour; `theoremI`, `cor41` assume `hdeg` (NOTE: "d_v must be invertible in 𝕜").
- ≅: `EvoIso` = a `𝕜`-linear equivalence `φ` with `φ(xy) = φ(x)φ(y)` (brief: "linear equivalence preserving the product").
- Regular: `∃ d, G.IsRegularOfDegree d` (Mathlib: every degree is `d`).
- Biregular (CRR l. 159: "a bipartite graph G = (V₁, V₂, E) for which every two vertices on the same side of the given bipartition have the same degree"; bipartite: "vertices can be divided into two disjoint sets V₁ and V₂ such that every edge connects a vertex in V₁ to one in V₂"): `IsBiregular G s t d₁ d₂` = Mathlib `IsBipartiteWith s t` (disjoint, every edge crosses) + `s ∪ t = V` + constant degrees `d₁` on `s`, `d₂` on `t`.
- Theorem A's condition 2, "d₁/d₂ ∈ (𝕜^×)³": `∃ c : 𝕜, c ≠ 0 ∧ (d₁ : 𝕜) / d₂ = c ^ 3` (symmetric in the two sides: `d₂/d₁ = (c⁻¹)³`).
- Twin quotient H, τ, blocks P_S (NOTE §1, §3): vertex-level. Classes `twin G u`; `τ(S) ∋ i` iff `(φ e_a)_i ≠ 0` for some `a ∈ S`; "τ is a bijection of classes" = columns and rows non-zero + `a ≡ b ↔ i ≡ j`; "τ ∈ Aut H" = `a ∼ b ↔ i ∼ j`; (a), (b) as sums over the blocks; `k_{τS} = k_S`, `d_{τS} = d_S` in `𝕜` (equal in ℕ in char 0). `c : V → 𝕜` is constant on classes (forced by (b)), i.e. a function on `V(H)`.
- Corollary 4.1 (i) "c_T = c_{T'} whenever T, T' have a common H-neighbour" = `∀ u v v', u ∼ v → u ∼ v' → c v = c v'`; (ii) "c_T² = c_S d_S whenever S ∼ T" = `u ∼ v → c v² = c u d_u` (`u ∈ S`, `v ∈ T`). Hypothesis: every twin-class size non-zero in `𝕜` (char 0: automatic).
- Lemma 2: square form with `ι`-indexed `P` (`P w a` = row `w`, column `a`) and with Mathlib `Matrix` (`Pᵀ * P = α • 1`, `P *ᵥ 1 = γ • 1`); the NOTE's identity is `lemma2_identity` verbatim.
- Theorem N (NOTE §2.2) in coordinates: a natural basis indexed by `V` is `f_a = ψ e_a` (`ψ` a linear equivalence: `(f_a)` a basis) with `f_af_b = 0`; "⇒" gives `f_a ∈ E_T`, `f_a² ≠ 0`, `b_T`-orthogonality and `|B'_T| = |T|` in `𝕜` (in ℕ for char 0; char p: congruence only); "⇐" gives naturality. Parts (a)–(c) of Theorem N are not formalized.
- Corollary A: `SurveyHolds 𝕜` is the survey's Question q:rw (l. 1351–1354) over `𝕜`, vertex types in `Type`. Over ℚ: `not_surveyHolds_rat` (K_{1,2}); K_{1,3}, K_{2,3} non-isomorphic, K_{1,8}, K_{2,16} isomorphic. `corA_iff` = NOTE's second item; `corA_primes` = first item for "cube root of every prime"; ℝ, ℂ: `corA_allCubes` with the field condition "every element is a cube" as hypothesis, not discharged in Lean (cap).

Not formalized: the characteristic-2 clauses of Lemma 2 and Corollary 4.1, Theorem P and Lemma 3 (char p), Proposition H, Remark 3.2, NOTE §8, the instances ℝ, ℂ.

## Reviewer (statement fidelity, verdict-only, statements pasted inline)

Three passes ([audit/review_prompt3.txt](Lean/audit/review_prompt3.txt) = NOTE §0–§5 statements + all Lean statements above; verdicts in audit/review1–3.txt):
- Pass 1: ISSUES — (blocking) necessity half of Corollary A missing; (minor) no ℝ/ℂ specialization; (minor) K_{2,16} missing. Fixed: CorA2 (`corA_necessity`, `corA_iff`, `corA_allCubes`, `K2_16_iso`).
- Pass 2: ISSUES — (minor) ℝ/ℂ conditional; (minor) no bridge from "cube root of every prime". Fixed: CorA3 (`cube_of_primes`, `corA_primes`).
- Pass 3: **FAITHFUL** — "Minor certification gap only: unconditional ℝ/ℂ specializations remain uncompiled. `corA_allCubes` states the correct abstract implication; `corA_primes` now covers the prime-root hypothesis. No remaining statement-fidelity defect found."

## FORMALIZATIONS.md entry

- EvoAlg (our discovery [evolution-algebras](discoveries/evolution-algebras.md)) machine-checked in Lean 4 + Mathlib: EvoAlg.theoremA — for char 𝕜 = 0 and G finite connected simple with ≥ 2 vertices, A(G) ≅ A_RW(G) ⟺ G regular or biregular with d₁/d₂ ∈ (𝕜^×)³, both directions; with Theorem I (1)+(2), Lemma 2, Corollary 4.1, Theorem N, and Corollary A (false over ℚ: K_{1,2}, K_{1,3}, K_{2,3}; holds over 𝕜 ⟺ every positive rational is a cube). No named hypotheses (BCS22 Thm 2.4(i) proved), sorry-free, standard axioms, 1,725 lines; ℝ/ℂ instances not compiled (290 MB cap). [LEAN.md](Riemann/Deep/EvoAlg/LEAN.md)

