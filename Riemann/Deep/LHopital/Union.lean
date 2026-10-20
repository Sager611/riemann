import Riemann.Fayed.Theorem1

/-!
# A correct L'Hôpital rule on a union of intervals

Fayed's Appendix Theorem 1 (formalized and refuted in \`Riemann/Fayed/Theorem1.lean\`)
claims a L'Hôpital rule along a disconnected domain \`D = ⋃ₘ (a M, b M)\`. The honest
version proved here: if \`F'/G' → L\` along \`D\`, the left-endpoint ratios
\`F (a M) / G (a M) → L\`, and \`G (a M) / G x\` stays bounded on each block, then
\`F/G → L\` along \`D\`. Proof: Cauchy's mean value theorem on \`[a M, x]\` gives
\`(F x - F (a M)) / (G x - G (a M)) = F' ξ / G' ξ\` with \`ξ ∈ (a M, x) ⊆ D\`, and
\`F x / G x = (F(a)/G(a))·(G(a)/G(x)) + (F' ξ/G' ξ)·(1 - G(a)/G(x))\`.

Note on the paper: Fayed's argument (step (55)) actually needs the CONVERSE direction,
\`F/G → L ⇒ F'/G' → L\`. That converse is false even on a single ray \`(c, ∞)\`
(e.g. \`F x = sin x\`, \`G x = x\`: \`F/G → 0\` but \`F'/G' = cos x\` has no limit), so no
hypotheses on the gaps between blocks can rescue it. This theorem shows what a
disconnected-domain L'Hôpital rule can honestly give: the forward direction only, and
only with the block-endpoint data \`hend\` and the uniform bound \`hratio\` that tie the
blocks together.
-/

open Filter Topology Set

namespace Fayed.LHopital

theorem lhopital_union {F G F' G' : ℝ → ℝ} {a b : ℕ → ℝ} {L : ℝ}
    (hab : ∀ M, a M < b M) (hbound : ∀ M, b M ≤ a (M+1)) (ha : Tendsto a atTop atTop)
    (hF : ∀ M, ∀ x ∈ Set.Ioo (a M) (b M), HasDerivAt F (F' x) x)
    (hG : ∀ M, ∀ x ∈ Set.Ioo (a M) (b M), HasDerivAt G (G' x) x)
    (hFc : ∀ M, ContinuousOn F (Set.Icc (a M) (b M)))
    (hGc : ∀ M, ContinuousOn G (Set.Icc (a M) (b M)))
    (hG'0 : ∀ M, ∀ x ∈ Set.Ioo (a M) (b M), G' x ≠ 0)
    (hGpos : ∀ M, ∀ x ∈ Set.Icc (a M) (b M), 0 < G x)
    (hratio : ∃ K, ∀ M, ∀ x ∈ Set.Icc (a M) (b M), G (a M) / G x ≤ K)
    (hend : Tendsto (fun M => F (a M) / G (a M)) atTop (nhds L))
    (hder : Tendsto (fun x => F' x / G' x) (along (⋃ M, Set.Ioo (a M) (b M))) (nhds L)) :
    Tendsto (fun x => F x / G x) (along (⋃ M, Set.Ioo (a M) (b M))) (nhds L) := by
  obtain ⟨K, hK⟩ := hratio
  set D := ⋃ M, Set.Ioo (a M) (b M) with hD
  have hmono : StrictMono a := strictMono_nat_of_lt_succ fun M => (hab M).trans_le (hbound M)
  rw [Metric.tendsto_nhds]
  intro ε hε
  set C : ℝ := |K| + 1 with hC
  have hCpos : 0 < C := by positivity
  set η : ℝ := ε / (2 * C) with hη
  have hηpos : 0 < η := by positivity
  obtain ⟨N1, hN1⟩ := eventually_atTop.1 (Metric.tendsto_nhds.1 hend η hηpos)
  have h2 := Metric.tendsto_nhds.1 hder η hηpos
  rw [along, eventually_inf_principal] at h2
  obtain ⟨X, hX⟩ := eventually_atTop.1 h2
  obtain ⟨N2, hN2⟩ := eventually_atTop.1 (tendsto_atTop.1 ha X)
  set N := max N1 N2
  rw [along, eventually_inf_principal]
  filter_upwards [eventually_ge_atTop (a (N + 1))] with x hx hxD
  obtain ⟨M, hxM⟩ := mem_iUnion.1 hxD
  have hMN : N < M := by
    by_contra h
    have h' : M + 1 ≤ N + 1 := by omega
    have : b M ≤ a (N + 1) := (hbound M).trans (hmono.monotone h')
    linarith [hxM.2]
  have hM1 : N1 ≤ M := by omega
  have hM2 : N2 ≤ M := by omega
  have hax : a M < x := hxM.1
  have hsubI : Icc (a M) x ⊆ Icc (a M) (b M) := Icc_subset_Icc le_rfl hxM.2.le
  have hsubO : Ioo (a M) x ⊆ Ioo (a M) (b M) := Ioo_subset_Ioo le_rfl hxM.2.le
  obtain ⟨ξ, hξ, hξeq⟩ := exists_ratio_hasDerivAt_eq_ratio_slope F F' hax ((hFc M).mono hsubI)
    (fun y hy => hF M y (hsubO hy)) G G' ((hGc M).mono hsubI) (fun y hy => hG M y (hsubO hy))
  obtain ⟨θ, hθ, hθeq⟩ := exists_hasDerivAt_eq_slope G G' hax ((hGc M).mono hsubI)
    (fun y hy => hG M y (hsubO hy))
  have hGθ : G' θ ≠ 0 := hG'0 M θ (hsubO hθ)
  have hGξ : G' ξ ≠ 0 := hG'0 M ξ (hsubO hξ)
  have hdiff : G x - G (a M) ≠ 0 := by
    intro h0
    rw [h0, zero_div] at hθeq
    exact hGθ hθeq
  have hGx : 0 < G x := hGpos M x ⟨hax.le, hxM.2.le⟩
  have hGa : 0 < G (a M) := hGpos M (a M) ⟨le_rfl, (hab M).le⟩
  set t := G (a M) / G x with ht
  have htpos : 0 < t := div_pos hGa hGx
  have htK : t ≤ K := hK M x ⟨hax.le, hxM.2.le⟩
  have hr : F' ξ / G' ξ = (F x - F (a M)) / (G x - G (a M)) := by
    rw [div_eq_div_iff hGξ hdiff]; linarith [hξeq]
  have hident : F x / G x - L
      = (F (a M) / G (a M) - L) * t + (F' ξ / G' ξ - L) * (1 - t) := by
    rw [hr, ht]
    field_simp
    ring
  have e1 : |F (a M) / G (a M) - L| < η := by
    have := hN1 M hM1; rwa [Real.dist_eq] at this
  have e2 : |F' ξ / G' ξ - L| < η := by
    have hξX : X ≤ ξ := (hN2 M hM2).trans hξ.1.le
    have := hX ξ hξX (mem_iUnion.2 ⟨M, hsubO hξ⟩)
    rwa [Real.dist_eq] at this
  have ht1 : |t| ≤ C := by rw [abs_of_pos htpos]; linarith [le_abs_self K]
  have ht2 : |1 - t| ≤ C := by
    rw [abs_le]; constructor <;> linarith [le_abs_self K]
  rw [Real.dist_eq, hident]
  calc |(F (a M) / G (a M) - L) * t + (F' ξ / G' ξ - L) * (1 - t)|
      ≤ |F (a M) / G (a M) - L| * |t| + |F' ξ / G' ξ - L| * |1 - t| := by
        rw [← abs_mul, ← abs_mul]; exact abs_add_le _ _
    _ ≤ |F (a M) / G (a M) - L| * C + |F' ξ / G' ξ - L| * C := by
        gcongr
    _ < η * C + η * C := by
        gcongr
    _ = ε := by rw [hη]; field_simp; ring

end Fayed.LHopital

#print axioms Fayed.LHopital.lhopital_union
