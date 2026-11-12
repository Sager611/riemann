import Mathlib.Basic.Real.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

/-!
# The certificate's A: piecewise linear on 78 rational nodes, constant tail (DERIVATION §5)

plin h L y₀: linear interpolation through the values L at y₀, y₀ + h, …, constant (= the last value) beyond
the last node. Acert = plin dg nodesA λ with λ = 96/125, dg = 192/9625 = (U − λ)/77, nodesA = certs/cert768.json
"A" (the values cert_aug.py interpolates, l.98–103; A(x) = A(U) for x ≥ U, l.128–129). Proved: Acert = 0 on
[U, ∞) (last node 0) and A_0 ≤ Acert ≤ 0 on [λ, ∞) (all 78 nodes in [A_0, 0], by norm_num on ℚ).
-/

namespace PPLean

/-- Piecewise-linear interpolation through the values L at y₀, y₀ + h, y₀ + 2h, …; the last value beyond. -/
noncomputable def plin (h : ℝ) : List ℚ → ℝ → ℝ → ℝ
  | [], _, _ => 0
  | [a], _, _ => (a : ℝ)
  | a :: b :: rest, y₀, y =>
      if y ≤ y₀ + h then (a : ℝ) + ((b : ℝ) - a) * ((y - y₀) / h) else plin h (b :: rest) (y₀ + h) y

/-- The last element (0 for the empty list). -/
def lastQ : List ℚ → ℚ
  | [] => 0
  | [a] => a
  | _ :: b :: rest => lastQ (b :: rest)

/-- Beyond the last node, plin is the last value. -/
lemma plin_tail {h : ℝ} (hh : 0 < h) : ∀ (L : List ℚ) (y₀ y : ℝ),
    y₀ + ((L.length : ℝ) - 1) * h ≤ y → plin h L y₀ y = (lastQ L : ℝ)
  | [], _, _, _ => by simp [plin, lastQ]
  | [a], _, _, _ => by simp [plin, lastQ]
  | a :: b :: rest, y₀, y, hy => by
    have hlen : ((a :: b :: rest).length : ℝ) = (rest.length : ℝ) + 2 := by
      simp only [List.length_cons]; push_cast; ring
    rw [hlen] at hy
    have hr0 : (0 : ℝ) ≤ (rest.length : ℝ) := Nat.cast_nonneg _
    simp only [plin, lastQ]
    split_ifs with h1
    · rcases rest with _ | ⟨c, rest'⟩
      · have hy' : y = y₀ + h := by simp at hy; linarith
        rw [hy', show (y₀ + h - y₀) / h = 1 by rw [show y₀ + h - y₀ = h by ring]; exact div_self hh.ne']
        simp [lastQ]
      · simp only [List.length_cons] at hy
        push_cast at hy
        have : (0 : ℝ) ≤ (rest'.length : ℝ) := Nat.cast_nonneg _
        nlinarith
    · have := plin_tail hh (b :: rest) (y₀ + h) y (by simp only [List.length_cons]; push_cast; linarith)
      rw [this]

/-- plin stays between the least and the largest node value (from the first node on). -/
lemma plin_mem {h m M : ℝ} (hh : 0 < h) (hm : m ≤ 0) (hM : 0 ≤ M) : ∀ (L : List ℚ) (y₀ y : ℝ),
    y₀ ≤ y → (∀ q ∈ L, m ≤ (q : ℝ) ∧ (q : ℝ) ≤ M) → m ≤ plin h L y₀ y ∧ plin h L y₀ y ≤ M
  | [], _, _, _, _ => by simp only [plin]; exact ⟨hm, hM⟩
  | [a], _, _, _, hL => by simp only [plin]; exact hL a (by simp)
  | a :: b :: rest, y₀, y, hy, hL => by
    simp only [plin]
    split_ifs with h1
    · have ha := hL a (by simp)
      have hb := hL b (by simp)
      have ht0 : 0 ≤ (y - y₀) / h := div_nonneg (by linarith) hh.le
      have ht1 : (y - y₀) / h ≤ 1 := by rw [div_le_one hh]; linarith
      generalize (y - y₀) / h = t at ht0 ht1 ⊢
      have e1 := mul_nonneg ht0 (by linarith [hb.1] : (0 : ℝ) ≤ (b : ℝ) - m)
      have e2 := mul_nonneg (by linarith : (0 : ℝ) ≤ 1 - t) (by linarith [ha.1] : (0 : ℝ) ≤ (a : ℝ) - m)
      have e3 := mul_nonneg ht0 (by linarith [hb.2] : (0 : ℝ) ≤ M - (b : ℝ))
      have e4 := mul_nonneg (by linarith : (0 : ℝ) ≤ 1 - t) (by linarith [ha.2] : (0 : ℝ) ≤ M - (a : ℝ))
      constructor
      · linarith
      · linarith
    · exact plin_mem hh hm hM (b :: rest) (y₀ + h) y (by linarith)
        (fun q hq => hL q (List.mem_cons_of_mem a hq))

/-- certs/cert768.json "A": the values of A at λ + k·dg, k = 0, …, 77. -/
def nodesA : List ℚ := [
    -378347742845/758342134484,
    -467732272337/969851355419,
    -216489104093/480910907194,
    -282065550921/627980265532,
    -223405774904/548019355385,
    -348165969697/896092847827,
    -58935669735/158981887277,
    -353460938441/918671047283,
    -28261456203/83474548727,
    -236319896615/667062578209,
    -238587194441/768007019913,
    -915707451/2811421249,
    -218612182274/699882247441,
    -123223408786/411131200267,
    -34357881072/128962023479,
    -250580553752/972549497515,
    -241962229347/968816139922,
    -92873173291/360931649823,
    -50279191574/202036539767,
    -185186679337/767574790709,
    -9304902305/41221918313,
    -201608969933/910940066850,
    -186491187977/834992391745,
    -48654150523/222266368629,
    -209762217555/974885547428,
    -152622374877/719545868807,
    -169561343344/819795398219,
    -72894961359/350638849568,
    -143490684009/701878607360,
    -68743423491/337016934589,
    -201037854241/985544582789,
    -96237439141/470743985768,
    -203478587544/985198814717,
    -90590505345/438787932763,
    -3705471884/17823433261,
    -208655358549/993290183096,
    -207910962211/982353214825,
    -151328135427/709306191631,
    -100953251855/469278687679,
    -143140734363/659815812526,
    -154734076677/707337391522,
    -216347236282/980901298391,
    -207567399047/933689147312,
    -85563592837/382001592169,
    -71013332508/314809399561,
    -221147655691/974039006564,
    -205367208122/899126526921,
    -187984969131/819363059890,
    -183791145996/796777317509,
    -201258815638/869194365031,
    -89707671332/386271828885,
    -134207098612/576642286259,
    -121020503143/519316352923,
    -226143095792/993305112231,
    -78660360855/337547654992,
    -119571645074/513749184817,
    -57066159671/245685126341,
    -214260048761/998381388599,
    -21784063841/94355634050,
    -100765690925/438141332737,
    -160921913587/702696691686,
    -160427983846/703743094645,
    -213843966111/942526362323,
    -68973667994/305475142627,
    -153656868451/683800609118,
    -112627477923/725848285558,
    -210652547973/946074245663,
    -194595221960/877654226651,
    -206362075801/934366958921,
    -128428344239/583570369546,
    -94158538071/968777612915,
    -34746129332/158834592113,
    -184774798411/846725547271,
    -35347850542/622466395889,
    -39309874364/920420722357,
    -27632983727/972277599383,
    -13841655406/990364696389,
    (0 : ℚ)]

lemma nodesA_length : nodesA.length = 78 := rfl

lemma lastQ_nodesA : lastQ nodesA = 0 := rfl

lemma nodesA_bounds : ∀ q ∈ nodesA, (-378347742845 / 758342134484 : ℚ) ≤ q ∧ q ≤ 0 := by
  simp only [nodesA, List.forall_mem_cons]
  norm_num

/-- The certificate's A. -/
noncomputable def Acert (y : ℝ) : ℝ := plin (192 / 9625) nodesA (96 / 125) y

/-- A(x) = A(U) = 0 for x ≥ U = 288/125 = λ + 77·dg. -/
lemma Acert_tail (y : ℝ) (hy : 288 / 125 ≤ y) : Acert y = 0 := by
  unfold Acert
  rw [plin_tail (by norm_num) nodesA (96 / 125) y (by rw [nodesA_length]; norm_num; linarith),
    lastQ_nodesA]
  simp

/-- A_0 ≤ A ≤ 0 on [λ, ∞). -/
lemma Acert_range (y : ℝ) (hy : 96 / 125 ≤ y) :
    (-378347742845 / 758342134484 : ℝ) ≤ Acert y ∧ Acert y ≤ 0 := by
  have h := plin_mem (h := 192 / 9625) (m := ((-378347742845 / 758342134484 : ℚ) : ℝ)) (M := 0)
    (by norm_num) (by norm_num) le_rfl nodesA (96 / 125) y hy (fun q hq => by
      have := nodesA_bounds q hq
      exact ⟨by exact_mod_cast this.1, by exact_mod_cast this.2⟩)
  unfold Acert
  push_cast at h
  exact h

end PPLean
