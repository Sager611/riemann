import Riemann.Deep.JerbyFals.Lean.Defs

/-!
# Jerby's Gram-point repulsion conjecture is false

The named hypotheses `certAL` … `certCR` are the Arb outputs of certify_cx.py (python-flint,
256 bits; logs/certify_cx.log, certify_cx2.log, certify_cx3.log): an interval-Newton enclosure of
the Gram point g_k (θ(g_k) = kπ), and enclosures of Z(g_k) and, for k = n, of Z′(g_n), as computed
there (`hardyZ'`). Each interval is the printed Arb ball [mid ± rad] (lower end mid − rad,
upper end mid + rad). Everything else (parity, good/bad, |Z′| < 4|Z|, uniqueness and
consecutiveness of the Gram points under `ThetaDerivPos`) is proved here.
-/

open Set
open scoped Real

namespace JerbyLean

/-- Arb record: a Gram point of index k in [a, b] at which Z ∈ [zl, zu]. -/
def GramCert (k : ℤ) (a b zl zu : ℝ) : Prop :=
  ∃ t ∈ Icc a b, theta t = k * π ∧ hardyZ t ∈ Icc zl zu

/-- Arb record with the derivative: also Z′ ∈ [dl, du] (Z′ = `hardyZ'`, Arb's formula). -/
def GramCertD (k : ℤ) (a b zl zu dl du : ℝ) : Prop :=
  ∃ t ∈ Icc a b, theta t = k * π ∧ hardyZ t ∈ Icc zl zu ∧ hardyZ' t ∈ Icc dl du

/-- With θ increasing, Gram points of index n − 1, n, n + 1 are unique and consecutive. -/
theorem unique_consecutive (hθ : ThetaDerivPos) {n : ℤ} {gl g gr : ℝ}
    (hl : IsGram (n - 1) gl) (hg : IsGram n g) (hr : IsGram (n + 1) gr) :
    gl < g ∧ g < gr ∧ (∀ m t, IsGram m t → gl < t → t < gr → m = n ∧ t = g) ∧
      (∀ t, IsGram (n - 1) t → t = gl) ∧ (∀ t, IsGram (n + 1) t → t = gr) := by
  obtain ⟨h1, h2, h3⟩ := consecutive hθ hl hg hr
  exact ⟨h1, h2, h3, fun t ht => isGram_unique hθ ht hl, fun t ht => isGram_unique hθ ht hr⟩

/-! ## n = 1,565,734,083 (odd; logs/certify_cx.log) -/

/-- logs/certify_cx.log line 2: k = n − 1. -/
def certAL : Prop := GramCert 1565734082
  568009639.356005286611718227430573 568009639.356005286611718227431427
  0.27634005895678409472 0.27634005895678410528

/-- logs/certify_cx.log lines 3 and 6: k = n, and Z′(g_n). -/
def certAM : Prop := GramCertD 1565734083
  568009639.698978187581685032968666 568009639.698978187581685032969334
  18.39714778558963574 18.39714778558964426
  73.00385040031973562 73.00385040031974438

/-- logs/certify_cx.log line 4: k = n + 1. -/
def certAR : Prop := GramCert 1565734084
  568009640.041951088540347534043632 568009640.041951088540347534044368
  4.076959407739047812 4.076959407739048188

theorem counterexampleA (hL : certAL) (hM : certAM) (hR : certAR) :
    ∃ gl g gr : ℝ, IsGram (1565734083 - 1) gl ∧ IsGram 1565734083 g ∧
      IsGram (1565734083 + 1) gr ∧ Good (1565734083 - 1) gl ∧
      (-1 : ℝ) ^ (1565734083 : ℤ) * hardyZ g < 0 ∧ Good (1565734083 + 1) gr ∧
      |deriv hardyZ g| < 4 * |hardyZ g| := by
  obtain ⟨gl, ⟨hl1, -⟩, hlθ, hlz, -⟩ := hL
  obtain ⟨g, ⟨hg1, -⟩, hgθ, ⟨hgz, -⟩, hgd1, hgd2⟩ := hM
  obtain ⟨gr, ⟨hr1, -⟩, hrθ, hrz, -⟩ := hR
  have hodd : Odd (1565734083 : ℤ) := ⟨782867041, by norm_num⟩
  have hel : Even ((1565734083 : ℤ) - 1) := ⟨782867041, by norm_num⟩
  have her : Even ((1565734083 : ℤ) + 1) := ⟨782867042, by norm_num⟩
  have hz : 0 < hardyZ g := lt_of_lt_of_le (by norm_num) hgz
  have hd : 0 < hardyZ' g := lt_of_lt_of_le (by norm_num) hgd1
  refine ⟨gl, g, gr, ⟨le_trans (by norm_num) hl1, by rw [hlθ]; norm_num⟩,
    ⟨le_trans (by norm_num) hg1, hgθ⟩, ⟨le_trans (by norm_num) hr1, by rw [hrθ]; norm_num⟩,
    good_of_even_pos hel (lt_of_lt_of_le (by norm_num) hlz), gramLaw_neg_of_odd_pos hodd hz,
    good_of_even_pos her (lt_of_lt_of_le (by norm_num) hrz), ?_⟩
  rw [deriv_hardyZ, abs_of_pos hd, abs_of_pos hz]
  have key : (73.00385040031974438 : ℝ) < 4 * 18.39714778558963574 := by norm_num
  linarith

/-- Jerby's conjecture, neighbours given, is false (no hypothesis on θ′). -/
theorem not_jerbyConjEx (hL : certAL) (hM : certAM) (hR : certAR) : ¬ JerbyConjEx := by
  obtain ⟨gl, g, gr, h1, h2, h3, h4, h5, h6, h7⟩ := counterexampleA hL hM hR
  exact not_jerbyConjEx_of h1 h2 h3 h4 (bad_of_gramLaw_neg h5) h6 h7.le

/-- **Main theorem.** Jerby's conjecture (2406.06548 conj:rep) is false, given θ′ > 0 on [7, ∞)
and the Arb certificates of n = 1,565,734,083. -/
theorem not_jerbyConj (hθ : ThetaDerivPos) (hL : certAL) (hM : certAM) (hR : certAR) :
    ¬ JerbyConj :=
  fun H => not_jerbyConjEx hL hM hR (jerbyConj_imp_ex hθ H)

/-- The full picture at n = 1,565,734,083: g_{n−1} < g_n < g_{n+1} are the unique Gram points of
their indices, consecutive, good–bad–good, and |Z′(g_n)| < 4|Z(g_n)|. -/
theorem counterexampleA_full (hθ : ThetaDerivPos) (hL : certAL) (hM : certAM) (hR : certAR) :
    ∃ gl g gr : ℝ, IsGram (1565734083 - 1) gl ∧ IsGram 1565734083 g ∧
      IsGram (1565734083 + 1) gr ∧ gl < g ∧ g < gr ∧
      (∀ m t, IsGram m t → gl < t → t < gr → m = 1565734083 ∧ t = g) ∧
      (∀ t, IsGram (1565734083 - 1) t → t = gl) ∧ (∀ t, IsGram (1565734083 + 1) t → t = gr) ∧
      Good (1565734083 - 1) gl ∧ Bad 1565734083 g ∧ Good (1565734083 + 1) gr ∧
      |deriv hardyZ g| < 4 * |hardyZ g| := by
  obtain ⟨gl, g, gr, h1, h2, h3, h4, h5, h6, h7⟩ := counterexampleA hL hM hR
  obtain ⟨c1, c2, c3, c4, c5⟩ := unique_consecutive hθ h1 h2 h3
  exact ⟨gl, g, gr, h1, h2, h3, c1, c2, c3, c4, c5, h4, bad_of_gramLaw_neg h5, h6, h7⟩

/-! ## n = 2,565,007,238 (even; logs/certify_cx2.log) -/

/-- logs/certify_cx2.log line 2: k = n − 1. -/
def certBL : Prop := GramCert 2565007237
  906089969.418246515681814579743688 906089969.418246515681814579744312
  (-4.26922836496586012) (-4.26922836496585988)

/-- logs/certify_cx2.log lines 3 and 6: k = n, and Z′(g_n). -/
def certBM : Prop := GramCertD 2565007238
  906089969.752693821200280999796854 906089969.752693821200280999797146
  (-13.98639483953291117) (-13.98639483953290883)
  55.17454072188090567 55.17454072188091433

/-- logs/certify_cx2.log line 4: k = n + 1. -/
def certBR : Prop := GramCert 2565007239
  906089970.087141126712176413105991 906089970.087141126712176413106009
  (-0.3215394489475987333) (-0.3215394489475986667)

theorem counterexampleB (hL : certBL) (hM : certBM) (hR : certBR) :
    ∃ gl g gr : ℝ, IsGram (2565007238 - 1) gl ∧ IsGram 2565007238 g ∧
      IsGram (2565007238 + 1) gr ∧ Good (2565007238 - 1) gl ∧
      (-1 : ℝ) ^ (2565007238 : ℤ) * hardyZ g < 0 ∧ Good (2565007238 + 1) gr ∧
      |deriv hardyZ g| < 4 * |hardyZ g| := by
  obtain ⟨gl, ⟨hl1, -⟩, hlθ, -, hlz⟩ := hL
  obtain ⟨g, ⟨hg1, -⟩, hgθ, ⟨-, hgz⟩, hgd1, hgd2⟩ := hM
  obtain ⟨gr, ⟨hr1, -⟩, hrθ, -, hrz⟩ := hR
  have hev : Even (2565007238 : ℤ) := ⟨1282503619, by norm_num⟩
  have hol : Odd ((2565007238 : ℤ) - 1) := ⟨1282503618, by norm_num⟩
  have hor : Odd ((2565007238 : ℤ) + 1) := ⟨1282503619, by norm_num⟩
  have hz : hardyZ g < 0 := lt_of_le_of_lt hgz (by norm_num)
  have hd : 0 < hardyZ' g := lt_of_lt_of_le (by norm_num) hgd1
  refine ⟨gl, g, gr, ⟨le_trans (by norm_num) hl1, by rw [hlθ]; norm_num⟩,
    ⟨le_trans (by norm_num) hg1, hgθ⟩, ⟨le_trans (by norm_num) hr1, by rw [hrθ]; norm_num⟩,
    good_of_odd_neg hol (lt_of_le_of_lt hlz (by norm_num)), gramLaw_neg_of_even_neg hev hz,
    good_of_odd_neg hor (lt_of_le_of_lt hrz (by norm_num)), ?_⟩
  rw [deriv_hardyZ, abs_of_pos hd, abs_of_neg hz]
  have key : (55.17454072188091433 : ℝ) < 4 * 13.98639483953290883 := by norm_num
  linarith

theorem not_jerbyConj_B (hθ : ThetaDerivPos) (hL : certBL) (hM : certBM) (hR : certBR) :
    ¬ JerbyConj := fun H => by
  obtain ⟨gl, g, gr, h1, h2, h3, h4, h5, h6, h7⟩ := counterexampleB hL hM hR
  exact not_jerbyConjEx_of h1 h2 h3 h4 (bad_of_gramLaw_neg h5) h6 h7.le (jerbyConj_imp_ex hθ H)

/-! ## n = 3,315,629,581 (odd; logs/certify_cx3.log) -/

/-- logs/certify_cx3.log line 2: k = n − 1. -/
def certCL : Prop := GramCert 3315629580
  1155454822.77741041408041789235762 1155454822.77741041408041789236238
  4.066892916417592641 4.066892916417593359

/-- logs/certify_cx3.log lines 3 and 6: k = n, and Z′(g_n). -/
def certCM : Prop := GramCertD 3315629581
  1155454823.10758508544277960702679 1155454823.10758508544277960703321
  12.80566445995212596 12.80566445995213404
  (-50.235085620125960222) (-50.235085620125959778)

/-- logs/certify_cx3.log line 4: k = n + 1. -/
def certCR : Prop := GramCert 3315629582
  1155454823.437759756800183414149812 1155454823.437759756800183414150188
  0.2817344137619704801 0.2817344137619705199

theorem counterexampleC (hL : certCL) (hM : certCM) (hR : certCR) :
    ∃ gl g gr : ℝ, IsGram (3315629581 - 1) gl ∧ IsGram 3315629581 g ∧
      IsGram (3315629581 + 1) gr ∧ Good (3315629581 - 1) gl ∧
      (-1 : ℝ) ^ (3315629581 : ℤ) * hardyZ g < 0 ∧ Good (3315629581 + 1) gr ∧
      |deriv hardyZ g| < 4 * |hardyZ g| := by
  obtain ⟨gl, ⟨hl1, -⟩, hlθ, hlz, -⟩ := hL
  obtain ⟨g, ⟨hg1, -⟩, hgθ, ⟨hgz, -⟩, hgd1, hgd2⟩ := hM
  obtain ⟨gr, ⟨hr1, -⟩, hrθ, hrz, -⟩ := hR
  have hodd : Odd (3315629581 : ℤ) := ⟨1657814790, by norm_num⟩
  have hel : Even ((3315629581 : ℤ) - 1) := ⟨1657814790, by norm_num⟩
  have her : Even ((3315629581 : ℤ) + 1) := ⟨1657814791, by norm_num⟩
  have hz : 0 < hardyZ g := lt_of_lt_of_le (by norm_num) hgz
  have hd : hardyZ' g < 0 := lt_of_le_of_lt hgd2 (by norm_num)
  refine ⟨gl, g, gr, ⟨le_trans (by norm_num) hl1, by rw [hlθ]; norm_num⟩,
    ⟨le_trans (by norm_num) hg1, hgθ⟩, ⟨le_trans (by norm_num) hr1, by rw [hrθ]; norm_num⟩,
    good_of_even_pos hel (lt_of_lt_of_le (by norm_num) hlz), gramLaw_neg_of_odd_pos hodd hz,
    good_of_even_pos her (lt_of_lt_of_le (by norm_num) hrz), ?_⟩
  rw [deriv_hardyZ, abs_of_neg hd, abs_of_pos hz]
  have key : (50.235085620125960222 : ℝ) < 4 * 12.80566445995212596 := by norm_num
  linarith

theorem not_jerbyConj_C (hθ : ThetaDerivPos) (hL : certCL) (hM : certCM) (hR : certCR) :
    ¬ JerbyConj := fun H => by
  obtain ⟨gl, g, gr, h1, h2, h3, h4, h5, h6, h7⟩ := counterexampleC hL hM hR
  exact not_jerbyConjEx_of h1 h2 h3 h4 (bad_of_gramLaw_neg h5) h6 h7.le (jerbyConj_imp_ex hθ H)

end JerbyLean
