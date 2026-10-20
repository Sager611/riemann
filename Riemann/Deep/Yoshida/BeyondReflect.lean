/-
Riemann/Deep/Yoshida/BeyondReflect.lean — the `ℤ`-ball checker for `W3` certificates and the certificate theorem.

At a grid point `r = i/256`: `θ = rℓ/8` (`ℓ = log 2`) and `θ_J = 2ar` are reduced mod `2π`; `cos/sin` by
`BeyondBall`; harmonics `cos(jθ)` by angle addition. The prime is harmonic 8 (`−√2 log 2 · cos(8θ)`), the free
terms are harmonics `j ≥ 9` (`x = jℓ/8 > 2a`) with coefficients `E/2³⁰`. Cells are checked in `ℚ` from these
lower bounds (`affOK3_of`); the moment matrices reuse `Reflect`.
-/
import Riemann.Deep.Yoshida.BeyondBall

noncomputable section

open Finset

namespace Fayed.Yoshida.Beyond.Reflect5

open Fayed.Yoshida Fayed.Yoshida.Beyond Fayed.Yoshida.Beyond.Ball Fayed.Yoshida.Reflect
  Fayed.Yoshida.AReflect Zeta23

/-! ### Harmonics -/

def rot (c1 s1 : Bl) (p : Bl × Bl) : Bl × Bl :=
  (add (mul p.1 c1) (neg (mul p.2 s1)), add (mul p.2 c1) (mul p.1 s1))

lemma rot_sound {θ φ : ℝ} {c1 s1 : Bl} {p : Bl × Bl} (hc : Mem (Real.cos θ) c1)
    (hs : Mem (Real.sin θ) s1) (h1 : Mem (Real.cos φ) p.1) (h2 : Mem (Real.sin φ) p.2) :
    Mem (Real.cos (φ + θ)) (rot c1 s1 p).1 ∧ Mem (Real.sin (φ + θ)) (rot c1 s1 p).2 := by
  rw [Real.cos_add, Real.sin_add, sub_eq_add_neg]
  exact ⟨mem_add (mem_mul h1 hc) (mem_neg (mem_mul h2 hs)), mem_add (mem_mul h2 hc) (mem_mul h1 hs)⟩

/-- Accumulate `Σ_k b_k cos((j+1+k)θ)` from the state `p = (cos jθ, sin jθ)`. -/
def tacc (c1 s1 : Bl) : List Bl → Bl × Bl → Bl → Bl
  | [], _, acc => acc
  | b :: bs, p, acc => tacc c1 s1 bs (rot c1 s1 p) (add acc (mul b (rot c1 s1 p).1))

def tsum (θ : ℝ) : List ℝ → ℕ → ℝ
  | [], _ => 0
  | c :: cs, j => c * Real.cos (j * θ) + tsum θ cs (j + 1)

def Rel : List ℝ → List Bl → Prop
  | [], [] => True
  | c :: cs, b :: bs => Mem c b ∧ Rel cs bs
  | _, _ => False

lemma tacc_sound {θ : ℝ} {c1 s1 : Bl} (hc : Mem (Real.cos θ) c1) (hs : Mem (Real.sin θ) s1) :
    ∀ (cs : List ℝ) (bs : List Bl) (j : ℕ) (p : Bl × Bl) (A : ℝ) (acc : Bl), Rel cs bs →
    Mem (Real.cos (j * θ)) p.1 → Mem (Real.sin (j * θ)) p.2 → Mem A acc →
    Mem (A + tsum θ cs (j + 1)) (tacc c1 s1 bs p acc)
  | [], [], _, _, A, acc, _, _, _, hA => by simpa [tacc, tsum] using hA
  | c :: cs, b :: bs, j, p, A, acc, hR, h1, h2, hA => by
    obtain ⟨hb, hR'⟩ := hR
    obtain ⟨k1, k2⟩ := rot_sound hc hs h1 h2
    rw [show (j : ℝ) * θ + θ = ((j + 1 : ℕ) : ℝ) * θ by push_cast; ring] at k1 k2
    have hA' := mem_add hA (mem_mul hb k1)
    have := tacc_sound hc hs cs bs (j + 1) (rot c1 s1 p) _ _ hR' k1 k2 hA'
    simp only [tacc, tsum]
    rw [show A + (c * Real.cos (((j + 1 : ℕ) : ℝ) * θ) + tsum θ cs (j + 1 + 1))
        = A + c * Real.cos (((j + 1 : ℕ) : ℝ) * θ) + tsum θ cs (j + 1 + 1) by ring]
    exact this
  | [], _ :: _, _, _, _, _, h, _, _, _ => absurd h (by simp [Rel])
  | _ :: _, [], _, _, _, _, h, _, _, _ => absurd h (by simp [Rel])

/-! ### Free terms -/

/-- Free terms `(E/2³⁰, jℓ/8)` for `j = j₀, j₀+1, …`. -/
def fsAux (j : ℕ) : List ℤ → List (ℝ × ℝ)
  | [] => []
  | e :: es => ((e : ℝ) / 2 ^ 30, (j : ℝ) * Real.log 2 / 8) :: fsAux (j + 1) es

lemma trigF_fsAux (r : ℝ) : ∀ (E : List ℤ) (j : ℕ),
    trigF (fsAux j E) r = tsum (r * Real.log 2 / 8) (E.map fun e => (e : ℝ) / 2 ^ 30) j
  | [], _ => rfl
  | e :: es, j => by
    simp only [fsAux, trigF, List.map_cons, tsum, trigF_fsAux r es (j + 1)]
    congr 2
    ring

lemma fsAux_freq : ∀ (E : List ℤ) (j : ℕ), ∀ p ∈ fsAux j E, (j : ℝ) * Real.log 2 / 8 ≤ p.2
  | [], _, p, hp => by simp [fsAux] at hp
  | e :: es, j, p, hp => by
    simp only [fsAux, List.mem_cons] at hp
    rcases hp with rfl | hp
    · exact le_refl _
    · have := fsAux_freq es (j + 1) p hp
      have hl := Real.log_pos (by norm_num : (1 : ℝ) < 2)
      push_cast at this
      nlinarith

def absQ : List ℤ → ℚ
  | [] => 0
  | e :: es => |(e : ℚ)| / 2 ^ 30 + absQ es

lemma absSum0_fsAux : ∀ (E : List ℤ) (j : ℕ), absSum0 (fsAux j E) = (absQ E : ℝ)
  | [], _ => by simp [fsAux, absSum0, absQ]
  | e :: es, j => by
    simp only [fsAux, absSum0, absQ, absSum0_fsAux es (j + 1)]
    push_cast
    rw [abs_div, abs_of_pos (by norm_num : (0 : ℝ) < 2 ^ 30)]

/-- `Σ |E|/2³⁰ (jℓ/8)²` with `ℓ ≤ 0.6931471808`. -/
def curvQ (j : ℕ) : List ℤ → ℚ
  | [] => 0
  | e :: es => |(e : ℚ)| / 2 ^ 30 * ((j : ℚ) * (6931471808 / 10 ^ 10) / 8) ^ 2 + curvQ (j + 1) es

lemma absSum_fsAux_le : ∀ (E : List ℤ) (j : ℕ), absSum (fsAux j E) ≤ (curvQ j E : ℝ)
  | [], _ => by simp [fsAux, absSum, curvQ]
  | e :: es, j => by
    simp only [fsAux, absSum, curvQ]
    have ih := absSum_fsAux_le es (j + 1)
    have hl := Real.log_two_lt_d9
    have hl0 := Real.log_pos (by norm_num : (1 : ℝ) < 2)
    push_cast at ih ⊢
    have h1 : ((j : ℝ) * Real.log 2 / 8) ^ 2 ≤ ((j : ℝ) * (6931471808 / 10 ^ 10) / 8) ^ 2 := by
      have hj : (0 : ℝ) ≤ j := Nat.cast_nonneg j
      apply pow_le_pow_left₀ (by positivity)
      have : Real.log 2 ≤ 6931471808 / 10 ^ 10 := by norm_num at hl ⊢; linarith
      gcongr
    rw [abs_div, abs_of_pos (by norm_num : (0 : ℝ) < 2 ^ 30)]
    have h2 : 0 ≤ |(e : ℝ)| / 2 ^ 30 := by positivity
    nlinarith [mul_le_mul_of_nonneg_left h1 h2]

lemma cP_le : cP ≤ 98027 / 100000 := by
  unfold cP
  have hl := Real.log_two_lt_d9
  have hs : Real.sqrt 2 ≤ 141421357 / 100000000 := by
    rw [Real.sqrt_le_left (by norm_num)]
    norm_num
  have hl0 := Real.log_pos (by norm_num : (1 : ℝ) < 2)
  norm_num at hl
  nlinarith [Real.sqrt_nonneg 2]

lemma curvP_fsAux_le (E : List ℤ) :
    curvP (fsAux 9 E) ≤ 98027 / 100000 * (6931471808 / 10 ^ 10) ^ 2 + (curvQ 9 E : ℝ) := by
  unfold curvP
  have h1 := absSum_fsAux_le E 9
  have hl := Real.log_two_lt_d9
  have hl0 := Real.log_pos (by norm_num : (1 : ℝ) < 2)
  have h2 : Real.log 2 ^ 2 ≤ (6931471808 / 10 ^ 10 : ℝ) ^ 2 := by
    apply pow_le_pow_left₀ hl0.le
    norm_num at hl ⊢
    linarith
  have h3 := cP_le
  have h4 := cP_nonneg
  nlinarith [mul_le_mul h3 h2 (by positivity) (by norm_num)]

/-! ### Coefficients: harmonics 1..7 zero, 8 the prime, 9.. free -/

def csR (E : List ℤ) : List ℝ :=
  0 :: 0 :: 0 :: 0 :: 0 :: 0 :: 0 :: (-cP) :: E.map fun e => (e : ℝ) / 2 ^ 30

def bsB (E : List ℤ) : List Bl :=
  (0, 0) :: (0, 0) :: (0, 0) :: (0, 0) :: (0, 0) :: (0, 0) :: (0, 0) :: neg cPB ::
    E.map fun e => (e * 4294967296, 0)

lemma mem_zero : Mem 0 (0, 0) := mem_of_eq (by simp)

lemma Mem.congr' {x y : ℝ} {b : Bl} (h : x = y) (hx : Mem x b) : Mem y b := h ▸ hx

lemma rel_map : ∀ E : List ℤ, Rel (E.map fun e => (e : ℝ) / 2 ^ 30) (E.map fun e => (e * 4294967296, 0))
  | [] => trivial
  | e :: es => ⟨mem_of_eq (by unfold Ball.S; push_cast; ring), rel_map es⟩

lemma rel_cs (E : List ℤ) : Rel (csR E) (bsB E) :=
  ⟨mem_zero, mem_zero, mem_zero, mem_zero, mem_zero, mem_zero, mem_zero, mem_neg mem_cP, rel_map E⟩

lemma tsum_cs (θ : ℝ) (E : List ℤ) :
    tsum θ (csR E) 1 = -cP * Real.cos (8 * θ) + tsum θ (E.map fun e => (e : ℝ) / 2 ^ 30) 9 := by
  simp only [csR, tsum]
  push_cast
  ring

/-! ### Angles -/

/-- Reduce an angle ball by `2kπ`. -/
def angB (t : Bl) : Bl := add t (smul (-2 * ((t.1 + piB.1) / (2 * piB.1))) piB)

lemma mem_angB {θ : ℝ} {t : Bl} (h : Mem θ t) :
    Mem (Real.cos θ) (cosB (angB t)) ∧ Mem (Real.sin θ) (sinB (angB t)) := by
  set k : ℤ := (t.1 + piB.1) / (2 * piB.1)
  have hy := mem_add h (mem_smul (-2 * k) mem_pi)
  have e : θ + ((-2 * k : ℤ) : ℝ) * Real.pi = θ - (k : ℝ) * (2 * Real.pi) := by push_cast; ring
  rw [e] at hy
  have hc := mem_cosB hy
  have hs := mem_sinB hy
  rw [Real.cos_sub_int_mul_two_pi] at hc
  rw [Real.sin_sub_int_mul_two_pi] at hs
  exact ⟨hc, hs⟩

/-- `θ = i ℓ / 2¹¹` (`= rℓ/8`, `r = i/256`). -/
def t1 (i : ℕ) : Bl := divN 2048 (smul i lB)

lemma mem_t1 (i : ℕ) : Mem ((i : ℝ) / 256 * Real.log 2 / 8) (t1 i) := by
  have h := mem_divN (by norm_num : (0 : ℤ) < 2048) (mem_smul (i : ℤ) mem_log2)
  have e : ((i : ℤ) : ℝ) * Real.log 2 / ((2048 : ℤ) : ℝ) = (i : ℝ) / 256 * Real.log 2 / 8 := by
    push_cast; ring
  rw [e] at h
  exact h

/-- `θ_J = 2 a r = num·i/(den·128)` for `a = num/den`. -/
def t2 (a : ℚ) (i : ℕ) : Bl := divN ((a.den : ℤ) * 128) (smul (a.num * i) (Ball.S, 0))

lemma mem_t2 (a : ℚ) (i : ℕ) : Mem (2 * (a : ℝ) * ((i : ℝ) / 256)) (t2 a i) := by
  have hd : (0 : ℤ) < (a.den : ℤ) * 128 := by have := a.den_pos; positivity
  have h := mem_divN hd (mem_smul (a.num * i) Ball.mem_one)
  have e : ((a.num * i : ℤ) : ℝ) * 1 / (((a.den : ℤ) * 128 : ℤ) : ℝ) = 2 * (a : ℝ) * ((i : ℝ) / 256) := by
    have ha : (a : ℝ) = (a.num : ℝ) / (a.den : ℝ) := by
      rw [← Rat.cast_intCast, ← Rat.cast_natCast, ← Rat.cast_div, Rat.num_div_den]
    have hd' : (0 : ℝ) < (a.den : ℝ) := by exact_mod_cast a.den_pos
    rw [ha]
    push_cast
    field_simp
    ring
  rw [e] at h
  exact h

/-- `1/(1/4 + r²)` at `r = i/256`. -/
def invB (i : ℕ) : Bl := (Ball.S * 262144 / (65536 + 4 * (i : ℤ) ^ 2), 1)

lemma mem_invB (i : ℕ) : Mem (1 / (1 / 4 + ((i : ℝ) / 256) ^ 2)) (invB i) := by
  have hq : (0 : ℤ) < 65536 + 4 * (i : ℤ) ^ 2 := by positivity
  have k1 := ediv_le (Ball.S * 262144) hq
  have k2 := lt_ediv (Ball.S * 262144) hq
  unfold Mem invB
  push_cast at k1 k2 ⊢
  have e : 1 / (1 / 4 + ((i : ℝ) / 256) ^ 2) * Ball.S = Ball.S * 262144 / (65536 + 4 * (i : ℝ) ^ 2) := by
    field_simp
    ring
  rw [e, abs_le]
  constructor <;> linarith

/-! ### The weight at a grid point -/

/-- Ball for `J(a, r) + trigP(fs, r)` at `r = i/256`. -/
def jtB (a : ℚ) (E : List ℤ) (shB chB : Bl) (i : ℕ) : Bl :=
  add (tacc (cosB (angB (t1 i))) (sinB (angB (t1 i))) (bsB E) ((Ball.S, 0), (0, 0)) (0, 0))
    (mul (add (smul 2 (mul shB (cosB (angB (t2 a i)))))
      (smul 4 (divN 256 (smul i (mul chB (sinB (angB (t2 a i)))))))) (invB i))

theorem jt_sound (a : ℚ) (E : List ℤ) {shB chB : Bl} (hsh : Mem (Real.sinh a) shB)
    (hch : Mem (Real.cosh a) chB) (i : ℕ) :
    Mem (J a ((i : ℝ) / 256) + trigP (fsAux 9 E) ((i : ℝ) / 256)) (jtB a E shB chB i) := by
  set r : ℝ := (i : ℝ) / 256 with hr
  obtain ⟨hc1, hs1⟩ := mem_angB (mem_t1 i)
  obtain ⟨hc2, hs2⟩ := mem_angB (mem_t2 a i)
  have hT := tacc_sound hc1 hs1 (csR E) (bsB E) 0 ((Ball.S, 0), (0, 0)) 0 (0, 0) (rel_cs E)
    (by simpa using Ball.mem_one) (by simpa using mem_zero) mem_zero
  rw [zero_add, tsum_cs] at hT
  have hJ := mem_mul (mem_add (mem_smul 2 (mem_mul hsh hc2))
    (mem_smul 4 (mem_divN (by norm_num : (0 : ℤ) < 256) (mem_smul (i : ℤ) (mem_mul hch hs2)))))
    (mem_invB i)
  have hsum := mem_add hT hJ
  have c8 : Real.cos (8 * (r * Real.log 2 / 8)) = Real.cos (r * Real.log 2) := by congr 1; ring
  refine Mem.congr' ?_ hsum
  unfold J trigP
  rw [trigF_fsAux, c8, hr]
  push_cast
  field_simp
  ring

/-! ### Cells and the walk -/

def cellB (s Ej lo hi A B mL jL mH jH : ℚ) : Bool :=
  decide (0 ≤ lo) && decide (lo < hi) && decide (mL ≤ mH) && decide (0 ≤ A + B * lo) &&
    decide (0 ≤ A + B * hi) &&
    decide (s - (A + B * lo) ≤ mL + jL - ((mH - mL) * ((hi - lo) / (4 * (hi + lo))) + Ej * (hi - lo) ^ 2)) &&
    decide (s - (A + B * hi) ≤ mH + jH - ((mH - mL) * ((hi - lo) / (4 * (hi + lo))) + Ej * (hi - lo) ^ 2))

/-- Certificate configuration. -/
structure Cfg where
  a : ℚ
  s : ℚ
  Ej : ℚ
  E : List ℤ
  shB : Bl
  chB : Bl
  M : ℕ
  Mp : ℕ

structure CfgValid (C : Cfg) : Prop where
  ha : 0 < (C.a : ℝ)
  sh : Mem (Real.sinh C.a) C.shB
  ch : Mem (Real.cosh C.a) C.chB
  M : C.M ≤ C.Mp
  Ej : 4 / 3 * (C.a : ℝ) ^ 3 * Real.cosh C.a + curvP (fsAux 9 C.E) / 8 ≤ C.Ej

def muQ5 (C : Cfg) (i : ℕ) : ℚ := (muLowZ 8 i C.M C.Mp : ℚ) / (Ball.S : ℚ)

def jtQ (C : Cfg) (i : ℕ) : ℚ := (((jtB C.a C.E C.shB C.chB i).1 - (jtB C.a C.E C.shB C.chB i).2 : ℤ) : ℚ) / (Ball.S : ℚ)

lemma muQ5_le (C : Cfg) (hC : CfgValid C) (i : ℕ) :
    ((muQ5 C i : ℚ) : ℝ) ≤ 2 * Real.pi * mu ((i : ℝ) / 256) := by
  have h := muLowZ_le 8 i C.M C.Mp hC.M
  unfold muQ5
  push_cast
  rw [div_le_iff₀ S_posR]
  norm_num at h ⊢
  linarith

lemma jtQ_le (C : Cfg) (hC : CfgValid C) (i : ℕ) :
    ((jtQ C i : ℚ) : ℝ) ≤ J C.a ((i : ℝ) / 256) + trigP (fsAux 9 C.E) ((i : ℝ) / 256) := by
  have h := low_le (jt_sound C.a C.E hC.sh hC.ch i)
  unfold jtQ
  push_cast at h ⊢
  exact h

def chkB (C : Cfg) : ℕ → ℚ → ℚ → List (ℕ × ℕ × ℤ × ℤ) → ℕ → Bool
  | x, _, _, [], R => decide (x = R)
  | x, mL, jL, c :: cs, R =>
    decide (c.1 = x) && cellB C.s C.Ej ((c.1 : ℚ) / 2 ^ 8) ((c.2.1 : ℚ) / 2 ^ 8) ((c.2.2.1 : ℚ) / 2 ^ 30)
      ((c.2.2.2 : ℚ) / 2 ^ 30) mL jL (muQ5 C c.2.1) (jtQ C c.2.1)
      && chkB C c.2.1 (muQ5 C c.2.1) (jtQ C c.2.1) cs R

theorem chkB_sound (C : Cfg) (hC : CfgValid C) :
    ∀ (cs : List (ℕ × ℕ × ℤ × ℤ)) (x : ℕ) (mL jL : ℚ) (R : ℕ),
    ((mL : ℚ) : ℝ) ≤ 2 * Real.pi * mu ((x : ℝ) / 256) →
    ((jL : ℚ) : ℝ) ≤ J C.a ((x : ℝ) / 256) + trigP (fsAux 9 C.E) ((x : ℝ) / 256) →
    chkB C x mL jL cs R = true →
    Chain4 (((x : ℚ) / 2 ^ 8 : ℚ) : ℝ) (toR4 (cs.map (toQ4 8 30))) (((R : ℚ) / 2 ^ 8 : ℚ) : ℝ) ∧
      ∀ c ∈ toR4 (cs.map (toQ4 8 30)), AffOK3 C.a C.s (fsAux 9 C.E) c
  | [], x, mL, jL, R, _, _, h => by
    simp only [chkB, decide_eq_true_eq] at h
    simp [toR4, Chain4, h]
  | c :: cs, x, mL, jL, R, hmL, hjL, h => by
    simp only [chkB, Bool.and_eq_true, decide_eq_true_eq] at h
    obtain ⟨⟨hcx, hcell⟩, hrest⟩ := h
    have ih := chkB_sound C hC cs c.2.1 _ _ R (muQ5_le C hC c.2.1) (jtQ_le C hC c.2.1) hrest
    refine ⟨?_, ?_⟩
    · simp only [toR4, List.map_cons, Chain4, toQ4]
      exact ⟨by rw [hcx], ih.1⟩
    · intro d hd
      simp only [toR4, List.map_cons, List.mem_cons] at hd
      rcases hd with rfl | hd
      · unfold cellB at hcell
        simp only [Bool.and_eq_true, decide_eq_true_eq] at hcell
        obtain ⟨⟨⟨⟨⟨⟨h0, h1⟩, h2⟩, h3⟩, h4⟩, h5⟩, h6⟩ := hcell
        have cv (u v : ℚ) (huv : u ≤ v) : (u : ℝ) ≤ (v : ℝ) := by exact_mod_cast huv
        have h5' := cv _ _ h5
        have h6' := cv _ _ h6
        have h3' := cv _ _ h3
        have h4' := cv _ _ h4
        push_cast at h5' h6' h3' h4'
        rw [← hcx] at hmL hjL
        have elo : ((c.1 : ℝ) / 256) = ((c.1 : ℚ) / 2 ^ 8 : ℚ) := by push_cast; norm_num
        have ehi : ((c.2.1 : ℝ) / 256) = ((c.2.1 : ℚ) / 2 ^ 8 : ℚ) := by push_cast; norm_num
        have hmH := muQ5_le C hC c.2.1
        have hjH := jtQ_le C hC c.2.1
        rw [elo] at hmL hjL
        rw [ehi] at hmH hjH
        have h1r := cv _ _ h1.le
        have h1s : (((c.1 : ℚ) / 2 ^ 8 : ℚ) : ℝ) < (((c.2.1 : ℚ) / 2 ^ 8 : ℚ) : ℝ) := by exact_mod_cast h1
        simp only [toQ4]
        push_cast at hmL hjL hmH hjH h1s ⊢
        exact affOK3_of hC.ha.le (by positivity) h1s hmL hmH hjL hjH hC.Ej
          (by exact_mod_cast h2) h3' h4' h5' h6'
      · exact ih.2 d hd

/-! ### The tail -/

def tailB (C : Cfg) (Mt : ℕ) (R : ℕ) (shhi chhi : ℚ) : Bool :=
  decide (256 ≤ R) && decide (C.s ≤ (muLowZ 8 R Mt C.Mp : ℚ) / (Ball.S : ℚ)
    - 4 * (shhi / (2 * ((R : ℚ) / 256) ^ 2) + chhi / ((R : ℚ) / 256)) - 98027 / 100000 - absQ C.E)

theorem tailB_sound (C : Cfg) (hC : CfgValid C) {Mt R : ℕ} {shhi chhi : ℚ} (hMt : Mt ≤ C.Mp)
    (hsh : Real.sinh C.a ≤ shhi) (hch : Real.cosh C.a ≤ chhi) (h : tailB C Mt R shhi chhi = true) :
    ∀ r : ℝ, (((R : ℚ) / 2 ^ 8 : ℚ) : ℝ) ≤ |r| → (C.s : ℝ) ≤ W3 C.a (fsAux 9 C.E) r := by
  intro r hr
  unfold tailB at h
  simp only [Bool.and_eq_true, decide_eq_true_eq] at h
  obtain ⟨hR, hs⟩ := h
  have hR1 : (1 : ℝ) ≤ (R : ℝ) / 256 := by
    have : (256 : ℝ) ≤ R := by exact_mod_cast hR
    linarith
  have hr' : (R : ℝ) / 256 ≤ |r| := by push_cast at hr; norm_num at hr ⊢; exact hr
  have hmu := muLowZ_le 8 R Mt C.Mp hMt
  have hmono := mu_abs_mono (by positivity) hr'
  have hJ := J_ge_tail hC.ha.le hR1 hr'
  have hcos : -cP ≤ -cP * Real.cos (r * Real.log 2) := by
    have := Real.cos_le_one (r * Real.log 2)
    have := Real.neg_one_le_cos (r * Real.log 2)
    nlinarith [cP_nonneg]
  have htF := (abs_le.1 (abs_trigF_le (fsAux 9 C.E) r)).1
  rw [absSum0_fsAux] at htF
  have hcP := cP_le
  have hs' : ((C.s : ℚ) : ℝ) ≤ (((muLowZ 8 R Mt C.Mp : ℚ) / (Ball.S : ℚ)
      - 4 * (shhi / (2 * ((R : ℚ) / 256) ^ 2) + chhi / ((R : ℚ) / 256)) - 98027 / 100000 - absQ C.E : ℚ) : ℝ) := by
    exact_mod_cast hs
  push_cast at hs'
  have hSp := S_posR
  have hm2 : ((muLowZ 8 R Mt C.Mp : ℤ) : ℝ) / Ball.S ≤ 2 * Real.pi * mu ((R : ℝ) / 2 ^ 8) := by
    rw [div_le_iff₀ hSp]; linarith
  have hR0 : (0 : ℝ) < (R : ℝ) / 256 := by linarith
  have k : Real.sinh C.a / (2 * ((R : ℝ) / 256) ^ 2) + Real.cosh C.a / ((R : ℝ) / 256)
      ≤ shhi / (2 * ((R : ℝ) / 256) ^ 2) + chhi / ((R : ℝ) / 256) := by gcongr
  unfold W3 W trigP
  norm_num at hm2 hs' ⊢
  have hπ := Real.pi_pos
  nlinarith [mul_le_mul_of_nonneg_left hmono (by positivity : (0 : ℝ) ≤ 2 * Real.pi)]

/-! ### The certificate theorem -/

/-- **Rational/`ℤ`-ball certificate for `W3` ⇒ Weil positivity on every window `L ≤ 2a`.** -/
theorem weilPos_of_bcert (C : Cfg) (hC : CfgValid C) (h3 : 2 * (C.a : ℝ) < Real.log 3)
    (hfs : 2 * (C.a : ℝ) < 9 * Real.log 2 / 8)
    (cN : List (ℕ × ℕ × ℤ × ℤ)) (R : ℕ) (K : ℕ) (momD : List ℚ)
    (Lm Gm LGm Mm U MU Dm Wm : List (List ℚ)) (α ε β : ℚ) (Mt : ℕ) (shhi chhi : ℚ)
    (hMt : Mt ≤ C.Mp) (hsh : Real.sinh C.a ≤ shhi) (hch : Real.cosh C.a ≤ chhi)
    (hcells : chkB C 0 (muQ5 C 0) (jtQ C 0) cN R = true)
    (htail : tailB C Mt R shhi chhi = true)
    (hK : kQ4 C.a K (cN.map (toQ4 8 30)) = true) (hmom : momCheckA cN 8 30 K momD = true)
    (hGm : glQ K C.a Lm Gm = true) (hLGm : tmulQ K Lm Gm LGm = true)
    (hM : mQ K α Lm LGm Mm momD = true) (hMU : mulQ K Mm U MU = true)
    (hD : tmulQ K U MU Dm = true) (hUW : invQ K U Wm = true) (hdom : domQ K Dm = true)
    (hnum : decide (0 < ε ∧ 0 ≤ α ∧ 0 ≤ C.s
      ∧ 8 * C.a ^ (2 * K + 1) * momD.getD (2 * K) 0 / ((K.factorial : ℚ)) ^ 2 ≤ β
      ∧ (1 + ε) * α + (1 + 1 / ε) * β ≤ 2 * piLo * C.s) = true)
    {L : ℝ} (hL1 : L ≤ 2 * (C.a : ℝ)) : Riemann.Deep.Weil.WeilPos L := by
  simp only [decide_eq_true_eq] at hnum
  obtain ⟨hε, hα, hs, hβ, hn⟩ := hnum
  set cs := cN.map (toQ4 8 30) with hcs
  have hmom' : ∀ p ≤ 2 * K, ((momD.getD p 0 : ℚ) : ℝ) = affSum (toR4 cs) p := by
    intro p hp
    have := List.all_eq_true.1 hmom p (List.mem_range.2 (by omega))
    simp only [decide_eq_true_eq] at this
    rw [this, ← momA_grid 8 30 p cN, ← hcs, cast_momA]
  have hGm' := all_range₂ hGm
  have hM' := all_range₂ hM
  have hpsd := hpsd_of K (α : ℝ) (entR Lm) (entR Gm) (entR LGm) (entR Mm) (Ghat (C.a : ℝ))
    (fun p => affSum (toR4 cs) p)
    (fun k hk j hj => by
      unfold entR
      rw [hGm' k hk j hj]
      push_cast
      simp_rw [cast_ghatQ])
    (tmulQ_sound hLGm)
    (fun i hi j hj => by
      unfold entR
      rw [hM' i hi j hj]
      push_cast
      rw [hmom' (i + j) (by omega)])
    (qf_nonneg_of_cert K (entR Mm) (entR U) (entR Wm) (entR MU) (entR Dm) (mulQ_sound hMU)
      (tmulQ_sound hD) (invQ_sound hUW) (domQ_sound hdom))
  have hβ' : 8 * (C.a : ℝ) ^ (2 * K + 1) * affSum (toR4 cs) (2 * K) / (K.factorial : ℝ) ^ 2 ≤ β := by
    rw [← hmom' (2 * K) le_rfl]
    have : ((8 * C.a ^ (2 * K + 1) * momD.getD (2 * K) 0 / ((K.factorial : ℚ)) ^ 2 : ℚ) : ℝ) ≤ (β : ℝ) := by
      exact_mod_cast hβ
    push_cast at this
    exact this
  have hn' : (1 + (ε : ℝ)) * α + (1 + 1 / (ε : ℝ)) * β ≤ 2 * Real.pi * C.s := by
    have h1 : (((1 + ε) * α + (1 + 1 / ε) * β : ℚ) : ℝ) ≤ ((2 * piLo * C.s : ℚ) : ℝ) := by
      exact_mod_cast hn
    push_cast at h1
    have h2 := piLo_lt
    have h3' : (0 : ℝ) ≤ C.s := by exact_mod_cast hs
    nlinarith
  have h0 := chkB_sound C hC cN 0 _ _ R (muQ5_le C hC 0) (jtQ_le C hC 0) hcells
  have e0 : ((((0 : ℕ) : ℚ) / 2 ^ 8 : ℚ) : ℝ) = 0 := by simp
  rw [e0] at h0
  have hfs' : ∀ p ∈ fsAux 9 C.E, 2 * (C.a : ℝ) < |p.2| := by
    intro p hp
    have := fsAux_freq C.E 9 p hp
    push_cast at this
    exact lt_of_lt_of_le hfs (this.trans (le_abs_self _))
  exact weilPos_of_affine3 hC.ha le_rfl h3 (fsAux 9 C.E) hfs' (toR4 cs) h0.1 h0.2
    (tailB_sound C hC hMt hsh hch htail) K (kQ4_sound hK) (entR Lm) hpsd (by exact_mod_cast hε)
    (by exact_mod_cast hα) hβ' hn' hL1

end Fayed.Yoshida.Beyond.Reflect5
