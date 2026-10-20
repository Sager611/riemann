/-
Riemann/Deep/Yoshida/BeyondBall.lean — fixed-point ball arithmetic over `ℤ` for fast kernel checks.

A ball `(m, ρ) : ℤ × ℤ` holds the reals `x` with `|x·2⁶² − m| ≤ ρ`. Sums, products (rounded down, radius
rounded up), integer scaling and division by positive integers are sound. On top: `cos`/`sin` by range
reduction and Taylor loops (remainder from `Complex.exp_bound'`), harmonics by angle addition, and a lower
bound for `2πμ` (digamma partial sum + telescoping tail, `two_pi_mu_ge_tel`). The kernel evaluates `ℤ`
arithmetic about 20× faster than `ℚ` arithmetic (`.scratch/yoshida5/Bench3.lean`).
-/
import Riemann.Deep.Yoshida.BeyondCore
import Riemann.Deep.Yoshida.AReflect

noncomputable section

open Finset

namespace Fayed.Yoshida.Beyond.Ball

open Fayed.Yoshida Fayed.Yoshida.Taylor Zeta23

/-- Fixed-point scale `2⁶²`. -/
def S : ℤ := 4611686018427387904

lemma S_pos : (0 : ℤ) < S := by norm_num [S]

lemma S_posR : (0 : ℝ) < (S : ℝ) := by exact_mod_cast S_pos

abbrev Bl := ℤ × ℤ

/-- `x ∈ (m, ρ)` iff `|x·S − m| ≤ ρ`. -/
def Mem (x : ℝ) (b : Bl) : Prop := |x * S - b.1| ≤ b.2

lemma Mem.rad {x : ℝ} {b : Bl} (h : Mem x b) : (0 : ℝ) ≤ b.2 := (abs_nonneg _).trans h

def add (b c : Bl) : Bl := (b.1 + c.1, b.2 + c.2)
def neg (b : Bl) : Bl := (-b.1, b.2)
def smul (k : ℤ) (b : Bl) : Bl := (k * b.1, |k| * b.2)
def mul (b c : Bl) : Bl := (b.1 * c.1 / S, (|b.1| * c.2 + |c.1| * b.2 + b.2 * c.2) / S + 2)
def divN (q : ℤ) (b : Bl) : Bl := (b.1 / q, b.2 / q + 2)

lemma ediv_le (a : ℤ) {q : ℤ} (hq : 0 < q) : ((a / q : ℤ) : ℝ) ≤ (a : ℝ) / q := by
  have hq' : (0 : ℝ) < q := by exact_mod_cast hq
  rw [le_div_iff₀ hq']
  exact_mod_cast Int.ediv_mul_le a hq.ne'

lemma lt_ediv (a : ℤ) {q : ℤ} (hq : 0 < q) : (a : ℝ) / q < ((a / q : ℤ) : ℝ) + 1 := by
  have hq' : (0 : ℝ) < q := by exact_mod_cast hq
  rw [div_lt_iff₀ hq']
  exact_mod_cast Int.lt_ediv_add_one_mul_self a hq

lemma mem_add {x y : ℝ} {b c : Bl} (hx : Mem x b) (hy : Mem y c) : Mem (x + y) (add b c) := by
  unfold Mem add at *
  push_cast
  calc |(x + y) * S - (b.1 + c.1)| = |(x * S - b.1) + (y * S - c.1)| := by congr 1; ring
    _ ≤ |x * S - b.1| + |y * S - c.1| := abs_add_le _ _
    _ ≤ b.2 + c.2 := add_le_add hx hy

lemma mem_neg {x : ℝ} {b : Bl} (hx : Mem x b) : Mem (-x) (neg b) := by
  unfold Mem neg at *
  push_cast
  rw [show -x * S - -(b.1 : ℝ) = -(x * S - b.1) by ring, abs_neg]
  exact hx

lemma mem_smul {x : ℝ} {b : Bl} (k : ℤ) (hx : Mem x b) : Mem (k * x) (smul k b) := by
  unfold Mem smul at *
  push_cast
  rw [show (k : ℝ) * x * S - k * b.1 = k * (x * S - b.1) by ring, abs_mul]
  exact mul_le_mul_of_nonneg_left hx (abs_nonneg _)

lemma mem_mul {x y : ℝ} {b c : Bl} (hx : Mem x b) (hy : Mem y c) : Mem (x * y) (mul b c) := by
  have h0 := hx.rad
  have h0' := hy.rad
  unfold Mem mul at *
  have hS := S_posR
  have k1 := ediv_le (b.1 * c.1) S_pos
  have k2 := lt_ediv (b.1 * c.1) S_pos
  have k3 := lt_ediv (|b.1| * c.2 + |c.1| * b.2 + b.2 * c.2) S_pos
  push_cast at k1 k2 k3 ⊢
  set e1 := x * S - b.1 with he1
  set e2 := y * S - c.1 with he2
  set q1 : ℝ := ((b.1 * c.1 / S : ℤ) : ℝ)
  set q3 : ℝ := (((|b.1| * c.2 + |c.1| * b.2 + b.2 * c.2) / S : ℤ) : ℝ)
  have hxy : x * y * S = ((b.1 : ℝ) * c.1 + b.1 * e2 + c.1 * e1 + e1 * e2) / S := by
    rw [he1, he2]; field_simp; ring
  have hD : |(b.1 : ℝ) * e2 + c.1 * e1 + e1 * e2| ≤ |(b.1 : ℝ)| * c.2 + |(c.1 : ℝ)| * b.2 + b.2 * c.2 := by
    calc |(b.1 : ℝ) * e2 + c.1 * e1 + e1 * e2| ≤ |(b.1 : ℝ) * e2| + |(c.1 : ℝ) * e1| + |e1 * e2| :=
          abs_add_three _ _ _
      _ = |(b.1 : ℝ)| * |e2| + |(c.1 : ℝ)| * |e1| + |e1| * |e2| := by rw [abs_mul, abs_mul, abs_mul]
      _ ≤ |(b.1 : ℝ)| * c.2 + |(c.1 : ℝ)| * b.2 + b.2 * c.2 := by
          have := abs_nonneg e1
          have := abs_nonneg e2
          gcongr
  rw [hxy]
  have hsplit : ((b.1 : ℝ) * c.1 + b.1 * e2 + c.1 * e1 + e1 * e2) / S - q1
      = ((b.1 : ℝ) * c.1 / S - q1) + ((b.1 : ℝ) * e2 + c.1 * e1 + e1 * e2) / S := by ring
  rw [hsplit]
  have ha1 : |(b.1 : ℝ) * c.1 / S - q1| ≤ 1 := abs_le.2 ⟨by linarith, by linarith⟩
  have ha2 : |((b.1 : ℝ) * e2 + c.1 * e1 + e1 * e2) / S|
      ≤ (|(b.1 : ℝ)| * c.2 + |(c.1 : ℝ)| * b.2 + b.2 * c.2) / S := by
    rw [abs_div, abs_of_pos hS]
    exact div_le_div_of_nonneg_right hD hS.le
  calc _ ≤ |(b.1 : ℝ) * c.1 / S - q1| + |((b.1 : ℝ) * e2 + c.1 * e1 + e1 * e2) / S| := abs_add_le _ _
    _ ≤ 1 + (|(b.1 : ℝ)| * c.2 + |(c.1 : ℝ)| * b.2 + b.2 * c.2) / S := add_le_add ha1 ha2
    _ ≤ q3 + 2 := by linarith

lemma mem_divN {x : ℝ} {b : Bl} {q : ℤ} (hq : 0 < q) (hx : Mem x b) : Mem (x / q) (divN q b) := by
  have h0 := hx.rad
  unfold Mem divN at *
  have hq' : (0 : ℝ) < q := by exact_mod_cast hq
  have k1 := ediv_le b.1 hq
  have k2 := lt_ediv b.1 hq
  have k3 := lt_ediv b.2 hq
  push_cast at k1 k2 k3 ⊢
  set q1 : ℝ := ((b.1 / q : ℤ) : ℝ)
  set q3 : ℝ := ((b.2 / q : ℤ) : ℝ)
  have e : x / q * S - q1 = ((b.1 : ℝ) / q - q1) + (x * S - b.1) / q := by field_simp; ring
  rw [e]
  have ha1 : |(b.1 : ℝ) / q - q1| ≤ 1 := abs_le.2 ⟨by linarith, by linarith⟩
  have ha2 : |(x * S - b.1) / q| ≤ (b.2 : ℝ) / q := by
    rw [abs_div, abs_of_pos hq']
    exact div_le_div_of_nonneg_right hx hq'.le
  calc _ ≤ |(b.1 : ℝ) / q - q1| + |(x * S - b.1) / q| := abs_add_le _ _
    _ ≤ 1 + (b.2 : ℝ) / q := add_le_add ha1 ha2
    _ ≤ q3 + 2 := by linarith

lemma mem_of_eq {x : ℝ} {m : ℤ} (h : x * S = m) : Mem x (m, 0) := by
  unfold Mem
  simp [h]

lemma low_le {x : ℝ} {b : Bl} (h : Mem x b) : ((b.1 - b.2 : ℤ) : ℝ) / S ≤ x := by
  unfold Mem at h
  rw [div_le_iff₀ S_posR]
  push_cast
  linarith [(abs_le.1 h).1]

lemma abs_le_of_mem {x : ℝ} {b : Bl} (h : Mem x b) : |x| ≤ (((|b.1| + b.2) / S + 1 : ℤ) : ℝ) := by
  unfold Mem at h
  have k := lt_ediv (|b.1| + b.2) S_pos
  push_cast at k ⊢
  have hS := S_posR
  have h1 : |x * S| ≤ |(b.1 : ℝ)| + b.2 := by
    have := abs_sub_abs_le_abs_sub (x * S) (b.1 : ℝ)
    linarith
  rw [abs_mul, abs_of_pos hS] at h1
  have h2 : |x| ≤ (|(b.1 : ℝ)| + b.2) / S := by rw [le_div_iff₀ hS]; exact h1
  linarith

/-! ### Taylor loops for `cos` and `sin` -/

/-- After `j` steps: `((−1)^j y^{2j}/(2j)!, cosT (j+1) y)`, from a ball `z ∋ y²`. -/
def cosLoop (z : Bl) : ℕ → Bl × Bl
  | 0 => ((S, 0), (S, 0))
  | j + 1 => (divN (((2 * j + 1) * (2 * j + 2) : ℕ) : ℤ) (neg (mul (cosLoop z j).1 z)),
      add (cosLoop z j).2 (divN (((2 * j + 1) * (2 * j + 2) : ℕ) : ℤ) (neg (mul (cosLoop z j).1 z))))

/-- After `j` steps: `((−1)^j y^{2j+1}/(2j+1)!, sinT (j+1) y)`. -/
def sinLoop (y z : Bl) : ℕ → Bl × Bl
  | 0 => (y, y)
  | j + 1 => (divN (((2 * j + 2) * (2 * j + 3) : ℕ) : ℤ) (neg (mul (sinLoop y z j).1 z)),
      add (sinLoop y z j).2 (divN (((2 * j + 2) * (2 * j + 3) : ℕ) : ℤ) (neg (mul (sinLoop y z j).1 z))))

lemma mem_one : Mem 1 (S, 0) := mem_of_eq (by ring)

lemma cosLoop_sound {y : ℝ} {z : Bl} (hz : Mem (y ^ 2) z) : ∀ j : ℕ,
    Mem ((-1) ^ j * y ^ (2 * j) / ((2 * j).factorial : ℝ)) (cosLoop z j).1 ∧
    Mem (cosT (j + 1) y) (cosLoop z j).2
  | 0 => by
    simp only [cosLoop, cosT, pow_zero, mul_zero, Nat.factorial_zero, Nat.cast_one, div_one,
      one_mul, zero_add, Finset.sum_range_one]
    exact ⟨mem_one, mem_one⟩
  | j + 1 => by
    obtain ⟨h1, h2⟩ := cosLoop_sound hz j
    have hq : (0 : ℤ) < (((2 * j + 1) * (2 * j + 2) : ℕ) : ℤ) := by positivity
    have ht := mem_divN hq (mem_neg (mem_mul h1 hz))
    have e : -((-1) ^ j * y ^ (2 * j) / ((2 * j).factorial : ℝ) * y ^ 2)
        / ((((2 * j + 1) * (2 * j + 2) : ℕ) : ℤ) : ℝ)
        = (-1) ^ (j + 1) * y ^ (2 * (j + 1)) / ((2 * (j + 1)).factorial : ℝ) := by
      rw [show 2 * (j + 1) = 2 * j + 1 + 1 by ring, Nat.factorial_succ, Nat.factorial_succ]
      push_cast
      field_simp
      ring
    rw [e] at ht
    refine ⟨ht, ?_⟩
    have hs := mem_add h2 ht
    have e2 : cosT (j + 1) y + (-1) ^ (j + 1) * y ^ (2 * (j + 1)) / ((2 * (j + 1)).factorial : ℝ)
        = cosT (j + 1 + 1) y := by
      unfold cosT
      rw [Finset.sum_range_succ (n := j + 1)]
    rw [e2] at hs
    exact hs

lemma sinLoop_sound {y : ℝ} {yb z : Bl} (hy : Mem y yb) (hz : Mem (y ^ 2) z) : ∀ j : ℕ,
    Mem ((-1) ^ j * y ^ (2 * j + 1) / ((2 * j + 1).factorial : ℝ)) (sinLoop yb z j).1 ∧
    Mem (sinT (j + 1) y) (sinLoop yb z j).2
  | 0 => by
    simp only [sinLoop, sinT, pow_zero, mul_zero, zero_add, pow_one, Nat.factorial_one, Nat.cast_one,
      div_one, one_mul, Finset.sum_range_one]
    exact ⟨hy, hy⟩
  | j + 1 => by
    obtain ⟨h1, h2⟩ := sinLoop_sound hy hz j
    have hq : (0 : ℤ) < (((2 * j + 2) * (2 * j + 3) : ℕ) : ℤ) := by positivity
    have ht := mem_divN hq (mem_neg (mem_mul h1 hz))
    have e : -((-1) ^ j * y ^ (2 * j + 1) / ((2 * j + 1).factorial : ℝ) * y ^ 2)
        / ((((2 * j + 2) * (2 * j + 3) : ℕ) : ℤ) : ℝ)
        = (-1) ^ (j + 1) * y ^ (2 * (j + 1) + 1) / ((2 * (j + 1) + 1).factorial : ℝ) := by
      rw [show 2 * (j + 1) + 1 = 2 * j + 1 + 1 + 1 by ring, Nat.factorial_succ (2 * j + 1 + 1),
        Nat.factorial_succ (2 * j + 1)]
      push_cast
      field_simp
      ring
    rw [e] at ht
    refine ⟨ht, ?_⟩
    have hs := mem_add h2 ht
    have e2 : sinT (j + 1) y + (-1) ^ (j + 1) * y ^ (2 * (j + 1) + 1) / ((2 * (j + 1) + 1).factorial : ℝ)
        = sinT (j + 1 + 1) y := by
      unfold sinT
      rw [Finset.sum_range_succ (n := j + 1)]
    rw [e2] at hs
    exact hs

lemma re_im_partial (y : ℝ) : ∀ N : ℕ,
    (∑ m ∈ range (2 * N), ((y : ℂ) * Complex.I) ^ m / (m.factorial : ℂ)).re = cosT N y ∧
    (∑ m ∈ range (2 * N), ((y : ℂ) * Complex.I) ^ m / (m.factorial : ℂ)).im = sinT N y
  | 0 => by simp [cosT, sinT]
  | N + 1 => by
    obtain ⟨h1, h2⟩ := re_im_partial y N
    rw [show 2 * (N + 1) = 2 * N + 1 + 1 by ring, Finset.sum_range_succ, Finset.sum_range_succ]
    have e1 : ((y : ℂ) * Complex.I) ^ (2 * N) = (((-1) ^ N * y ^ (2 * N) : ℝ) : ℂ) := by
      rw [mul_pow, pow_mul Complex.I, Complex.I_sq]
      push_cast
      ring
    have e2 : ((y : ℂ) * Complex.I) ^ (2 * N + 1)
        = (((-1) ^ N * y ^ (2 * N + 1) : ℝ) : ℂ) * Complex.I := by
      rw [pow_succ, e1]
      push_cast
      ring
    have f1 : ((2 * N).factorial : ℂ) = (((2 * N).factorial : ℝ) : ℂ) := by norm_cast
    have f2 : ((2 * N + 1).factorial : ℂ) = (((2 * N + 1).factorial : ℝ) : ℂ) := by norm_cast
    rw [e1, e2, f1, f2]
    unfold cosT sinT
    rw [Finset.sum_range_succ, Finset.sum_range_succ]
    refine ⟨?_, ?_⟩
    · rw [Complex.add_re, Complex.add_re, h1, Complex.div_ofReal_re, Complex.div_ofReal_re,
        Complex.ofReal_re, Complex.re_ofReal_mul, Complex.I_re]
      unfold cosT
      ring
    · rw [Complex.add_im, Complex.add_im, h2, Complex.div_ofReal_im, Complex.div_ofReal_im,
        Complex.ofReal_im, Complex.im_ofReal_mul, Complex.I_im]
      unfold sinT
      ring

lemma cos_sin_rem (N : ℕ) {y : ℝ} (hy : |y| ≤ N) :
    |Real.cos y - cosT N y| ≤ 2 * |y| ^ (2 * N) / ((2 * N).factorial : ℝ) ∧
    |Real.sin y - sinT N y| ≤ 2 * |y| ^ (2 * N) / ((2 * N).factorial : ℝ) := by
  have hn : ‖(y : ℂ) * Complex.I‖ / ((2 * N).succ : ℝ) ≤ 1 / 2 := by
    rw [norm_mul, Complex.norm_I, mul_one, Complex.norm_real, Real.norm_eq_abs,
      div_le_iff₀ (by positivity)]
    push_cast
    linarith
  have hb := Complex.exp_bound' hn
  rw [norm_mul, Complex.norm_I, mul_one, Complex.norm_real, Real.norm_eq_abs] at hb
  obtain ⟨hre, him⟩ := re_im_partial y N
  set P := ∑ m ∈ range (2 * N), ((y : ℂ) * Complex.I) ^ m / (m.factorial : ℂ)
  have k1 : |Real.cos y - cosT N y| ≤ ‖Complex.exp ((y : ℂ) * Complex.I) - P‖ := by
    have := Complex.abs_re_le_norm (Complex.exp ((y : ℂ) * Complex.I) - P)
    rw [Complex.sub_re, hre, Complex.exp_ofReal_mul_I_re] at this
    exact this
  have k2 : |Real.sin y - sinT N y| ≤ ‖Complex.exp ((y : ℂ) * Complex.I) - P‖ := by
    have := Complex.abs_im_le_norm (Complex.exp ((y : ℂ) * Complex.I) - P)
    rw [Complex.sub_im, him, Complex.exp_ofReal_mul_I_im] at this
    exact this
  have k3 : ‖Complex.exp ((y : ℂ) * Complex.I) - P‖ ≤ 2 * |y| ^ (2 * N) / ((2 * N).factorial : ℝ) := by
    refine hb.trans (le_of_eq ?_)
    ring
  exact ⟨k1.trans k3, k2.trans k3⟩

/-- Taylor order. -/
def NT : ℕ := 16

/-- Remainder radius for `|y| ≤ Y`. -/
def remZ (Y : ℤ) : ℤ := 2 * S * Y ^ (2 * NT) / ((2 * NT).factorial : ℤ) + 1

lemma rem_le {y : ℝ} {Y : ℤ} (hY : |y| ≤ Y) :
    2 * |y| ^ (2 * NT) / ((2 * NT).factorial : ℝ) * S ≤ (remZ Y : ℝ) := by
  have hf : (0 : ℤ) < ((2 * NT).factorial : ℤ) := by positivity
  have k := lt_ediv (2 * S * Y ^ (2 * NT)) hf
  unfold remZ
  push_cast at k ⊢
  have hS := S_posR
  have hp : |y| ^ (2 * NT) ≤ (Y : ℝ) ^ (2 * NT) := pow_le_pow_left₀ (abs_nonneg y) hY _
  have hF : (0 : ℝ) < ((2 * NT).factorial : ℝ) := by positivity
  have : 2 * |y| ^ (2 * NT) / ((2 * NT).factorial : ℝ) * S ≤ 2 * S * (Y : ℝ) ^ (2 * NT) / ((2 * NT).factorial : ℝ) := by
    rw [div_mul_eq_mul_div, div_le_div_iff_of_pos_right hF]
    nlinarith
  linarith

/-- `cos` and `sin` balls from a ball `yb ∋ y` with `|y|` small. -/
def cosB (yb : Bl) : Bl :=
  if (|yb.1| + yb.2) / S + 1 ≤ (NT : ℤ) then
    add (cosLoop (mul yb yb) (NT - 1)).2 (0, remZ ((|yb.1| + yb.2) / S + 1))
  else (0, S)

def sinB (yb : Bl) : Bl :=
  if (|yb.1| + yb.2) / S + 1 ≤ (NT : ℤ) then
    add (sinLoop yb (mul yb yb) (NT - 1)).2 (0, remZ ((|yb.1| + yb.2) / S + 1))
  else (0, S)

lemma mem_rem {x c : ℝ} {b : Bl} {R : ℤ} (hb : Mem c b) (h : |x - c| * S ≤ R) :
    Mem x (add b (0, R)) := by
  unfold Mem add at *
  push_cast
  have hS := S_posR
  calc |x * S - (b.1 + 0)| = |(x - c) * S + (c * S - b.1)| := by congr 1; ring
    _ ≤ |(x - c) * S| + |c * S - b.1| := abs_add_le _ _
    _ ≤ R + b.2 := by rw [abs_mul, abs_of_pos hS]; linarith
    _ = b.2 + R := by ring

lemma mem_unit {x : ℝ} (h : |x| ≤ 1) : Mem x (0, S) := by
  unfold Mem
  push_cast
  rw [sub_zero, abs_mul, abs_of_pos S_posR]
  nlinarith [S_posR]

theorem mem_cosB {y : ℝ} {yb : Bl} (hy : Mem y yb) : Mem (Real.cos y) (cosB yb) := by
  unfold cosB
  split_ifs with h
  · have hY := abs_le_of_mem hy
    have hN : |y| ≤ (NT : ℝ) := by
      have : ((((|yb.1| + yb.2) / S + 1 : ℤ)) : ℝ) ≤ ((NT : ℤ) : ℝ) := by exact_mod_cast h
      push_cast at this hY
      linarith
    have hz : Mem (y ^ 2) (mul yb yb) := by rw [sq]; exact mem_mul hy hy
    have hl := (cosLoop_sound hz (NT - 1)).2
    rw [show NT - 1 + 1 = NT by norm_num [NT]] at hl
    refine mem_rem hl ?_
    have hr := (cos_sin_rem NT hN).1
    have := rem_le hY
    have hS := S_posR
    nlinarith
  · exact mem_unit (Real.abs_cos_le_one y)

theorem mem_sinB {y : ℝ} {yb : Bl} (hy : Mem y yb) : Mem (Real.sin y) (sinB yb) := by
  unfold sinB
  split_ifs with h
  · have hY := abs_le_of_mem hy
    have hN : |y| ≤ (NT : ℝ) := by
      have : ((((|yb.1| + yb.2) / S + 1 : ℤ)) : ℝ) ≤ ((NT : ℤ) : ℝ) := by exact_mod_cast h
      push_cast at this hY
      linarith
    have hz : Mem (y ^ 2) (mul yb yb) := by rw [sq]; exact mem_mul hy hy
    have hl := (sinLoop_sound hy hz (NT - 1)).2
    rw [show NT - 1 + 1 = NT by norm_num [NT]] at hl
    refine mem_rem hl ?_
    have hr := (cos_sin_rem NT hN).2
    have := rem_le hY
    have hS := S_posR
    nlinarith
  · exact mem_unit (Real.abs_sin_le_one y)

/-! ### Constants -/

/-- `log 2` (from `Real.log_two_near_10`). -/
def lB : Bl := (3196577161290467479, 461168603)

theorem mem_log2 : Mem (Real.log 2) lB := by
  have h := Real.log_two_near_10
  unfold Mem lB
  have hq : |(287209 / 414355 : ℝ) * S - 3196577161290467479| ≤ 1 := by
    norm_num [S, abs_le]
  have e : Real.log 2 * S - 3196577161290467479
      = (Real.log 2 - 287209 / 414355) * S + ((287209 / 414355 : ℝ) * S - 3196577161290467479) := by ring
  push_cast
  rw [e]
  have hS := S_posR
  calc _ ≤ |(Real.log 2 - 287209 / 414355) * S| + |(287209 / 414355 : ℝ) * S - 3196577161290467479| :=
        abs_add_le _ _
    _ ≤ 1 / 10 ^ 10 * S + 1 := by
        rw [abs_mul, abs_of_pos hS]
        exact add_le_add (mul_le_mul_of_nonneg_right h hS.le) hq
    _ ≤ 461168603 := by norm_num [S]

/-- `π` (20 digits). -/
def piB : Bl := (14488038916154245684, 2)

theorem mem_pi : Mem Real.pi piB := by
  have h1 := Real.pi_gt_d20
  have h2 := Real.pi_lt_d20
  unfold Mem piB
  push_cast
  rw [abs_le]
  constructor <;> nlinarith [S_posR, show (S : ℝ) = 4611686018427387904 by norm_num [S]]

/-- `√2`. -/
def s2B : Bl := (6521908912666391106, 1)

theorem mem_sqrt2 : Mem (Real.sqrt 2) s2B := by
  unfold Mem s2B
  have hS : (S : ℝ) = 4611686018427387904 := by norm_num [S]
  rw [hS]
  push_cast
  have hlo : (6521908912666391106 : ℝ) / 4611686018427387904 ≤ Real.sqrt 2 := by
    rw [Real.le_sqrt (by norm_num) (by norm_num)]
    norm_num
  have hhi : Real.sqrt 2 ≤ (6521908912666391107 : ℝ) / 4611686018427387904 := by
    rw [Real.sqrt_le_left (by norm_num)]
    norm_num
  rw [abs_le]
  constructor
  · have := mul_le_mul_of_nonneg_right hlo (by norm_num : (0 : ℝ) ≤ 4611686018427387904)
    rw [div_mul_cancel₀ _ (by norm_num)] at this
    linarith
  · have := mul_le_mul_of_nonneg_right hhi (by norm_num : (0 : ℝ) ≤ 4611686018427387904)
    rw [div_mul_cancel₀ _ (by norm_num)] at this
    linarith

def cPB : Bl := mul s2B lB

theorem mem_cP : Mem Beyond.cP cPB := mem_mul mem_sqrt2 mem_log2

/-! ### The digamma lower bound -/

/-- `Σ_{n<M}` lower bounds of `S·fT(i/2^e, n)`. -/
def muSumZ (e i : ℕ) : ℕ → ℤ
  | 0 => 0
  | n + 1 => muSumZ e i n + (S / ((n : ℤ) + 1)
      + (-(S * (4 * (4 * (n : ℤ) + 5) * 4 ^ e))) / ((4 * (n : ℤ) + 5) ^ 2 * 4 ^ e + 4 * (i : ℤ) ^ 2))

lemma fT_grid (e i n : ℕ) :
    S * fT ((i : ℝ) / 2 ^ e) n = S / ((n : ℝ) + 1)
      - S * (4 * (4 * (n : ℝ) + 5) * 4 ^ e) / ((4 * (n : ℝ) + 5) ^ 2 * 4 ^ e + 4 * (i : ℝ) ^ 2) := by
  unfold fT
  have h4 : (0 : ℝ) < 4 ^ e := by positivity
  have hD : (0 : ℝ) < (4 * (n : ℝ) + 5) ^ 2 * 4 ^ e + 4 * (i : ℝ) ^ 2 := by positivity
  have e4 : (4 : ℝ) ^ e = (2 ^ e) ^ 2 := by rw [← pow_mul, mul_comm, pow_mul]; norm_num
  rw [e4] at hD ⊢
  field_simp
  ring

lemma muSumZ_le (e i : ℕ) : ∀ M : ℕ,
    (muSumZ e i M : ℝ) ≤ S * ∑ n ∈ range M, fT ((i : ℝ) / 2 ^ e) n
  | 0 => by simp [muSumZ]
  | n + 1 => by
    rw [muSumZ, Finset.sum_range_succ, mul_add (S : ℝ), fT_grid]
    have ih := muSumZ_le e i n
    have hq1 : (0 : ℤ) < (n : ℤ) + 1 := by positivity
    have hq2 : (0 : ℤ) < (4 * (n : ℤ) + 5) ^ 2 * 4 ^ e + 4 * (i : ℤ) ^ 2 := by positivity
    have k1 := ediv_le S hq1
    have k2 := ediv_le (-(S * (4 * (4 * (n : ℤ) + 5) * 4 ^ e))) hq2
    push_cast at k1 k2 ⊢
    have e2 : -((S : ℝ) * (4 * (4 * (n : ℝ) + 5) * 4 ^ e)) / ((4 * (n : ℝ) + 5) ^ 2 * 4 ^ e + 4 * (i : ℝ) ^ 2)
        = -((S : ℝ) * (4 * (4 * (n : ℝ) + 5) * 4 ^ e) / ((4 * (n : ℝ) + 5) ^ 2 * 4 ^ e + 4 * (i : ℝ) ^ 2)) := by
      ring
    rw [e2] at k2
    linarith

/-- Integer lower bound of `S · 2πμ(i/2^e)` (`γ < 0.57725`, `log π < 1.14474`). -/
def muLowZ (e i M Mp : ℕ) : ℤ :=
  (-(172199 * S)) / 100000 + (-(S * (4 * 4 ^ e))) / (4 ^ e + 4 * (i : ℤ) ^ 2) + muSumZ e i M
    + (S / (4 * (M : ℤ) + 5) + (-S) / (4 * (Mp : ℤ) + 5))
    + (if (4 : ℤ) ^ (e + 1) ≤ 2 * (i : ℤ) ^ 2 then
        (S * (8 * (i : ℤ) ^ 2)) / ((4 * (M : ℤ) + 3) ^ 2 * 4 ^ (e + 1) + 16 * (i : ℤ) ^ 2)
        + (-(S * (8 * (i : ℤ) ^ 2))) / ((4 * (Mp : ℤ) + 3) ^ 2 * 4 ^ (e + 1) + 16 * (i : ℤ) ^ 2)
      else 0)

theorem muLowZ_le (e i M Mp : ℕ) (hM : M ≤ Mp) :
    ((muLowZ e i M Mp : ℤ) : ℝ) ≤ S * (2 * Real.pi * mu ((i : ℝ) / 2 ^ e)) := by
  have hg : Real.eulerMascheroniConstant ≤ ((57725 / 100000 : ℚ) : ℝ) := by
    have := euler_lt'
    push_cast
    norm_num at this ⊢
    linarith
  have hl : Real.log Real.pi ≤ ((114474 / 100000 : ℚ) : ℝ) := by
    have := AReflect.log_pi_lt'
    push_cast
    norm_num at this ⊢
    linarith
  have ht := two_pi_mu_ge_tel ((i : ℝ) / 2 ^ e) M Mp hM hg hl
  have hS := S_posR
  have hs := muSumZ_le e i M
  set r : ℝ := (i : ℝ) / 2 ^ e with hr
  have e4 : (4 : ℝ) ^ (e + 1) = ((2 : ℝ) ^ e) ^ 2 * 4 := by
    rw [pow_succ, ← pow_mul, mul_comm e 2, pow_mul]
    norm_num
  have hr2 : (r / 2) ^ 2 = (i : ℝ) ^ 2 / 4 ^ (e + 1) := by
    rw [hr, e4]
    field_simp
    ring
  have h4 : (0 : ℝ) < 4 ^ e := by positivity
  -- head terms
  have q1 : (0 : ℤ) < 100000 := by norm_num
  have q2 : (0 : ℤ) < 4 ^ e + 4 * (i : ℤ) ^ 2 := by positivity
  have q3 : (0 : ℤ) < 4 * (M : ℤ) + 5 := by positivity
  have q4 : (0 : ℤ) < 4 * (Mp : ℤ) + 5 := by positivity
  have k1 := ediv_le (-(172199 * S)) q1
  have k2 := ediv_le (-(S * (4 * 4 ^ e))) q2
  have k3 := ediv_le S q3
  have k4 := ediv_le (-S) q4
  have a1 : (S : ℝ) * (1 / 4 / ((1 / 4) ^ 2 + (r / 2) ^ 2)) = S * (4 * 4 ^ e) / (4 ^ e + 4 * (i : ℝ) ^ 2) := by
    rw [hr2, pow_succ]
    field_simp
    ring
  have a3 : (S : ℝ) * h1R M = S / (4 * (M : ℝ) + 5) := by unfold h1R; field_simp
  have a4 : (S : ℝ) * h1R Mp = S / (4 * (Mp : ℝ) + 5) := by unfold h1R; field_simp
  have key : ((muLowZ e i M Mp : ℤ) : ℝ) ≤ S * (-((57725 / 100000 : ℚ) : ℝ) - ((114474 / 100000 : ℚ) : ℝ)
      - (1 / 4) / ((1 / 4) ^ 2 + (r / 2) ^ 2) + ∑ n ∈ range M, fT r n
      + (h1R M - h1R Mp + (if 1 / 2 ≤ (r / 2) ^ 2 then h2R ((r / 2) ^ 2) M - h2R ((r / 2) ^ 2) Mp
        else 0))) := by
    have hcond : ((4 : ℤ) ^ (e + 1) ≤ 2 * (i : ℤ) ^ 2) ↔ (1 / 2 ≤ (r / 2) ^ 2) := by
      rw [hr2, le_div_iff₀ (by positivity)]
      constructor
      · intro h
        have : ((4 ^ (e + 1) : ℤ) : ℝ) ≤ ((2 * (i : ℤ) ^ 2 : ℤ) : ℝ) := by exact_mod_cast h
        push_cast at this
        linarith
      · intro h
        have : ((4 : ℝ) ^ (e + 1)) ≤ 2 * (i : ℝ) ^ 2 := by linarith
        exact_mod_cast this
    unfold muLowZ
    push_cast at k1 k2 k3 k4 ⊢
    have base : -(172199 * (S : ℝ)) / 100000 + -((S : ℝ) * (4 * 4 ^ e)) / (4 ^ e + 4 * (i : ℝ) ^ 2)
        + ((S : ℝ) / (4 * (M : ℝ) + 5) + -(S : ℝ) / (4 * (Mp : ℝ) + 5))
        = S * (-((57725 / 100000 : ℚ) : ℝ) - ((114474 / 100000 : ℚ) : ℝ))
          - S * (1 / 4 / ((1 / 4) ^ 2 + (r / 2) ^ 2)) + S * (h1R M - h1R Mp) := by
      rw [a1, mul_sub (S : ℝ) (h1R M), a3, a4]
      push_cast
      ring
    by_cases hc : (4 : ℤ) ^ (e + 1) ≤ 2 * (i : ℤ) ^ 2
    · rw [if_pos hc, if_pos (hcond.1 hc)]
      have q5 : (0 : ℤ) < (4 * (M : ℤ) + 3) ^ 2 * 4 ^ (e + 1) + 16 * (i : ℤ) ^ 2 := by positivity
      have q6 : (0 : ℤ) < (4 * (Mp : ℤ) + 3) ^ 2 * 4 ^ (e + 1) + 16 * (i : ℤ) ^ 2 := by positivity
      have k5 := ediv_le (S * (8 * (i : ℤ) ^ 2)) q5
      have k6 := ediv_le (-(S * (8 * (i : ℤ) ^ 2))) q6
      push_cast at k5 k6
      have a5 : ∀ n : ℕ, (S : ℝ) * h2R ((r / 2) ^ 2) n
          = S * (8 * (i : ℝ) ^ 2) / ((4 * (n : ℝ) + 3) ^ 2 * 4 ^ (e + 1) + 16 * (i : ℝ) ^ 2) := by
        intro n
        unfold h2R
        rw [hr2]
        have : (0 : ℝ) < 4 ^ (e + 1) := by positivity
        field_simp
        ring
      have e6 : -((S : ℝ) * (8 * (i : ℝ) ^ 2)) / ((4 * (Mp : ℝ) + 3) ^ 2 * 4 ^ (e + 1) + 16 * (i : ℝ) ^ 2)
          = -((S : ℝ) * h2R ((r / 2) ^ 2) Mp) := by rw [a5]; ring
      rw [e6] at k6
      rw [← a5] at k5
      nlinarith
    · rw [if_neg hc, if_neg (fun h => hc (hcond.2 h))]
      nlinarith
  calc ((muLowZ e i M Mp : ℤ) : ℝ) ≤ _ := key
    _ ≤ S * (2 * Real.pi * mu r) := mul_le_mul_of_nonneg_left ht hS.le

end Fayed.Yoshida.Beyond.Ball
