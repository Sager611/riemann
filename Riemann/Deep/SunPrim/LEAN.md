PARTIAL: 11 of the 13 table rows of discoveries/sun-counterexamples.md are refuted in Lean 4 + Mathlib (all the finite ones: 1309.1679 Conj. 3.2, 3.4(i), 3.10(ii); 1811.10503 Conj. 4.9 (n = 3, 4 in F₁₁, n = 4 in F₁₃); MathOverflow 316143; OEIS A228766, A228728; Wang–Sun 2210.12044 Conj. 1.1 (C); A229827(i) in F₈ and F₉; 1208.2683 Conj. 3.1; 1211.1588 Conj. 4.30(iii)), each from the source's own statement; left: 1208.2683 Conj. 2.1 (real analysis: p^α for α ≈ 10⁻³) and 1402.6641 Conj. 2.25 (prime ideals of ℤ[i]); the two ⁿ√ conjectures use an abstract root function, since Mathlib's real powers exceed the memory cap.

# Lean refutations of Sun's conjectures (sunlean, 2026-11-12)

Files: `Riemann/Deep/SunPrim/Lean/` (namespace `SunLean`). One `Prop` per conjecture (`Conj_*`), one theorem `not_conj_* : ¬ Conj_*` per row. Sorry-free; axioms `propext, Classical.choice, Quot.sound` only; no `native_decide`, no `skipKernelTC`. Finite facts by `decide` / `decide +kernel` (kernel evaluation, standard axioms). Sources: arXiv TeX e-prints (latest versions: 1309.1679v9 = `196a.tex`, 1811.10503v6, 2210.12044v1, 1211.1588v29, 1208.2683v9), OEIS text format (`https://oeis.org/search?q=id:A228766&fmt=text`), Stack Exchange API (questions 316142, 316143); downloads deleted after use.

## Rows, statements, theorems

| Row | Source (line) | Lean statement | Theorem | Instance |
|---|---|---|---|---|
| 1309.1679 Conj. 3.2 | TeX label Conj3.3, l. 512–516 | `Conj_3_2` (Groups) | `not_conj_3_2`, `not_conj_3_2_with` (any Sylow predicate) | ℤ/6, (0, 1, 3, 5) |
| 1309.1679 Conj. 3.4(i) | Conj3.5, l. 541–548 | `Conj_3_4_i` (Groups) | `not_conj_3_4_i` | ℤ/7, {0, …, 4} |
| 1309.1679 Conj. 3.10(ii) | Conj3.12, l. 667–672 | `Conj_3_10_ii R` (Conj310) | `not_conj_3_10_ii` (every char-0 field), `not_conj_3_10_ii_real` | {2, 3, 4, 6} |
| 1811.10503 Conj. 4.9 | l. 767–770 | `Conj_4_9` (Conj49) | `not_conj_4_9_n3` (Conj49), `_n4_F11` (Conj49b), `_n4_F13` (Conj49c) | F₁₁ {0, 1, 2, 5}; F₁₃ {0, 1, 3, 4} |
| MathOverflow 316143 | question body (2018-11-25) | `Conj_MO316143` (MO316143) | `not_conj_MO316143` | λ = (1, 1, 2), {0, 1, 2, 3} ⊂ F₁₁ |
| OEIS A228766 | %C, 2nd sentence | `Conj_A228766` (Groups) | `not_conj_A228766` | ℤ/6, {0, 1, 3, 4} |
| OEIS A228728 | %C "Conjecture: …" | `Conj_A228728 R` (A228728) | `not_conj_A228728` (every ordered field), `_real` | (0, 1, 2, 3, 4, 5, 8) |
| Wang–Sun 2210.12044 Conj. 1.1 | label ConjG, l. 174–184; p(G) l. 153 | `Conj_WS`, `Conj_WS_C` (WangSun) | `not_conj_WS`, `not_conj_WS_C` | ℤ/7, ({0,1}, {0,1,2,4,5}, {2,4,5}, {0,1,4,5}) |
| OEIS A229827 (i) | %C "Moreover, …" | `Conj_A229827_i` (A229827) | `not_conj_A229827_i_F8`, `_F9` | F₈, x²; F₉, x³ |
| 1208.2683 Conj. 3.1 | l. 462–470 | `Conj_3_1 K rt` (Roots) | `not_conj_3_1`, `_real` | (A, B) = (3, −26), n = 5 |
| 1211.1588 Conj. 4.30(iii) | l. 1339–1340 (p^{(m)}), 1349–1352 | `Conj_4_30_iii K rt` (Roots) | `not_conj_4_30_iii`, `_real` | m = 1, n = 2 |

Skipped: 1208.2683 Conj. 2.1 (needs (1/9)Σ p_k^α vs 23^α/(α+1) at α ≈ 0.001 to 10⁻⁸: real powers, which do not load under the cap, plus certified bounds) and 1402.6641 Conj. 2.25 (counts of prime ideals of ℤ[i]; no light faithful encoding).

## Statements next to the sources

**1309.1679 Conj. 3.2** (l. 512–516): "Let a_1, …, a_n be n distinct elements of a finite additive abelian group G. Suppose that n ∤ |G|, or n is even and the Sylow 2-subgroup of G is cyclic. Then there exists a permutation (b_1, …, b_n) of a_1, …, a_n with b_1 = a_1 such that the n − 1 elements b_i − b_{i+1} (0 < i < n) are pairwise distinct."
```lean
def Syl2Cyclic (G : Type) [AddCommGroup G] : Prop :=
  ∃ g : G, (∃ k : ℕ, 2 ^ k • g = 0) ∧ ∀ h : G, (∃ k : ℕ, 2 ^ k • h = 0) → ∃ m : ℕ, h = m • g
def Conj_3_2_with (syl : (G : Type) → [AddCommGroup G] → Prop) : Prop :=
  ∀ (G : Type) [AddCommGroup G] [Fintype G] (a : List G), a.Nodup →
    (¬ a.length ∣ Fintype.card G ∨ (Even a.length ∧ syl G)) →
    ∃ b : List G, b.Perm a ∧ b.head? = a.head? ∧ (List.zipWith (· - ·) b b.tail).Nodup
def Conj_3_2 : Prop := Conj_3_2_with Syl2Cyclic
```

**1309.1679 Conj. 3.4(i)** (l. 541–548): "Let A be a finite subset of an additive abelian group G with |A| = n > 3. (i) If G is finite with |G| ≢ 0 (mod 3), then there is a numbering a_1, …, a_n of all the elements of A such that the n sums a_1 + 2a_2, a_2 + 2a_3, …, a_{n−1} + 2a_n, a_n + 2a_1 are pairwise distinct."
```lean
def Conj_3_4_i : Prop :=
  ∀ (G : Type) [AddCommGroup G] [Fintype G] (A : Finset G), 3 < A.card → ¬ 3 ∣ Fintype.card G →
    ∃ l : List G, l.Nodup ∧ (∀ x, x ∈ l ↔ x ∈ A) ∧
      (List.zipWith (fun x y => x + 2 • y) l (l.rotate 1)).Nodup
```

**1309.1679 Conj. 3.10(ii)** (l. 667–672): "For any set A of n > 3 distinct nonzero real numbers, there is a circular permutation (a_1, …, a_n) of all the elements of A such that the n adjacent differences a_1 − a_2, …, a_{n−1} − a_n, a_n − a_1 are pairwise distinct, and that the n adjacent products a_1a_2, …, a_{n−1}a_n, a_na_1 are also pairwise distinct, unless |A| = 4 and A has the form {±s, ±t}."
```lean
def Conj_3_10_ii (R : Type) [Field R] : Prop :=
  ∀ A : Finset R, 3 < A.card → (∀ x ∈ A, x ≠ 0) →
    ¬ (A.card = 4 ∧ ∃ s t : R, (A : Set R) = {s, -s, t, -t}) →
    ∃ l : List R, l.Nodup ∧ (∀ x, x ∈ l ↔ x ∈ A) ∧
      (List.zipWith (· - ·) l (l.rotate 1)).Nodup ∧ (List.zipWith (· * ·) l (l.rotate 1)).Nodup
-- not_conj_3_10_ii (R) [Field R] [CharZero R] : ¬ Conj_3_10_ii R;  not_conj_3_10_ii_real : ¬ Conj_3_10_ii ℝ
```

**1811.10503 Conj. 4.9** (l. 767–770): "Let n ∈ ℤ⁺ and let F be a field with p(F) > n + 1, where p(F) = p if the characteristic of F is a prime p, and p(F) = +∞ if the characteristic of F is zero. Let A be any finite subset of F with |A| ≥ n + δ_{n,3} … Then, for the set S(A) := {Σ_{k=1}^n k a_k : a_1, …, a_n are distinct elements of A}, we have |S(A)| ≥ min{p(F), (|A| − n) n(n+1)/2 + n(n² − 1)/6 + 1}."
```lean
noncomputable def pF (F : Type) [Field F] : ℕ∞ := if ringChar F = 0 then ⊤ else (ringChar F : ℕ∞)
def WSumSet {F : Type} [Field F] (n : ℕ) (w : ℕ → ℕ) (A : Finset F) : Set F :=
  {s | ∃ a : ℕ → F, (∀ k ∈ Icc 1 n, a k ∈ A) ∧ (∀ i ∈ Icc 1 n, ∀ j ∈ Icc 1 n, a i = a j → i = j) ∧
    s = ∑ k ∈ Icc 1 n, (w k : F) * a k}
def Conj_4_9 : Prop :=
  ∀ n : ℕ, 0 < n → ∀ (F : Type) [Field F], ((n + 1 : ℕ) : ℕ∞) < pF F →
    ∀ A : Finset F, n + (if n = 3 then 1 else 0) ≤ A.card →
      min (pF F) (((A.card - n) * (n * (n + 1) / 2) + n * (n ^ 2 - 1) / 6 + 1 : ℕ) : ℕ∞) ≤
        (WSumSet n (fun k => k) A).encard
```

**MathOverflow 316143** (body, "Conjecture."): "Let λ_1, …, λ_n (n ≥ 3) be positive integers with gcd(λ_1, …, λ_n) = 1 and λ_{k+1} − λ_k ∈ {0, 1} for all k = 1, …, n − 1. Let F be a field with p(F) > n + 1. Then, for any finite subset A of F with |A| ≥ n + δ_{n,3} we have |{Σ_{k=1}^n λ_k a_k : a_1, …, a_n are distinct elements of A}| ≥ min{p(F), (λ_1 + … + λ_n)(|A| − n) + Σ_{k=1}^{⌊n/2⌋} (n + 1 − 2k)(λ_{n+1−k} − λ_k) + 1}."
```lean
def bound316143 (n : ℕ) (lam : ℕ → ℕ) (m : ℕ) : ℤ :=
  (∑ k ∈ Icc 1 n, (lam k : ℤ)) * ((m : ℤ) - n) +
    (∑ k ∈ Icc 1 (n / 2), ((n : ℤ) + 1 - 2 * k) * ((lam (n + 1 - k) : ℤ) - lam k)) + 1
def Conj_MO316143 : Prop :=
  ∀ (n : ℕ) (lam : ℕ → ℕ), 3 ≤ n → (∀ k ∈ Icc 1 n, 0 < lam k) →
    (∀ d : ℕ, (∀ k ∈ Icc 1 n, d ∣ lam k) → d = 1) →
    (∀ k ∈ Icc 1 (n - 1), (lam (k + 1) : ℤ) - lam k = 0 ∨ (lam (k + 1) : ℤ) - lam k = 1) →
    ∀ (F : Type) [Field F], ((n + 1 : ℕ) : ℕ∞) < pF F →
    ∀ A : Finset F, n + (if n = 3 then 1 else 0) ≤ A.card →
      min (pF F) ((bound316143 n lam A.card).toNat : ℕ∞) ≤ (WSumSet n lam A).encard
```

**OEIS A228766** (%C): "In general, if a_1, …, a_n are n > 2 distinct elements of a finite additive abelian group G with n odd or |G| not divisible by n, then there exists a circular permutation b_1, …, b_n of a_1, …, a_n such that b_1+b_2, b_2+b_3, …, b_{n−1}+b_n, b_n+b_1 are pairwise distinct."
```lean
def Conj_A228766 : Prop :=
  ∀ (G : Type) [AddCommGroup G] [Fintype G] (a : List G), 2 < a.length → a.Nodup →
    (a.length % 2 = 1 ∨ ¬ a.length ∣ Fintype.card G) →
    ∃ b : List G, b.Perm a ∧ (List.zipWith (· + ·) b (b.rotate 1)).Nodup
```

**OEIS A228728** (%C): "For any n distinct real numbers a_1 < a_2 < … < a_n, if there is a permutation b_1, …, b_n of a_1, …, a_n with |b_1−b_2|, |b_2−b_3|, …, |b_{n−1}−b_n|, |b_n−b_1| pairwise distinct, then there exists a permutation c_1, …, c_n of a_1, …, a_n with c_1 = a_1 and c_n = a_n such that the n numbers |c_1−c_2|, …, |c_{n−1}−c_n|, |c_n−c_1| are pairwise distinct."
```lean
def Conj_A228728 (R : Type) [Field R] [LinearOrder R] : Prop :=
  ∀ a : List R, a.Pairwise (· < ·) →
    (∃ b : List R, b.Perm a ∧ (List.zipWith (fun x y => |x - y|) b (b.rotate 1)).Nodup) →
    ∃ c : List R, c.Perm a ∧ c.head? = a.head? ∧ c.getLast? = a.getLast? ∧
      (List.zipWith (fun x y => |x - y|) c (c.rotate 1)).Nodup
-- not_conj_A228728 (R) [Field R] [LinearOrder R] [IsStrictOrderedRing R];  not_conj_A228728_real
```

**Wang–Sun 2210.12044 Conj. 1.1** (l. 174–184; L, C l. 161–166; p(G) l. 153): "Let G be an additive group with |G| > 1, and let A_1, …, A_n (n > 1) be finite subsets of G with |A_i| > 1 for all i. Then |L(A_1, …, A_n)| ≥ min{p(G), |A_1| + ⋯ + |A_n| − 2n + 1 + {n}_2} and |C(A_1, …, A_n)| ≥ min{p(G), |A_1| + ⋯ + |A_n| − 2n + (−1)^n (1 + {n}_2)}", L: sums with a_i ≠ a_{i+1} (0 < i < n); C: also a_n ≠ a_1; p(G): "the minimum of orders of nonzero elements of G (which is +∞ if each nonzero element of G has infinite order)".
```lean
noncomputable def pG (G : Type) [AddCommGroup G] : ℕ∞ :=
  ⨅ (g : G) (_ : g ≠ 0) (m : ℕ) (_ : 0 < m) (_ : m • g = 0), (m : ℕ∞)
def Lset {G : Type} [AddCommGroup G] (n : ℕ) (A : Fin n → Finset G) : Set G :=
  {s | ∃ a : Fin n → G, (∀ i, a i ∈ A i) ∧ (∀ i j : Fin n, j.val = i.val + 1 → a i ≠ a j) ∧
    s = ∑ i, a i}
def Cset {G : Type} [AddCommGroup G] (n : ℕ) (A : Fin n → Finset G) : Set G :=
  {s | ∃ a : Fin n → G, (∀ i, a i ∈ A i) ∧ (∀ i j : Fin n, j.val = i.val + 1 → a i ≠ a j) ∧
    (∀ i j : Fin n, i.val + 1 = n → j.val = 0 → a i ≠ a j) ∧ s = ∑ i, a i}
def Conj_WS : Prop :=
  ∀ (G : Type) [AddCommGroup G] [Nontrivial G] (n : ℕ), 1 < n → ∀ A : Fin n → Finset G,
    (∀ i, 1 < (A i).card) →
      min (pG G) (((∑ i, ((A i).card : ℤ)) - 2 * n + 1 + (n % 2 : ℕ)).toNat : ℕ∞) ≤ (Lset n A).encard ∧
      min (pG G) (((∑ i, ((A i).card : ℤ)) - 2 * n + (-1) ^ n * (1 + (n % 2 : ℕ))).toNat : ℕ∞) ≤
        (Cset n A).encard
-- Conj_WS_C: the second inequality alone
```

**OEIS A229827 (i)** (%C): "Moreover, for any finite field F(q) with q elements and a polynomial P(x) over F(q) of degree smaller than q−1, if P(x) is not of the form c−x, then there is a circular permutation a_1, …, a_q of all the elements of F(q) with {P(a_1)+a_2, P(a_2)+a_3, …, P(a_{q−1})+a_q, P(a_q)+a_1} equal to F(q)."
```lean
def Conj_A229827_i : Prop :=
  ∀ (F : Type) [Field F] [Fintype F] (P : F[X]),
    P.degree < ((Fintype.card F - 1 : ℕ) : WithBot ℕ) → (∀ c : F, P ≠ C c - X) →
    ∃ a : List F, a.Nodup ∧ (∀ x, x ∈ a) ∧
      ∀ y : F, y ∈ List.zipWith (fun x z => P.eval x + z) a (a.rotate 1)
```

**1208.2683 Conj. 3.1** (l. 462–470): "The sequence (ⁿ√F_n)_{n≥2} is strictly increasing, and moreover the sequence (ⁿ⁺¹√F_{n+1}/ⁿ√F_n)_{n≥4} is strictly decreasing. Also, for any integers A > 1 and B ≠ 0 with A² > 4B and (A > 2 or B ≥ −9), the sequence (ⁿ⁺¹√u_{n+1}/ⁿ√u_n)_{n≥4} is strictly decreasing with limit 1, where u_0 = 0, u_1 = 1, and u_{n+1} = Au_n − Bu_{n−1} (n = 1, 2, 3, …)."
```lean
def IsRootFn (K : Type) [Field K] [LinearOrder K] (rt : ℕ → K → K) : Prop :=
  ∀ n : ℕ, 0 < n → ∀ x : K, 0 < x → 0 < rt n x ∧ rt n x ^ n = x
def lucasU (A B : ℤ) : ℕ → ℤ
  | 0 => 0
  | 1 => 1
  | n + 2 => A * lucasU A B (n + 1) - B * lucasU A B n
def Conj_3_1 (K : Type) [Field K] [LinearOrder K] (rt : ℕ → K → K) : Prop :=
  StrictMonoOn (fun n : ℕ => rt n (Nat.fib n : K)) {n | 2 ≤ n} ∧
  StrictAntiOn (fun n : ℕ => rt (n + 1) (Nat.fib (n + 1) : K) / rt n (Nat.fib n : K)) {n | 4 ≤ n} ∧
  ∀ A B : ℤ, 1 < A → B ≠ 0 → 4 * B < A ^ 2 → (2 < A ∨ -9 ≤ B) →
    StrictAntiOn (fun n : ℕ => rt (n + 1) (lucasU A B (n + 1) : K) / rt n (lucasU A B n : K))
        {n | 4 ≤ n} ∧
      ∀ ε : K, 0 < ε → ∃ N : ℕ, ∀ n ≥ N,
        |rt (n + 1) (lucasU A B (n + 1) : K) / rt n (lucasU A B n : K) - 1| < ε
-- not_conj_3_1 (K) [Field K] [LinearOrder K] [IsStrictOrderedRing K] (rt) (hrt : IsRootFn K rt)
```

**1211.1588 Conj. 4.30(iii)** (l. 1339–1340: "p^{(1)}_n = p_n, and p^{(m+1)}_n = p^{(m)}_{p_n}"; l. 1349–1352): "For any m, n ∈ ℤ⁺, we have ⁿ⁺¹√(p^{(m+1)}_{n+1}) / ⁿ√(p^{(m+1)}_n) < ⁿ⁺¹√(p^{(m)}_{n+1}) / ⁿ√(p^{(m)}_n) < 1."
```lean
noncomputable def pr (n : ℕ) : ℕ := Nat.nth Nat.Prime (n - 1)
noncomputable def iterPrime : ℕ → ℕ → ℕ
  | 0, n => n
  | m + 1, n => iterPrime m (pr n)
def Conj_4_30_iii (K : Type) [Field K] [LinearOrder K] (rt : ℕ → K → K) : Prop :=
  ∀ m n : ℕ, 0 < m → 0 < n →
    rt (n + 1) (iterPrime (m + 1) (n + 1) : K) / rt n (iterPrime (m + 1) n : K) <
        rt (n + 1) (iterPrime m (n + 1) : K) / rt n (iterPrime m n : K) ∧
      rt (n + 1) (iterPrime m (n + 1) : K) / rt n (iterPrime m n : K) < 1
```

## Fidelity

- Numberings, permutations, circular orders: duplicate-free lists (`List.Nodup`, `List.Perm`); "all the elements of A": `∀ x, x ∈ l ↔ x ∈ A`; circular neighbours: `l.rotate 1` (last ↦ first); "pairwise distinct": `Nodup` of the list of the n values. Differences in 3.2: `zipWith (· - ·) b b.tail`, the n − 1 values b_i − b_{i+1}.
- Groups: "additive (abelian) group" = `AddCommGroup`, quantified over `G : Type`; "2a" = `2 • a`; "n odd" = `n % 2 = 1`; |G| > 1 = `Nontrivial G` (G may be infinite). Refuting the abelian reading also refutes any wider reading.
- 3.2: "Sylow 2-subgroup is cyclic" is `Syl2Cyclic`, the 2-torsion subgroup (the unique Sylow 2-subgroup of a finite abelian group) has one generator. `not_conj_3_2_with` refutes the statement for every predicate `syl`, since 4 ∤ 6 uses the first disjunct; Mathlib's `Sylow` (279 MB import alone) is not needed.
- Reals: 3.10(ii) is refuted in every field of characteristic 0, A228728 in every linearly ordered field (casts ℤ → R injective; the `decide` runs on ℤ); `RealInst.lean` instantiates both with Mathlib's `ℝ`. "{±s, ±t}" = `(A : Set R) = {s, -s, t, -t}`; "a_1 < … < a_n" = `a.Pairwise (· < ·)`; c_1 = a_1, c_n = a_n = `head?`, `getLast?`.
- p(F) = `pF` via `ringChar` (0 ↦ ⊤); `|S| ≥ min{…}` in `ℕ∞` with `Set.encard` (finite or not). 4.9: (|A| − n) n(n+1)/2 + n(n² − 1)/6 + 1 in ℕ is exact (|A| ≥ n by hypothesis, both quotients integral). MO 316143: λ_k for k = 1, …, n as `lam : ℕ → ℕ` on `Icc 1 n`; gcd = 1 as "1 is the only common divisor"; bound in ℤ, read through `toNat` (a negative bound is vacuous either way). Refuted in F₁₁ (the six sums {3, 4, 5, 7, 8, 9} are the integer sums, as over ℚ); ℚ itself is not used.
- Wang–Sun: p(G) = `pG`, infimum over nonzero g of the least m > 0 with m • g = 0 (`⊤` if none) = "minimum of orders of nonzero elements"; `Conj_WS` is (L) ∧ (C); `not_conj_WS_C` refutes (C) alone. Indices 0-based in `Fin n`.
- A229827(i): `F[X]` is Mathlib's `Polynomial`; degree < q − 1 as `WithBot ℕ`; "not of the form c − x" = `∀ c, P ≠ C c - X`; "equal to F(q)" = every element is a sum. The fields: `F8`, `F9` are `Fin 8`, `Fin 9` with packed tables (F₂[t]/(t³+t+1), F₃[i]/(i²+1)) and a Mathlib `Field` instance whose axioms are checked by `decide +kernel` (`F8Add … F8DistR`); every field with q elements is isomorphic to these, so the instance is a finite field F(q). `GaloisField 2 3` is not used (heavy import, noncomputable). Proof: rotate the order to start at 0 (`CircAux.exists_rotate_head`, `sums_rotate`), then a pruned exhaustive search below 0 (F₈: 1091 nodes; F₉: 5895), split into 57 + 473 `decide +kernel` leaves and reassembled by `match` (`F8Top`, `F9Top`).
- ⁿ√ (3.1, 4.30(iii)): `Mathlib.Analysis.SpecialFunctions.Pow.Real` was killed at 318 MB on import, so ⁿ√ is any `rt` with `IsRootFn K rt` (rt n x > 0, (rt n x)ⁿ = x for n ≥ 1, x > 0) on any linearly ordered field; positive n-th roots are unique, so for K = ℝ every such rt equals x ↦ x^{1/n} on x > 0 and the theorems are the real refutations; `RealInst` states them for ℝ and every such rt (the rpow instance itself is not compiled). "With limit 1" in ε–N form. p_n = `Nat.nth Nat.Prime (n − 1)`; p^{(m)} = `iterPrime m` (p^{(0)}_n = n only as base case; it yields p^{(1)} = p and p^{(m+1)}_n = p^{(m)}_{p_n}).

## Compile

Recipe (`Lean/gen/mk.sh`; all modules in dependency order: `Lean/gen/rebuild.sh`), Lean v4.35.0-rc2, Mathlib v4.35.0-rc2 oleans from the repo `.lake`, one process at a time, `lean -j1`, oleans to `(olean build output, deleted)`:
```sh
LP="$PWD/(olean build output, deleted):$(lake env printenv LEAN_PATH)"     # computed once
LEAN_PATH="$LP" bash Lean/audit/memcap.sh 290 900 /usr/bin/time -l lean -j1 \
  -o (olean build output, deleted) Riemann/Deep/SunPrim/Lean/M.lean
```
Order: Groups ListAux Conj310 A228728 FieldsBase Conj49 Conj49b Conj49c MO316143 WangSun Roots RealInst CircAux F8Defs F8Add F8Mul F8DistL F8DistR F8 F8Search F8Leaf0–3 F8Top F9Defs F9Add F9Mul F9DistL F9DistR F9 F9Search F9Leaf0–17 F9Top A229827. Generated files: `F8*`, `F9*` by `gen/genfield.py`, `gen/gensearch.py F8|F9`; `A229827.lean` by `gen/genconj.py` (run from the repo root).

Output (final clean rebuild of all 52 modules from an empty olean directory, 12:40–12:51Z, `Lean/audit/rebuild.log`): every module exit 0, 0 errors, 0 warnings; no `sorry`, `native_decide`, `skipKernelTC` or `axiom` in the sources (`rg`). `#print axioms` of every theorem below: `[propext, Classical.choice, Quot.sound]`:
`not_conj_3_2`, `not_conj_3_4_i`, `not_conj_A228766`, `not_conj_3_10_ii`, `not_conj_A228728`, `not_conj_4_9_n3`, `not_conj_4_9_n4_F11`, `not_conj_4_9_n4_F13`, `not_conj_MO316143`, `not_conj_WS_C`, `not_conj_WS`, `not_conj_4_30_iii`, `not_conj_3_1`, the four `*_real`, `F8.instField`, `F9.instField`, `F8.no_order_zero`, `F9.no_order_zero`, `not_conj_A229827_i_F8`, `not_conj_A229827_i_F9`.

Peak memory footprint (`/usr/bin/time -l`; the watchdog's `top` samples are coarser), MB: Groups 260, ListAux 221, Conj310 233, A228728 267, FieldsBase 251, Conj49 262, Conj49b 278, Conj49c 284, MO316143 263, WangSun 256, Roots 251, RealInst 240, CircAux 219, F8Defs 213, F8Add 227, F8Mul 233, F8DistL 227, F8DistR 227, F8 214, F8Search 219, F8Leaf0–3 241/237/254/233, F8Top 218, F9Defs 213, F9Add 235, F9Mul 236, F9DistL 235, F9DistR 236, F9 213, F9Search 217, F9Leaf0–17 235–269, F9Top 263, A229827 270; maximum 284. Time per module 7–77 s.

Memory-cap record (290 MB watchdog; `top -l 1` took 5.5 s under load, so samples were ~7 s apart). Killed by the watchdog (not kept): `Mathlib.Analysis.SpecialFunctions.Pow.Real` import 318 MB; a 5040-permutation `decide +kernel` 586 MB; single-file F8 304/349 MB; the first search layouts (one declaration per first entry) 438–702 MB. Above 290 MB without a kill, found afterwards in `/usr/bin/time`: probe H4 374 MB (25 s) and probe K3 296 MB (41 s), both search-layout probes, and the first single `Fields` module 294 MB (split into `FieldsBase`, `Conj49`, `Conj49b`, `Conj49c`, `MO316143`). Imports measured alone: ZMod.Basic 213, Set.Card 191, List.Permutation 140, NormNum.Prime 175, PrimeCounting 209, GroupTheory.Sylow 279, Polynomial Degree+Eval 207.

## Review

Reviewer subagent (statement fidelity; every Lean statement and the source wording pasted inline): 10 rows FAITHFUL; MathOverflow 316143 ISSUES: in `bound316143` the final `+ 1` was inside the binder of the second `∑` (adding ⌊n/2⌋, not 1; equal at n = 3, so the refutation was unaffected). Fixed by parentheses, with `bound316143_check` (the instance's bound is 7; for λ_k = k it equals Conj. 4.9's bound, 11 and 51); re-check: FAITHFUL. Final: 11/11 FAITHFUL. Remarks: the real and root verdicts rely on the ℝ specialisation and `IsRootFn`; the abelian reading suffices for Wang–Sun.

## FORMALIZATIONS.md entry

- Sun's conjectures refuted in Lean 4 + Mathlib (our discovery [sun-counterexamples](discoveries/sun-counterexamples.md)), PARTIAL: 11 of 13 rows (1309.1679 Conj. 3.2, 3.4(i), 3.10(ii); 1811.10503 Conj. 4.9; MO 316143; OEIS A228766, A228728, A229827(i) in F₈, F₉; Wang–Sun (C); 1208.2683 Conj. 3.1; 1211.1588 Conj. 4.30(iii)), each a Prop from the source and its negation (SunLean.not_conj_*), sorry-free, standard axioms, no native_decide; ⁿ√ abstract (Pow.Real exceeds the cap); left: 1208.2683 Conj. 2.1, 1402.6641 Conj. 2.25. Reviewer: FAITHFUL (one bound fixed). [LEAN.md](Riemann/Deep/SunPrim/LEAN.md)
