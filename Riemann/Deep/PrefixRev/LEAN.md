PARTIAL: proved sorry-free (standard axioms only): Lemma 1 (P = r_k r_n r_m explicitly and as the first return of x ↦ x + α), Lemma 2 (cyclic coordinates), B.1 (P is an n-cycle for gcd(m+1,k+1) = 1, Mathlib IsCycle; the classes x ≡ ±c mod g are G-invariant; transitive iff gcd = 1), B.2 (n even: two blocks), B.3's interval blocks (α | n), the parity criterion, Lemma 3 (blocks are residue classes of the cyclic coordinate) with (★_d), every necessity direction of Theorem A on the line, and kernel-checked G = S_n for (7,4,2), (11,6,4), (13,10,2) (the literal Conjecture 4 is false). Theorem A and B.4 are proved from four named hypotheses: 2610.30892 Thm 3.1 (off the line, published), Jones 2014 Thm 1.2 case k = 0 (published), Lemma 7 (NOTE, collineation fixed-point bound) and Prop5Core. LEFT: Prop5Core = Steps 1-4 of our Proposition 5 ((★_d) for some divisor 1 < d < n ⇒ (k+1) | m, pure arithmetic, open obligation); not formalized: single orbits on each O_c (exact orbit count), the isomorphism G ≅ D_{2α} ≀ D_{2n/α} and the block sizes of B.3, Lemma 7's proof.

# Lean formalization of the prefix-reversal classification (side result)

Sources: discoveries/prefix-reversals.md, NOTE.md (Lemmas 1-7, Prop 5, Prop 8, Theorems A, B), REPORT.md, NOVELTY.md, verify.py.
Agent: prefixlean (Orca), 2026-11-12, 17:02-19:10 UTC. Lean 4 v4.35.0-rc2, Mathlib v4.35.0-rc2.

## Files (Riemann/Deep/PrefixRev/Lean/, 17 files)

| File | Content | Imports |
|---|---|---|
| Defs | rev n i (prefix reversal on Fin n, 0-based), G n m k = closure {r_n, r_m, r_k}, G_le, fibreFix, fibreStab | Mathlib.GroupTheory.Perm.Sign |
| CycleA | P = r_k r_n r_m, Lemma 1 (P_val, P_firstReturn), modular helpers, phi | Defs, Mathlib.Data.Int.GCD |
| CycleB | Lemma 2: phi_P (shift), phi_rn (j ↦ n+1-j), phi_rk (j ↦ n+2-j on E) | CycleA |
| CycleC | P_transitive, exists_theta, P_cycle (B.1, n-cycle) | CycleB |
| IsCycleP | P.IsCycle, P.support = univ, G_eq_top_of_swap (via Perm.closure_cycle_adjacent_swap) | CycleC, Mathlib.GroupTheory.Perm.Closure |
| Small | coded words, kernel checks (decide +kernel), G_7_4_2, G_11_6_4, G_13_10_2, Conjecture4Literal, not_conjecture4Literal | IsCycleP, Parity |
| Parity | sign_rev (sign r_i = (-1)^(i/2)), G_le_alt_iff, eq_top_of_alt, parity_line | Defs |
| BlocksA | G_intransitive / G_ne_top_of_gcd (gcd ≠ 1), G_le_intervals / G_ne_top_of_dvd (B.3 blocks) | Defs |
| BlocksB | rm_eq, G_ne_top_of_even (B.2), G_le_alt_of_mod4, line_necessary | CycleC, BlocksA, Parity |
| Main | IsBlockOf, IsPrimitive, hypotheses OffLineLZZW, Prop5Primitive, JonesLemma7; line_alt_le, line_B4, theoremA, theoremA_necessary | BlocksB |
| PrimA | block_eq_of_mem, P_pow_eq_iff, P_pow_n, T_add/T_sub/T_mul/T_mod (block translates) | Main |
| PrimB1 | dvd_shift_iff, StarD, block_classes (Lemma 3) | PrimA |
| PrimB2 | star_key, star_of_block ((★_d)), Prop5Core, prop5_of_core, theoremA_core | PrimB1 |
| Jones | nfix, JonesK0 Coll, Lemma7 Coll, line_alt_le_jones, theoremA_jones, line_B4_jones | PrimB2, Mathlib.Data.Nat.Prime.Defs |
| Orbits | cls (orbit label), cls_rev, G_le_cls (B.1 classes invariant), cls_surj | BlocksA, CycleA |
| CheckA, CheckB | #print axioms for 49 theorems | Jones + Orbits; Small |

## Main Lean statements (verbatim; @-entries are Lean's elaborated #check output, which shows section hypotheses)

### Definitions (Defs, CycleA, Main, PrimB1, Jones, Small, Orbits)

~~~lean
/-- The underlying function of the prefix reversal r_i on Fin n. -/
def revFun (n i : ℕ) (x : Fin n) : Fin n :=
  if h : x.val < i ∧ i ≤ n then ⟨i - 1 - x.val, by omega⟩ else x

/-- The prefix reversal r_i as a permutation of Fin n. -/
def rev (n i : ℕ) : Perm (Fin n) where
  toFun := revFun n i
  invFun := revFun n i
  left_inv := revFun_revFun n i
  right_inv := revFun_revFun n i

/-- G(n; m, k) = ⟨r_n, r_m, r_k⟩ ≤ S_n. -/
def G (n m k : ℕ) : Subgroup (Perm (Fin n)) :=
  Subgroup.closure {rev n n, rev n m, rev n k}

/-- P = r_k r_n r_m (r_m acts first). -/
def P (n m k : ℕ) : Perm (Fin n) := rev n k * rev n n * rev n m

/-- Cyclic coordinate (1-based, as in the NOTE): phi x = θ (x+1) mod (n+1). -/
def phi (n θ : ℕ) (x : Fin n) : ℕ := θ * (x.val + 1) % (n + 1)

/-- B is a block of H: every h ∈ H maps B onto B or onto a set disjoint from B. -/
def IsBlockOf {n : ℕ} (H : Subgroup (Perm (Fin n))) (B : Set (Fin n)) : Prop :=
  ∀ h ∈ H, (⇑h '' B = B) ∨ Disjoint (⇑h '' B) B

/-- H is primitive: transitive, and every block is trivial. -/
def IsPrimitive {n : ℕ} (H : Subgroup (Perm (Fin n))) : Prop :=
  (∀ x y : Fin n, ∃ h ∈ H, h x = y) ∧
    ∀ B : Set (Fin n), IsBlockOf H B → B.Subsingleton ∨ B = Set.univ

/-- (★_d) in the 0-based cyclic coordinate psi = phi - 1: each class psi ≡ c (mod d) with
c ≠ 0 lies inside E = {x : x < k} or misses it. (NOTE: classes C_c, c ≢ 1, 1-based.) -/
def StarD (n k θ d : ℕ) : Prop :=
  ∀ x y : Fin n, (phi n θ x - 1) % d = (phi n θ y - 1) % d → (phi n θ x - 1) % d ≠ 0 →
    (x.val < k ↔ y.val < k)

/-- Number of fixed points. -/
def nfix {n : ℕ} (h : Perm (Fin n)) : ℕ := (Finset.univ.filter (fun x => h x = x)).card

/-- Conjecture 4 of Blanco et al. (arXiv 2511.16959v3), literal form:
2 ≤ k < m < n and m + k < n imply that r_n, r_m, r_k do not generate S_n. -/
def Conjecture4Literal : Prop :=
  ∀ n m k : ℕ, 2 ≤ k → k < m → m < n → m + k < n → G n m k ≠ ⊤

/-- The orbit label: c ∈ [0, g/2] with x + 1 ≡ ±c (mod g) (1-based position x + 1). -/
def cls (g : ℕ) {n : ℕ} (x : Fin n) : ℕ :=
  if 2 * ((x.val + 1) % g) ≤ g then (x.val + 1) % g else g - (x.val + 1) % g

~~~

### Lemma 1, Lemma 2, B.1 n-cycle

~~~lean
/-- Lemma 1 (0-based): P x = x + α (x < m), P m = k, P x = x - β (x > m). -/
theorem P_val {n m k : ℕ} (hline : n = m + k + 1) (hkm : k < m) (x : Fin n) :
    (P n m k x).val =
      if x.val < m then x.val + (k + 1) else if x.val = m then k else x.val - (m + 1) := by

/-- Lemma 1, first-return form. Identify x with y = x.val + 1 in Z/(n+1) minus 0 and let
R y = y + α. Then P y = R y if R y ≠ 0, and P y = R (R y) if R y = 0. -/
theorem P_firstReturn {n m k : ℕ} (hline : n = m + k + 1) (hkm : k < m) (x : Fin n) :
    (P n m k x).val + 1 =
      if (x.val + 1 + (k + 1)) % (n + 1) = 0 then (x.val + 1 + 2 * (k + 1)) % (n + 1)
      else (x.val + 1 + (k + 1)) % (n + 1) := by

@PrefixRev.phi_P : ∀ {n m k θ : ℕ},
  n = m + k + 1 →
    k < m →
      (k + 1) * θ % (n + 1) = 1 →
        ∀ (x : Fin n),
          PrefixRev.phi n θ ((PrefixRev.P n m k) x) = if PrefixRev.phi n θ x < n then PrefixRev.phi n θ x + 1 else 1

@PrefixRev.phi_rn : ∀ {n m k θ : ℕ},
  n = m + k + 1 →
    k < m →
      (k + 1) * θ % (n + 1) = 1 → ∀ (x : Fin n), PrefixRev.phi n θ ((PrefixRev.rev n n) x) = n + 1 - PrefixRev.phi n θ x

@PrefixRev.phi_rk : ∀ {n m k θ : ℕ},
  n = m + k + 1 →
    k < m →
      (k + 1) * θ % (n + 1) = 1 →
        ∀ (x : Fin n),
          PrefixRev.phi n θ ((PrefixRev.rev n k) x) =
            if ↑x < k then n + 2 - PrefixRev.phi n θ x else PrefixRev.phi n θ x

/-- Theorem B.1 (n-cycle): on the line n = m + k + 1 with gcd(m+1, k+1) = 1, the powers of
P = r_k r_n r_m move every point to every point, so P is an n-cycle. -/
theorem P_cycle {n m k : ℕ} (hline : n = m + k + 1) (hkm : k < m)
    (hg : Nat.Coprime (k + 1) (m + 1)) (x y : Fin n) : ∃ t : ℕ, (P n m k ^ t) x = y := by

/-- Theorem B.1 (n-cycle, Mathlib form). -/
theorem P_isCycle {n m k : ℕ} (hline : n = m + k + 1) (hkm : k < m)
    (hg : Nat.Coprime (k + 1) (m + 1)) : (P n m k).IsCycle := by

theorem P_support {n m k : ℕ} (hline : n = m + k + 1) (hkm : k < m) :
    (P n m k).support = Finset.univ :=

/-- Criterion: if G contains the transposition (x, P x) for some x, then G = S_n. -/
theorem G_eq_top_of_swap {n m k : ℕ} (hline : n = m + k + 1) (hkm : k < m)
    (hg : Nat.Coprime (k + 1) (m + 1)) (x : Fin n) (hs : swap x (P n m k x) ∈ G n m k) :
    G n m k = ⊤ := by

~~~

### Parity

~~~lean
/-- sign r_i = (-1)^(i/2) for i ≤ n. -/
theorem sign_rev (n i : ℕ) (h : i ≤ n) : Perm.sign (rev n i) = (-1) ^ (i / 2) := by

/-- Parity: G ≤ A_n iff r_n, r_m, r_k are all even iff n/2, m/2, k/2 are all even. -/
theorem G_le_alt_iff {n m k : ℕ} (hm : m ≤ n) (hk : k ≤ n) :
    (∀ g ∈ G n m k, Perm.sign g = 1) ↔ (n / 2 % 2 = 0 ∧ m / 2 % 2 = 0 ∧ k / 2 % 2 = 0) := by

/-- A subgroup containing A_n and an odd permutation is S_n. -/
theorem eq_top_of_alt {n : ℕ} {H : Subgroup (Perm (Fin n))}
    (hA : ∀ σ : Perm (Fin n), Perm.sign σ = 1 → σ ∈ H) (g : Perm (Fin n)) (hg : g ∈ H)
    (hodd : Perm.sign g = -1) : H = ⊤ := by

~~~

### Necessity on the line (B.1 intransitivity, B.2, B.3 blocks, parity)

~~~lean
/-- g = gcd(k+1, m+1) ≠ 1 ⇒ G ≠ S_n. -/
theorem G_ne_top_of_gcd {n m k : ℕ} (hline : n = m + k + 1) (hkm : k < m)
    (hg : Nat.gcd (k + 1) (m + 1) ≠ 1) : G n m k ≠ ⊤ := by

/-- Lemma 6 / B.3: if α = k + 1 divides n (equivalently divides m) then the intervals of
length α are blocks of G, and G ≠ S_n. -/
theorem G_le_intervals {n m k : ℕ} (hline : n = m + k + 1) (hkm : k < m) (hdiv : (k + 1) ∣ n) :
    G n m k ≤ fibreStab (fun x : Fin n => x.val / (k + 1)) := by

theorem G_ne_top_of_dvd {n m k : ℕ} (hk : 2 ≤ k) (hline : n = m + k + 1) (hkm : k < m)
    (hdiv : (k + 1) ∣ m) : G n m k ≠ ⊤ := by

/-- B.2: n even on the line with gcd = 1 ⇒ G ≠ S_n. -/
theorem G_ne_top_of_even {n m k θ : ℕ} (hk2 : 2 ≤ k) (hline : n = m + k + 1) (hkm : k < m)
    (hθ : (k + 1) * θ % (n + 1) = 1) (hev : n % 2 = 0) : G n m k ≠ ⊤ := by

/-- n ≡ 1 (mod 4), 4 | k on the line ⇒ G ≤ A_n, so G ≠ S_n. -/
theorem G_le_alt_of_mod4 {n m k : ℕ} (hline : n = m + k + 1) (_hkm : k < m)
    (h4 : n % 4 = 1 ∧ k % 4 = 0) : ∀ g ∈ G n m k, Perm.sign g = 1 :=

/-- Theorem A, necessity on the line n = m + k + 1 (2 ≤ k < m):
G = S_n ⇒ n odd, gcd(m+1, k+1) = 1, (k+1) ∤ m, and not (n ≡ 1 mod 4 and 4 | k). -/
theorem line_necessary {n m k : ℕ} (hk : 2 ≤ k) (hkm : k < m) (hline : n = m + k + 1)
    (htop : G n m k = ⊤) :
    n % 2 = 1 ∧ Nat.Coprime (m + 1) (k + 1) ∧ ¬ (k + 1) ∣ m ∧ ¬ (n % 4 = 1 ∧ k % 4 = 0) := by

/-- Theorem A, the direction that needs no Prop 5 and no Jones: G = S_n ⇒ the conditions
(only the published off-line result is assumed). -/
theorem theoremA_necessary (hOff : OffLineLZZW) {n m k : ℕ} (hk : 2 ≤ k) (hkm : k < m)
    (hmn : m < n) (hsum : m + k < n) (htop : G n m k = ⊤) :
    n % 2 = 1 ∧ m + k = n - 1 ∧ Nat.Coprime (m + 1) (k + 1) ∧ ¬ (k + 1) ∣ m ∧
      ¬ (n % 4 = 1 ∧ k % 4 = 0) := by

@PrefixRev.G_le_cls : ∀ {n m k : ℕ},
  n = m + k + 1 → k < m → PrefixRev.G n m k ≤ PrefixRev.fibreFix fun x => PrefixRev.cls ((k + 1).gcd (m + 1)) x

/-- Every label 0 ≤ c ≤ g/2 occurs. -/
theorem cls_surj {n m k : ℕ} (hline : n = m + k + 1) (c : ℕ)
    (hc : 2 * c ≤ Nat.gcd (k + 1) (m + 1)) :
    ∃ x : Fin n, cls (Nat.gcd (k + 1) (m + 1)) x = c := by

~~~

### Counterexamples (kernel-checked)

~~~lean
/-- (7, 4, 2): r_7, r_4, r_2 generate S_7 although 4 + 2 < 7. -/
theorem G_7_4_2 : G 7 4 2 = ⊤ := by

/-- (11, 6, 4): r_11, r_6, r_4 generate S_11 although 6 + 4 < 11. -/
theorem G_11_6_4 : G 11 6 4 = ⊤ := by

/-- (13, 10, 2): r_13, r_10, r_2 generate S_13 although 10 + 2 < 13. -/
theorem G_13_10_2 : G 13 10 2 = ⊤ := by

/-- The literal Conjecture 4 is false: (7, 4, 2). -/
theorem not_conjecture4Literal : ¬ Conjecture4Literal := fun h =>

~~~

### Lemma 3 and (★_d)

~~~lean
@PrefixRev.block_classes : ∀ {n m k θ : ℕ},
  n = m + k + 1 →
    k < m →
      (k + 1) * θ % (n + 1) = 1 →
        ∀ {B : Set (Fin n)},
          PrefixRev.IsBlockOf (PrefixRev.G n m k) B →
            ∀ {x0 : Fin n},
              x0 ∈ B →
                B.Subsingleton ∨
                  B = Set.univ ∨
                    ∃ d,
                      1 < d ∧
                        d < n ∧
                          d ∣ n ∧ ∀ (y : Fin n), y ∈ B ↔ (PrefixRev.phi n θ y - 1) % d = (PrefixRev.phi n θ x0 - 1) % d

@PrefixRev.star_of_block : ∀ {n m k θ : ℕ},
  n = m + k + 1 →
    k < m →
      (k + 1) * θ % (n + 1) = 1 →
        n % 2 = 1 →
          ∀ {B : Set (Fin n)},
            PrefixRev.IsBlockOf (PrefixRev.G n m k) B →
              ∀ {x0 : Fin n} {d : ℕ},
                1 < d →
                  d < n →
                    d ∣ n →
                      (∀ (y : Fin n), y ∈ B ↔ (PrefixRev.phi n θ y - 1) % d = (PrefixRev.phi n θ x0 - 1) % d) →
                        PrefixRev.StarD n k θ d

/-- Lemma 3 + (★_d): Prop5Core implies the primitivity statement Prop5Primitive. -/
theorem prop5_of_core (hC : Prop5Core) : Prop5Primitive := by

~~~

### Theorem A and B.4

~~~lean
/-- Theorem A from: the published off-line result, the arithmetic core of Prop 5 (open),
Jones's theorem (k = 0) and Lemma 7. -/
theorem theoremA_jones {Coll : ∀ {n : ℕ}, ℕ → ℕ → Subgroup (Perm (Fin n)) → Prop}
    (hOff : OffLineLZZW) (hC : Prop5Core) (hJ : JonesK0 Coll) (hL : Lemma7 Coll) {n m k : ℕ}
    (hk : 2 ≤ k) (hkm : k < m) (hmn : m < n) (hsum : m + k < n) :
    G n m k = ⊤ ↔ (n % 2 = 1 ∧ m + k = n - 1 ∧ Nat.Coprime (m + 1) (k + 1) ∧ ¬ (k + 1) ∣ m ∧
      ¬ (n % 4 = 1 ∧ k % 4 = 0)) := by

/-- B.4, the alternating branch, from the same hypotheses: if n ≡ 1 (mod 4) and 4 | k then
G = A_n (the S_n branch is theoremA_jones). -/
theorem line_B4_jones {Coll : ∀ {n : ℕ}, ℕ → ℕ → Subgroup (Perm (Fin n)) → Prop}
    (hC : Prop5Core) (hJ : JonesK0 Coll) (hL : Lemma7 Coll) {n m k : ℕ} (hk : 2 ≤ k)
    (hkm : k < m) (hline : n = m + k + 1) (hodd : n % 2 = 1) (hg : Nat.Coprime (m + 1) (k + 1))
    (hdiv : ¬ (k + 1) ∣ m) (h4 : n % 4 = 1 ∧ k % 4 = 0) :
    ∀ σ : Perm (Fin n), σ ∈ G n m k ↔ Perm.sign σ = 1 := fun σ =>

theorem theoremA (hOff : OffLineLZZW) (hP5 : Prop5Primitive) (hJ : JonesLemma7) {n m k : ℕ}
    (hk : 2 ≤ k) (hkm : k < m) (hmn : m < n) (hsum : m + k < n) :
    G n m k = ⊤ ↔ (n % 2 = 1 ∧ m + k = n - 1 ∧ Nat.Coprime (m + 1) (k + 1) ∧ ¬ (k + 1) ∣ m ∧
      ¬ (n % 4 = 1 ∧ k % 4 = 0)) := by

/-- B.4: on the odd line with gcd = 1 and (k+1) ∤ m: G = A_n if n ≡ 1 (mod 4) and 4 | k,
and G = S_n otherwise. -/
theorem line_B4 (hP5 : Prop5Primitive) (hJ : JonesLemma7) {n m k : ℕ} (hk : 2 ≤ k)
    (hkm : k < m) (hline : n = m + k + 1) (hodd : n % 2 = 1) (hg : Nat.Coprime (m + 1) (k + 1))
    (hdiv : ¬ (k + 1) ∣ m) :
    ((n % 4 = 1 ∧ k % 4 = 0) → ∀ σ : Perm (Fin n), σ ∈ G n m k ↔ Perm.sign σ = 1) ∧
    (¬ (n % 4 = 1 ∧ k % 4 = 0) → G n m k = ⊤) := by

/-- Theorem A with Prop 5 reduced to its arithmetic core. -/
theorem theoremA_core (hOff : OffLineLZZW) (hC : Prop5Core) (hJ : JonesLemma7) {n m k : ℕ}
    (hk : 2 ≤ k) (hkm : k < m) (hmn : m < n) (hsum : m + k < n) :
    G n m k = ⊤ ↔ (n % 2 = 1 ∧ m + k = n - 1 ∧ Nat.Coprime (m + 1) (k + 1) ∧ ¬ (k + 1) ∣ m ∧
      ¬ (n % 4 = 1 ∧ k % 4 = 0)) :=

~~~

## Named hypotheses (Props passed as arguments; no axiom declarations)

### OffLineLZZW

Source: Lei, Zhao, Zhu, Wu, arXiv 2610.30892v1, Theorem 3.1, first sentence (TeX l. 221-223: "If m+k<n-1, then <r_n,r_m,r_k> < S_n"), standing assumption 2 <= k < m < n (TeX l. 109-111). Published, not formalized here (its graph/edge-count proof is NOTE section 1). Lean's natural n - 1 is exact (n >= 4).

~~~lean
/-- Named hypothesis: 2610.30892, Theorem 3.1, first sentence. -/
def OffLineLZZW : Prop :=
  ∀ n m k : ℕ, 2 ≤ k → k < m → m < n → m + k < n - 1 → G n m k ≠ ⊤
~~~

### Prop5Core

Source: OUR Proposition 5, Steps 1-4 (NOTE section 4), restricted to the odd coprime line. OPEN OBLIGATION (new mathematics, not formalized; verify.py checks it for n <= 101). StarD is the NOTE's (★_d) in psi = phi - 1 coordinates.

~~~lean
def Prop5Core : Prop :=
  ∀ n m k θ d : ℕ, 2 ≤ k → k < m → n = m + k + 1 → n % 2 = 1 → Nat.Coprime (m + 1) (k + 1) →
    (k + 1) * θ % (n + 1) = 1 → 1 < d → d < n → d ∣ n → StarD n k θ d → (k + 1) ∣ m
~~~

### JonesK0

Source: G. A. Jones, Primitive permutation groups containing a cycle, Bull. Aust. Math. Soc. 89 (2014) 159-165, Theorem 1.2, case k = 0 (arXiv 1209.5169, TeX l. 54-72; rests on Jones, J. Group Theory 5 (2002) and CFSG). Stated with weakenings that Jones implies: (a) C_p <= G <= AGL_1(p) replaced by its consequence nfix <= 1 for nonidentity elements; (b) through the abstract predicate Coll q d H, intended meaning PGL_d(q) <= H <= PΓL_d(q) on the points of PG(d-1,q) (Mathlib: Projectivization; not instantiated: too heavy for the 290 MB cap); (c) G = L_2(11), M_11, M_23 replaced by the fixed-point bounds 3 (degree 11, M_11 sharply 4-transitive) and 7 (M_23 on 23 points).

~~~lean
/-- Named hypothesis: Jones 2014, Theorem 1.2, case k = 0 (see the module doc). -/
def JonesK0 (Coll : ∀ {n : ℕ}, ℕ → ℕ → Subgroup (Perm (Fin n)) → Prop) : Prop :=
  ∀ (n : ℕ) (H : Subgroup (Perm (Fin n))), IsPrimitive H →
    ¬ (∀ σ : Perm (Fin n), Perm.sign σ = 1 → σ ∈ H) →
    (∃ c ∈ H, ∀ x y : Fin n, ∃ t : ℕ, (c ^ t) x = y) →
    (n.Prime ∧ ∀ h ∈ H, h ≠ 1 → nfix h ≤ 1) ∨
    (∃ q d : ℕ, (∃ p e : ℕ, p.Prime ∧ 0 < e ∧ q = p ^ e) ∧ 2 ≤ d ∧ n * (q - 1) = q ^ d - 1 ∧
      Coll q d H) ∨
    (n = 11 ∧ ∀ h ∈ H, h ≠ 1 → nfix h ≤ 3) ∨ (n = 23 ∧ ∀ h ∈ H, h ≠ 1 → nfix h ≤ 7)
~~~

### Lemma7

Source: NOTE Lemma 7 (fixed points of collineations; proof in NOTE section 5, brute-force check in verify.py). Named, not proved here.

~~~lean
/-- Named hypothesis: the NOTE's Lemma 7 (a nonidentity collineation of PG(d-1, q),
(d, q) ≠ (2, 4), fixes at most n/2 points). -/
def Lemma7 (Coll : ∀ {n : ℕ}, ℕ → ℕ → Subgroup (Perm (Fin n)) → Prop) : Prop :=
  ∀ (n q d : ℕ) (H : Subgroup (Perm (Fin n))), Coll q d H → 2 ≤ d → n * (q - 1) = q ^ d - 1 →
    ¬ (d = 2 ∧ q = 4) → ∀ h ∈ H, h ≠ 1 → 2 * nfix h ≤ n
~~~

### Prop5Primitive

Used by Main.theoremA (older interface): primitivity on the odd line. Proved from Prop5Core by PrimB2.prop5_of_core (Lemma 3 + (★_d)).

~~~lean
/-- Named hypothesis (OPEN OBLIGATION, our Proposition 5): primitivity on the odd line. -/
def Prop5Primitive : Prop :=
  ∀ n m k : ℕ, 2 ≤ k → k < m → n = m + k + 1 → n % 2 = 1 → Nat.Coprime (m + 1) (k + 1) →
    ¬ (k + 1) ∣ m → IsPrimitive (G n m k)
~~~

### JonesLemma7

Used by Main.theoremA (older interface): the merged consequence of JonesK0 + Lemma7. theoremA_jones (Jones.lean) uses JonesK0 and Lemma7 separately instead; prefer it.

~~~lean
/-- Named hypothesis: Jones 2014 Thm 1.2 (k = 0) + Lemma 7 (+ fixed points in the
sporadic and affine cases). -/
def JonesLemma7 : Prop :=
  ∀ (n : ℕ) (H : Subgroup (Perm (Fin n))), IsPrimitive H →
    (∃ c ∈ H, ∀ x y : Fin n, ∃ t : ℕ, (c ^ t) x = y) →
    (∃ h ∈ H, h ≠ 1 ∧ n < 2 * (Finset.univ.filter (fun x => h x = x)).card) →
    ∀ σ : Perm (Fin n), Perm.sign σ = 1 → σ ∈ H
~~~

Main entry point: theoremA_jones (hOff : OffLineLZZW) (hC : Prop5Core) (hJ : JonesK0 Coll) (hL : Lemma7 Coll). theoremA_necessary needs only OffLineLZZW; line_necessary needs nothing.

## Build

Clean rebuild 2026-11-12 18:18-18:57 UTC (Riemann/Deep/PrefixRev/Lean/audit/verify.sh; Main's 1650 s wall time is a host stall, 4.5 s user). One module at a time, no lake process, oleans to .scratch/prefixlean/olean (build output, deleted at the end):

    cp Riemann/Deep/DPCUncond/Lean/audit/leanpath.txt Riemann/Deep/PrefixRev/Lean/audit/leanpath.txt   # lake's LEAN_PATH
    # for m in Defs CycleA CycleB CycleC Parity IsCycleP Small BlocksA BlocksB Main PrimA PrimB1 PrimB2 Jones Orbits CheckA CheckB:
    LEAN_PATH=.scratch/prefixlean/olean:$(cat Riemann/Deep/PrefixRev/Lean/audit/leanpath.txt) bash .scratch/memcap.sh 290 1800 \
      /usr/bin/time -l ~/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean -j1 \
      -o .scratch/prefixlean/olean/Riemann/Deep/PrefixRev/Lean/$m.olean Riemann/Deep/PrefixRev/Lean/$m.lean

| Module | rc | peak footprint (time -l) | memcap peak | wall |
|---|---|---|---|---|
| Defs | 0 | 246 MB | 217 MB | 27.58 s |
| CycleA | 0 | 275 MB | 230 MB | 60.16 s |
| CycleB | 0 | 263 MB | 180 MB | 47.67 s |
| CycleC | 0 | 253 MB | 232 MB | 8.92 s |
| Parity | 0 | 279 MB | 231 MB | 177.34 s |
| IsCycleP | 0 | 267 MB | 243 MB | 7.85 s |
| Small | 0 | 281 MB | 252 MB | 62.03 s |
| BlocksA | 0 | 261 MB | 181 MB | 18.54 s |
| BlocksB | 0 | 272 MB | 235 MB | 9.79 s |
| Main | 0 | 269 MB | 244 MB | 1650.64 s |
| PrimA | 0 | 252 MB | 199 MB | 34.49 s |
| PrimB1 | 0 | 253 MB | 222 MB | 24.28 s |
| PrimB2 | 0 | 270 MB | 199 MB | 33.83 s |
| Jones | 0 | 274 MB | 253 MB | 28.26 s |
| Orbits | 0 | 256 MB | 230 MB | 28.37 s |
| CheckA | 0 | 241 MB | 183 MB | 25.76 s |
| CheckB | 0 | 264 MB | 231 MB | 28.47 s |

All 17 modules rc = 0; maximum footprint 281 MB (Small). Kernel checks: decide +kernel on Bool evaluations (Small.lean).

## Static scan

    rg -n -w 'sorry|admit|native_decide|skipKernelTC' Riemann/Deep/PrefixRev/Lean   # no matches
    rg -n '^\s*axiom\s' Riemann/Deep/PrefixRev/Lean                              # no matches

## #print axioms

CheckA.lean (40) and CheckB.lean (9), clean rebuild:

~~~
'PrefixRev.rev_val' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.rev_rev' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.P_val' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.P_firstReturn' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.phi_ne_zero' depends on axioms: [propext, Quot.sound]
'PrefixRev.phi_inj' depends on axioms: [propext, Quot.sound]
'PrefixRev.phi_P' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.phi_rn' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.phi_rk' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.P_transitive' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.exists_theta' depends on axioms: [propext, Quot.sound]
'PrefixRev.P_cycle' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.sign_rev' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.G_le_alt_iff' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.eq_top_of_alt' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.parity_line' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.G_intransitive' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.G_ne_top_of_gcd' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.G_le_intervals' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.G_ne_top_of_dvd' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.rm_eq' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.G_ne_top_of_even' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.G_le_alt_of_mod4' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.G_ne_top_of_mod4' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.line_necessary' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.theoremA_necessary' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.fix_rk_card' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.line_alt_le' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.line_B4' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.theoremA' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.block_classes' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.star_of_block' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.prop5_of_core' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.theoremA_core' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.line_alt_le_jones' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.theoremA_jones' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.line_B4_jones' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.cls_rev' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.G_le_cls' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.cls_surj' depends on axioms: [propext, Quot.sound]
'PrefixRev.P_mem' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.P_support' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.P_isCycle' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.G_eq_top_of_swap' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.codePerm_eq_swap' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.G_7_4_2' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.G_11_6_4' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.G_13_10_2' depends on axioms: [propext, Classical.choice, Quot.sound]
'PrefixRev.not_conjecture4Literal' depends on axioms: [propext, Classical.choice, Quot.sound]
~~~

## Fidelity

Definitions against the sources (2511.16959v3 = DMTCS 28:1 TeX; 2610.30892v1 TeX; NOTE.md):

- rev n i (Defs): r_i of 2511.16959 TeX l. 135-142 ("reverse the order of the first i elements", positions x ↦ i+1-x for x ≤ i, 1 < i ≤ n). Lean positions are 0-based (Lean x = NOTE x - 1), so x ↦ i-1-x for x < i. For i > n Lean's rev is the identity; never used (all theorems have i ∈ {n, m, k} with k < m < n). As a permutation of positions; the subgroup generated does not depend on the composition convention.
- G n m k = Subgroup.closure {r_n, r_m, r_k} = ⟨r_n, r_m, r_k⟩ (2610.30892 TeX l. 109-111; NOTE Notation). "= S_n" is = ⊤, "< S_n" is ≠ ⊤. A_n ≤ H is written ∀ σ, sign σ = 1 → σ ∈ H (Mathlib's alternatingGroup is the kernel of sign); G = A_n is ∀ σ, σ ∈ G ↔ sign σ = 1.
- P = rev n k * rev n n * rev n m: (σ * τ) x = σ (τ x), so r_m acts first, as the NOTE's r_k r_n r_m (products right to left).
- phi n θ x = θ (x+1) mod (n+1) is the NOTE's φ at the 1-based position x+1 (values 1..n); psi = phi - 1 is 0-based; "x.val < k" is the NOTE's x ∈ [1, k] = φ⁻¹(E).
- Conjecture4Literal: Conjecture 4, first sentence (TeX l. 884-886) under the standing 2 ≤ k < m < n (TeX l. 872). Its "Moreover" sentence (parity-split form) is not encoded.
- OffLineLZZW: Thm 3.1 of 2610.30892, first sentence (TeX l. 221-223). Its second sentence (n even, m+k < n) is not assumed: on the line it is proved here (G_ne_top_of_even), off the line it is covered by the first sentence.
- IsBlockOf/IsPrimitive (Main): B is a block iff every h maps B onto B or onto a set disjoint from B; primitive = transitive and every block subsingleton or univ. This is Mathlib's MulAction.IsPreprimitive (IsPretransitive + every IsBlock is IsTrivialBlock; for groups, IsBlock ⟺ the "g•B = B or disjoint" form). Not bridged in Lean (the Primitive import was not attempted under the cap).
- fibreFix b: permutations fixing every fibre of b setwise; fibreStab b: permutations preserving the partition into fibres (fibres may be exchanged). Reviewer round 2 confirmed both readings.
- nfix h = card {x | h x = x}; "2 * nfix h ≤ n" is nfix(h) ≤ n/2 exactly.

Theorems against NOTE.md:

- Lemma 1: P_val (0-based form of P(x) = x+α, x ≤ m; P(m+1) = α; P(x) = x-β, x ≥ m+2) and P_firstReturn (P y = R y if R y ≠ 0, else R(R y), R y = y + α on Z/(n+1), y = x+1). Exact.
- Lemma 2: phi_P, phi_rn, phi_rk (all with n = m+k+1, k < m and (k+1)θ ≡ 1 mod n+1, which forces gcd = 1). Exact. The Beatty description E = {⌈sN/α⌉} is not formalized.
- B.1: P_cycle (∀ x y, ∃ t, P^t x = y) and P_isCycle + P_support (Mathlib IsCycle, support = univ) for gcd(m+1, k+1) = 1; G_le_cls + cls_surj: every class O_c = {x : x ≡ ±c mod g}, 0 ≤ c ≤ g/2, is G-invariant and nonempty; G_ne_top_of_gcd / G_intransitive: g ≠ 1 ⇒ intransitive. So "transitive iff g = 1" is proved; "each O_c is one orbit" (hence exactly ⌊g/2⌋+1 orbits) is proved only for g = 1.
- B.2: G_ne_top_of_even (blocks = parity classes of phi). B.3: G_le_intervals (the intervals of length α are blocks when α | n), G_ne_top_of_dvd. The isomorphism G ≅ D_{2α} ≀ D_{2q}, |G| and the block sizes are not formalized.
- Lemma 3: block_classes (a nonempty block is a singleton, everything, or a residue class of psi mod d, d | n, 1 < d < n) and star_of_block (then (★_d) holds; n odd). The converse ((★_d) ⇒ 𝒞_d invariant) is not needed and not formalized.
- Prop 5: Prop5Core is its arithmetic content (Steps 1-4) restricted to the odd coprime line (reviewer: faithful restricted form); prop5_of_core turns it into primitivity. OPEN.
- Lemma 7: named (Lemma7 Coll). Prop 8: line_alt_le_jones (A_n ⊆ G) from Prop5Core, JonesK0, Lemma7: r_k ≠ 1 fixes n - k = m+1 > n/2 points (fix_rk_card); (d, q) = (2, 4) is excluded because it forces n = 5 < 7.
- Parity: sign_rev (sign r_i = (-1)^(i/2)), G_le_alt_iff, parity_line. B.4: line_B4 / line_B4_jones (A_n branch) and theoremA_jones (S_n branch).
- Theorem A: theoremA_jones (both directions; hypotheses 2 ≤ k < m < n, m+k < n; conditions n % 2 = 1, m+k = n-1, Nat.Coprime (m+1) (k+1), ¬ (k+1) ∣ m, ¬ (n % 4 = 1 ∧ k % 4 = 0)). theoremA_necessary: the "only if" direction from OffLineLZZW alone; line_necessary: on the line, from nothing.
- Counterexamples: G_7_4_2, G_11_6_4, G_13_10_2 (kernel evaluation of one word equal to the transposition (0, P 0) on all points, then Perm.closure_cycle_adjacent_swap); not_conjecture4Literal.

## Memory cap and incidents

Every compile ran alone, single-threaded (lean -j1), under bash .scratch/memcap.sh 290 (kills the tree above 290 MB, samples about every 3 s) wrapping /usr/bin/time -l (exact peak footprint). The import floor is the problem: Mathlib.GroupTheory.Perm.Sign alone is 232 MB, Perm.Closure 255 MB, Data.ZMod.Basic 251 MB (5-line test files), so files were kept small and ZMod was avoided (all modular arithmetic is on ℕ).

- A first one-file Cycle.lean reached 287.6 MB (time -l): split into CycleA/B/C.
- A first PrimB.lean was killed by memcap at 303 MB: split into PrimB1/PrimB2, the key step isolated in small arithmetic lemmas.
- Small23.lean (kernel words for the n = 23 triples, case (c) of Jones): the first version was killed by memcap at 402 MB; a chunked version finished (rc 0) but /usr/bin/time reported a 466 MB peak footprint that the 3-s sampling missed, i.e. a short cap breach. The file was deleted and is not part of the deliverable; case (c) is handled through fixed-point bounds in JonesK0 instead. The deliverable's compiles all stay below 290 MB (table below).
- lean's "maximum resident set size" (1.0-1.7 GB) counts mapped olean pages; the cap's measure is top's MEM = the physical footprint.

## Review

reviewer subagent on statement fidelity (definitions, main statements and the published statements pasted inline). Round 1: ISSUES; blocking only as an evidence gap (the extracted excerpts hid the section hypotheses of phi_P, phi_rn, phi_rk, block_classes, star_of_block), minor: line_B4_jones's docstring said "exactly when" but it proves only the A_n branch; supply fibreFix/fibreStab. Everything else faithful: rev/G, P, Conjecture4Literal, OffLineLZZW (n-1 harmless), JonesK0 (all three weakenings implied by Jones; the bounds 3 and 7 hold for all nonidentity elements of L_2(11), M_11, M_23), Lemma7, StarD/Prop5Core (faithful restricted form), IsPrimitive, theoremA*; no vacuity found. Fix: docstring rewritten; Lean's elaborated #check signatures and #print definitions supplied. Round 2: VERDICT: FAITHFUL (Lemma 2 signatures carry the line, order and inverse hypotheses; Lemma 3 is stated for every nonempty block without the unneeded converse; star_of_block assumes n odd; fibreFix is a setwise fibre fixer, fibreStab a partition stabilizer).

## Not done

- Prop5Core (Steps 1-4 of Proposition 5: lap counting, Lemma 4's classes in x-coordinates). This is the remaining open obligation of the new result.
- Lemma 7's proof (eigenspaces, Hilbert 90) and an instance of Coll via Mathlib's Projectivization (import too heavy for the cap).
- Exact orbits (transitivity on each O_c for g > 1), the wreath-product isomorphism, |G| and the block sizes of B.3.
- A bridge from IsPrimitive to Mathlib's MulAction.IsPreprimitive.

The repo-root .lake (a real Mathlib cache directory, 7.4 GB) was left as it is: no lake build, no cache get, no placeholder restored (the orchestrator restores it).

## FORMALIZATIONS.md entry

Prefix reversals, line m + k = n − 1 (side result): Lean, sorry-free — Lemmas 1–3, B.1 n-cycle, B.2, B.3 blocks, parity, Theorem A's necessity, kernel-checked G = S_n for (7,4,2), (11,6,4), (13,10,2) (literal Conjecture 4 false); Theorem A/B.4 from 2610.30892 Thm 3.1, Jones Thm 1.2, Lemma 7 and open Prop5Core. PARTIAL. [LEAN.md](Riemann/Deep/PrefixRev/LEAN.md)
