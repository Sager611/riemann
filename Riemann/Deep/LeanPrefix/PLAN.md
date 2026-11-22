# T4 plan (preliminary; blocked before the first compile, see STATUS.md)

Page: discoveries/prefix-reversals.md (Claim); proofs: ../PrefixRev/NOTE.md; reusable Lean: ../PrefixRev/Lean (copy, namespace Riemann.Deep.LeanPrefix).

## Main statement (sketch)

~~~lean
def rev (n i : ℕ) : Equiv.Perm (Fin n)          -- copy of PrefixRev.rev: x ↦ i-1-x (x < i ≤ n), 0-based
def G (n m k : ℕ) : Subgroup (Equiv.Perm (Fin n)) := Subgroup.closure {rev n n, rev n m, rev n k}
def dihedralPerm (q : ℕ) : DihedralGroup q →* Equiv.Perm (ZMod q)   -- r i • x = i + x, sr i • x = i - x
abbrev Wreath (α q : ℕ) := (ZMod q → DihedralGroup α) ⋊[wreathAct α q] DihedralGroup q  -- D_{2α} ≀ D_{2q}, q-gon action

theorem prefix_reversals (hOff : OffLineLZZW) :
  (∀ n m k, 2 ≤ k → k < m → m < n → m + k < n →                       -- Theorem A
    (G n m k = ⊤ ↔ n % 2 = 1 ∧ m + k = n - 1 ∧ Nat.gcd (m+1) (k+1) = 1 ∧ ¬ (k+1) ∣ m ∧ ¬ (n % 4 = 1 ∧ 4 ∣ k))) ∧
  (∀ n m k, 2 ≤ k → k < m → m + k = n - 1 →                            -- B.1, any parity; g = gcd(m+1, k+1)
    (∀ x y : Fin n, y ∈ MulAction.orbit (G n m k) x ↔
        ((y:ℤ)+1 ≡ (x:ℤ)+1 [ZMOD g] ∨ (y:ℤ)+1 ≡ -((x:ℤ)+1) [ZMOD g])) ∧
    Nat.card (MulAction.orbitRel.Quotient (G n m k) (Fin n)) = g / 2 + 1 ∧
    (MulAction.IsPretransitive (G n m k) (Fin n) ↔ g = 1)) ∧
  (∀ n m k, 2 ≤ k → k < m → m + k = n - 1 → n % 2 = 1 → (k+1) ∣ n →   -- B.3 (α | n forces g = 1)
    ¬ MulAction.IsPreprimitive (G n m k) (Fin n) ∧ Nonempty (G n m k ≃* Wreath (k+1) (n/(k+1)))) ∧
  (∀ n m k, 2 ≤ k → k < m → m + k = n - 1 → n % 2 = 1 → Nat.gcd (m+1) (k+1) = 1 → ¬ (k+1) ∣ n →   -- B.4
    MulAction.IsPreprimitive (G n m k) (Fin n) ∧ alternatingGroup (Fin n) ≤ G n m k ∧
    (G n m k = alternatingGroup (Fin n) ↔ n % 4 = 1 ∧ 4 ∣ k))
~~~
Deviation to argue: the Claim's "when α | n ... otherwise primitive ⊇ A_n" carries Theorem B's hypotheses (n odd; g = 1 for B.4); without them it is false (g > 1: intransitive; n even: two blocks, B.2). Optional extras: B.2; B.3 block sizes (NOTE, not in the Claim).

## Lemma DAG
| # | Lemma | Status |
|---|---|---|
| 1 | Lemma 1, Lemma 2, P n-cycle, parity, interval blocks, Lemma 3 + (★_d), necessity on the line | proved in ../PrefixRev/Lean (copy) |
| 2 | Orbits for any g: P is transitive on each class x+1 ≡ c (mod g) (first return of x ↦ x+α on Z/(n+1)); r_n joins c, -c; orbit count g/2+1 via cls | new |
| 3 | Prop5Core: (★_d), 1 < d < n, d ∣ n ⇒ (k+1) ∣ m (NOTE Steps 1-4, ℕ arithmetic, lap counting) | new, ours |
| 4 | Bridge: PrefixRev.IsPrimitive ↔ MulAction.IsPreprimitive (IsBlock ↔ g•B = B ∨ Disjoint) | new |
| 5 | A_n ≤ G (primitive case), route C or J below | new |
| 6 | G = A_n ↔ n ≡ 1 (4) ∧ 4 ∣ k: sign_rev + alternatingGroup = sign.ker | easy |
| 7 | B.3: interval blocks ⇒ ¬IsPreprimitive; Lemma 6: G = image of Wreath under its imprimitive action on ZMod q × ZMod α transported by x ↦ (x / α, x % α); faithful ⇒ G ≃* Wreath | new (Mathlib has only RegularWreathProduct, order (2α)^{2q}·2q: wrong object) |
| 8 | Off the line: literature input OffLineLZZW = 2610.30892 Thm 3.1, first sentence (as in ../PrefixRev/LEAN.md), or prove NOTE §1 | input |

## A_n ≤ G
- Route J (NOTE): Jones 2014 Thm 1.2, k = 0 (literature input; case (b) needs Mathlib's Projectivization and semilinear maps over a finite field) + Lemma 7 (ours: eigenspaces, Hilbert 90; must be proved). Statement cone grows to 683 MB .olean (memprobe.txt).
- Route C (no CFSG, no Lemma 7): a 3-cycle in G, then Mathlib's Equiv.Perm.alternatingGroup_le_of_isPreprimitive_of_isThreeCycle_mem (GroupAction/Jordan.lean; the prime-cycle version is a Mathlib TODO). If h ∈ G and |h(I) ∩ I| = 1 for I = supp r_k = [1, k], then (r_k · h r_k h⁻¹)² is a 3-cycle. Numerics (n ≤ 301, 4259 primitive line triples): h = P^t works in 3453; the 806 failures (smallest (29,16,12)) need longer words, e.g. h = (r_m r_k)^4 r_m r_n r_m for (33,18,14), r_m r_n r_m r_k r_m r_k r_m r_k r_m r_n r_m r_n r_m for (29,16,12). A uniform construction is open: research before choosing route C.

## Memory plan
Statement file imports Alternating, Primitive, Dihedral, SemidirectProduct (cone 555 MB .olean, ~320 MB footprint est.); proof files Perm.Sign/Closure only (246-281 MB footprint on 12 Nov); one lemma per small file; lean --threads=1, no lake.
