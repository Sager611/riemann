FORMALIZED: Theorem 1 (the exact set of perfect N for each of the 24 orders of + × − ÷, start {2}) and the Corollary (an order has perfect sequences for all large N iff − is third or fourth and it is not ÷ × − +), from no named hypotheses: sorry-free, axioms propext, Classical.choice, Quot.sound only; the computer-assisted parts (132 patches, 381 base sequences, the small-N search, the 13 dead ends) are Bool checkers proved sound and evaluated by the kernel (decide +kernel, no native_decide).

# Lean formalization: which cyclic orders of + × − ÷ make every number

Answer to the open question of arXiv 2610.10095 (v2 l. 702; v1 l. 922–925), as in [NOTE.md](NOTE.md). Core Lean 4 (v4.35.0-rc2), no Mathlib, no lake build. Directory Lean/:

| file | content |
|---|---|
| Defs.lean | rules, `Valid`, `Perfect`; Bool checker `check`, `check_iff` (check = true ↔ Perfect) |
| Ext.lean | `Steps` (operands made, right operation); `perfect_extend`: perfect for P + moves whose results are exactly (P, N] ⇒ perfect for N |
| Bulk.lean | bulk rounds j = 3..J (NOTE.md Lemma 2): `bulk_steps`, `bulk_vals` (makes exactly X ∪ Y), `bulk_nodup` |
| Patch.lean | values a·q + b, patch checker `patchOK`; Lemma 3 `patch_sound` (P = 4q + r, q ≥ q* ⇒ 3P + c) |
| Induct.lean | Theorem 2 scheme `works_of` (strong induction on N = 3P + c), `base_of` |
| D_<order>.lean ×11 | generated: σ, its 12 patches, base sequences N0 ≤ N < N1, kernel checks, `works` (all N ≥ N0) |
| Search.lean | exhaustive search over reachable collections: `reach`, `noPerfect`; `reach_sound`, `noPerfect_sound`, `bound_of_reach` |
| Fail.lean | generated: the 13 failing orders, perfect N ≤ 5 (no legal 1-, 2- or 5-move prefix) |
| Small1–6.lean | generated: the exact perfect N below N0 (working) or below 6 (failing) |
| Main.lean | generated: Theorem 1 (`thm1_<order>`, 24 statements), `corollary` |
| Check.lean | #print axioms for 82 theorems |
| gen.py, gen2.py | generators (deterministic: the clean rebuild regenerates all files and compares hashes) |

Order names: + p, × x, − m, ÷ d (D_pxmd = + × − ÷, the paper's cycle).

## Hypotheses

None. Every input is checked inside Lean: the 132 patches of patches.json (`patchOK`, 132 kernel evaluations), the base sequences N0 ≤ N < N1 and the sporadic witnesses found by verify.py's `search` (`check`, 411 kernel evaluations), the non-existence for small N (`noPerfect`) and the dead ends of the 13 failing orders (`reach … = []`). The Python files only produce data; nothing they claim is trusted.

## Main statements (verbatim)

Definitions (Defs.lean; docstrings dropped):

```lean
inductive Op where
  | add | mul | sub | div
  deriving DecidableEq, Repr

def Op.Rel : Op → Nat → Nat → Nat → Prop
  | .add, a, b, v => v = a + b
  | .mul, a, b, v => v = a * b
  | .sub, a, b, v => v + b = a
  | .div, a, b, v => 0 < b ∧ v * b = a

structure Order where
  s1 : Op
  s2 : Op
  s3 : Op
  s4 : Op
  deriving DecidableEq, Repr

def Order.IsPerm (σ : Order) : Prop :=
  σ.s1 ≠ σ.s2 ∧ σ.s1 ≠ σ.s3 ∧ σ.s1 ≠ σ.s4 ∧ σ.s2 ≠ σ.s3 ∧ σ.s2 ≠ σ.s4 ∧ σ.s3 ≠ σ.s4

def Order.get (σ : Order) : Nat → Op
  | 0 => σ.s1
  | 1 => σ.s2
  | 2 => σ.s3
  | _ => σ.s4

def Order.op (σ : Order) (i : Nat) : Op := σ.get (i % 4)

structure Move where
  a : Nat
  b : Nat
  v : Nat
  deriving DecidableEq, Repr

def Valid (σ : Order) (N : Nat) : List Nat → Nat → List Move → Prop
  | _, _, [] => True
  | S, i, x :: xs => x.a ∈ S ∧ x.b ∈ S ∧ (σ.op i).Rel x.a x.b x.v ∧ x.v ∉ S ∧
      1 ≤ x.v ∧ x.v ≤ N ∧ Valid σ N (x.v :: S) (i + 1) xs

def made (ms : List Move) : List Nat := 2 :: ms.map Move.v

def Perfect (σ : Order) (N : Nat) (ms : List Move) : Prop :=
  ms.length + 1 = N ∧ Valid σ N [2] 0 ms ∧ ∀ x, 1 ≤ x → x ≤ N → x ∈ made ms

-- validB, memB, relB, allB: Bool versions; check:
def check (σ : Order) (N : Nat) (ms : List Move) : Bool :=
  (ms.length + 1 == N) && validB σ N [2] 0 ms && allB (fun x => memB x (made ms)) N
theorem check_iff (σ : Order) (N : Nat) (ms : List Move) :
    check σ N ms = true ↔ Perfect σ N ms
```

Corollary and Theorem 1 (Main.lean):

```lean
theorem corollary (σ : Order) (hσ : σ.IsPerm) :
    (∃ N0, ∀ N, N0 ≤ N → ∃ ms, Perfect σ N ms) ↔
      ((σ.s3 = .sub ∨ σ.s4 = .sub) ∧ σ ≠ ⟨.div, .mul, .sub, .add⟩)

-- Theorem 1, all 24 orders (fields s1 s2 s3 s4):
theorem thm1_pxmd (N : Nat) : (∃ ms, Perfect (⟨.add, .mul, .sub, .div⟩ : Order) N ms) ↔ N ∈ [] ∨ 9 ≤ N
theorem thm1_pxdm (N : Nat) : (∃ ms, Perfect (⟨.add, .mul, .div, .sub⟩ : Order) N ms) ↔ N ∈ [] ∨ 10 ≤ N
theorem thm1_pmxd (N : Nat) : (∃ ms, Perfect (⟨.add, .sub, .mul, .div⟩ : Order) N ms) ↔ N ∈ []
theorem thm1_pmdx (N : Nat) : (∃ ms, Perfect (⟨.add, .sub, .div, .mul⟩ : Order) N ms) ↔ N ∈ []
theorem thm1_pdxm (N : Nat) : (∃ ms, Perfect (⟨.add, .div, .mul, .sub⟩ : Order) N ms) ↔ N ∈ [] ∨ 9 ≤ N
theorem thm1_pdmx (N : Nat) : (∃ ms, Perfect (⟨.add, .div, .sub, .mul⟩ : Order) N ms) ↔ N ∈ [4, 6, 10] ∨ 12 ≤ N
theorem thm1_xpmd (N : Nat) : (∃ ms, Perfect (⟨.mul, .add, .sub, .div⟩ : Order) N ms) ↔ N ∈ [] ∨ 9 ≤ N
theorem thm1_xpdm (N : Nat) : (∃ ms, Perfect (⟨.mul, .add, .div, .sub⟩ : Order) N ms) ↔ N ∈ [] ∨ 8 ≤ N
theorem thm1_xmpd (N : Nat) : (∃ ms, Perfect (⟨.mul, .sub, .add, .div⟩ : Order) N ms) ↔ N ∈ []
theorem thm1_xmdp (N : Nat) : (∃ ms, Perfect (⟨.mul, .sub, .div, .add⟩ : Order) N ms) ↔ N ∈ []
theorem thm1_xdpm (N : Nat) : (∃ ms, Perfect (⟨.mul, .div, .add, .sub⟩ : Order) N ms) ↔ N ∈ [4, 5, 6, 8, 9, 10] ∨ 12 ≤ N
theorem thm1_xdmp (N : Nat) : (∃ ms, Perfect (⟨.mul, .div, .sub, .add⟩ : Order) N ms) ↔ N ∈ [4, 5, 6] ∨ 12 ≤ N
theorem thm1_mpxd (N : Nat) : (∃ ms, Perfect (⟨.sub, .add, .mul, .div⟩ : Order) N ms) ↔ N ∈ []
theorem thm1_mpdx (N : Nat) : (∃ ms, Perfect (⟨.sub, .add, .div, .mul⟩ : Order) N ms) ↔ N ∈ []
theorem thm1_mxpd (N : Nat) : (∃ ms, Perfect (⟨.sub, .mul, .add, .div⟩ : Order) N ms) ↔ N ∈ []
theorem thm1_mxdp (N : Nat) : (∃ ms, Perfect (⟨.sub, .mul, .div, .add⟩ : Order) N ms) ↔ N ∈ []
theorem thm1_mdpx (N : Nat) : (∃ ms, Perfect (⟨.sub, .div, .add, .mul⟩ : Order) N ms) ↔ N ∈ []
theorem thm1_mdxp (N : Nat) : (∃ ms, Perfect (⟨.sub, .div, .mul, .add⟩ : Order) N ms) ↔ N ∈ []
theorem thm1_dpxm (N : Nat) : (∃ ms, Perfect (⟨.div, .add, .mul, .sub⟩ : Order) N ms) ↔ N ∈ [2, 3, 4] ∨ 9 ≤ N
theorem thm1_dpmx (N : Nat) : (∃ ms, Perfect (⟨.div, .add, .sub, .mul⟩ : Order) N ms) ↔ N ∈ [2, 3, 4] ∨ 12 ≤ N
theorem thm1_dxpm (N : Nat) : (∃ ms, Perfect (⟨.div, .mul, .add, .sub⟩ : Order) N ms) ↔ N ∈ [2, 4, 5, 6, 8, 9, 10] ∨ 12 ≤ N
theorem thm1_dxmp (N : Nat) : (∃ ms, Perfect (⟨.div, .mul, .sub, .add⟩ : Order) N ms) ↔ N ∈ [2, 4, 5]
theorem thm1_dmpx (N : Nat) : (∃ ms, Perfect (⟨.div, .sub, .add, .mul⟩ : Order) N ms) ↔ N ∈ [2]
theorem thm1_dmxp (N : Nat) : (∃ ms, Perfect (⟨.div, .sub, .mul, .add⟩ : Order) N ms) ↔ N ∈ [2]
```
Key lemmas:

```lean
theorem perfect_extend {σ : Order} {P N : Nat} {pre ext : List Move}
    (hpre : Perfect σ P pre) (hPN : P ≤ N)
    (hs : Steps σ (fun y => 1 ≤ y ∧ y ≤ P) (P - 1) ext)
    (hnd : (ext.map Move.v).Nodup)
    (hmem : ∀ y, y ∈ ext.map Move.v ↔ P < y ∧ y ≤ N)
    (hlen : ext.length = N - P) :
    Perfect σ N (pre ++ ext)

theorem bulk_steps {σ : Order} (hσ : σ.IsPerm) {i P d : Nat} (hmul : σ.op i = .mul)
    (hJ : 6 + 2 * d ≤ P) {S : Nat → Prop} {i' : Nat} (hi : i' % 4 = i % 4)
    (hS : ∀ y, 1 ≤ y → y ≤ P → S y) (hn3 : S (nj P 3)) : Steps σ S i' (bulk σ i P d)

theorem patch_sound {σ : Order} (hσ : σ.IsPerm) {r c : Nat} {p : Patch}
    (hok : patchOK σ r c p = true) {q : Nat} (hq : p.qstar ≤ q) {pre : List Move}
    (hpre : Perfect σ (4 * q + r) pre) :
    Perfect σ (3 * (4 * q + r) + c) (pre ++ ext σ r p q)

theorem works_of {σ : Order} (hσ : σ.IsPerm) (r N0 N1 : Nat) (h0 : 1 ≤ N0) (pt : Nat → Patch)
    (offs : List Nat) (hok : ∀ c ∈ offs, patchOK σ r c (pt c) = true)
    (hbase : ∀ N, N0 ≤ N → N < N1 → ∃ ms, Perfect σ N ms)
    (hcov : ∀ N, N1 ≤ N → ∃ c ∈ offs, ∃ q, N = 3 * (4 * q + r) + c ∧ (pt c).qstar ≤ q ∧
      N0 ≤ 4 * q + r) :
    ∀ N, N0 ≤ N → ∃ ms, Perfect σ N ms

theorem noPerfect_sound {σ : Order} {N : Nat} (h : noPerfect σ N = true) :
    ∀ ms, ¬ Perfect σ N ms

theorem bound_of_reach {σ : Order} {k : Nat} (h : reach σ (fun _ => true) 0 [[2]] k = []) :
    ∀ N ms, Perfect σ N ms → N ≤ k

theorem works : ∀ N, 9 ≤ N → ∃ ms, Perfect ord N ms

theorem fail_dxmp : ∀ N ms, Perfect (⟨.div, .mul, .sub, .add⟩ : Order) N ms → N ≤ 5
```
## Fidelity

- Paper versions. The rule lines cited in NOTE.md (l. 110–129, 151–153, 921–925) are v1 numbering (e-print 2610.10095v1, 8 Oct). v2 (31 Oct; amsd_v2_revised.tex) has the same rules, schedule and definition of perfect word for word (whitespace-normalized diff), at l. 70–86 and l. 101.
- v2 still states the question as open, verbatim (v2 l. 702 = v1 l. 922–925); it answers no part of it (no other order is discussed). v2 adds a patch search (l. 640–660): within its own patch shape (only the last d ∈ {1, 2} rounds rewritten, 2 ≤ q ≤ 40) no c = 3 patch exists; our c = 3 patch also rewrites head rounds 0–2, so it is outside that shape and not contradicted, but the side-result wording "the paper found none" should now cite v2 l. 651, 660.

| Lean | paper (v2) / NOTE.md | note |
|---|---|---|
| `Valid σ N [2] 0 ms` | start with 2 (l. 86 "all later results use s = 2"); NOTE.md §1 | collection = list S of numbers made so far |
| `x.a ∈ S ∧ x.b ∈ S` | "choose two numbers already made (possibly the same one twice)" (l. 72) | a = b allowed |
| `(σ.op i).Rel`, `σ.op i = σ.get (i % 4)`, i = m − 1 | "apply the operation required by the cycle" (l. 73), schedule ∘_m (l. 76–85) with σ1σ2σ3σ4 for + × − ÷ | move m uses σ_k, k ≡ m (mod 4): the question l. 702 asks about other orders, NOTE.md §1 |
| `.sub: v + b = a` | a − b | operand order a − b; v natural, and v ≥ 1 by the range condition, so the difference is positive |
| `.div: 0 < b ∧ v * b = a` | "a division must give a positive integer" (l. 73) | a ÷ b = v exactly; v ≥ 1 by the range condition |
| `x.v ∉ S ∧ 1 ≤ x.v ∧ x.v ≤ N` | "every single move must make a previously unconstructed number between 1 and N" (l. 101) | |
| `ms.length + 1 = N` | "in exactly N − 1 moves" (l. 101) | forces N ≥ 1; N = 0 is never perfect |
| `∀ x, 1 ≤ x → x ≤ N → x ∈ made ms` | "makes all of 1, …, N" (l. 101) | N = 1 is not perfect (1 ∉ {2}) |
| `Order.IsPerm` | "orders of the four operations" (l. 702) | four distinct entries |
| `corollary` | NOTE.md Corollary | "works" = perfect sequences for all large N (∃ N0 ∀ N ≥ N0) |
| `thm1_<order>` | NOTE.md §1 table (Theorem 1) | sets generated from verify.py's CLAIM, now proved in Lean |

Numbers are naturals: every number made in a perfect sequence is 2 or a result in [1, N], so integers or fractions never occur and the Nat encoding loses nothing. Lemma 2 and Lemma 3 of NOTE.md §4 are proved in general (all P, all q ≥ q*), not instance by instance; only the finite data (patches, base and small cases) is checked by evaluation.

## Proof route (differences from NOTE.md)

- Failing orders (NOTE.md §2): `bound_of_reach` plus the kernel replaces the hand arguments: the collections reachable by legal moves that each make a new number (no bound on the numbers; `reach`, repeats merged via sorted insertion) run out after 1 move (− first), 2 moves (+ −, × −, ÷ − first) or 5 moves (÷ × − +), so a perfect N has N ≤ 1, 2 or 5: the dead ends of NOTE.md §2 and verify.py part F.
- Working orders (NOTE.md §3–§5): `perfect_extend` reduces an extension to operands made and operations right (`Steps`), distinct results, results = (P, N], length N − P. The bulk is Lemma 2 in general (`bulk_steps`, `bulk_vals`, `bulk_nodup`; the numerator n_j comes from round ⌊(2j + 1 + P mod 2)/3⌋, NOTE.md's k(j)). `patchOK` checks Lemma 3's conditions: move P is × (σ.op (r + 3) = ×), 12 head moves, c + 4d + 2(r mod 2) tail moves, operations by the schedule, identities in q, every operand an earlier result or in [1, P] (head) or [1, P] ∪ X ∪ Y (tail) for all q ≥ q* (slope and value at q*, equivalent to verify.py's thresholds), n_3 a head result, J ≥ 2 at q*, head + tail results a permutation of the pool L ∪ M ∪ U (List.isPerm). `patch_sound` proves Lemma 3 for every patch that passes. `works_of` is the strong induction of Theorem 2 with N1 from patches.json; the base range N0 ≤ N < N1 uses verify.py's search sequences, each accepted by `check`.
- Small N (NOTE.md §6): `noPerfect` is the same search with results ≤ N (no other pruning), proved sound; sporadic N have witnesses accepted by `check`.
- Theorem 1 = small cases + Theorem 2 (working orders) or the dead-end bound (failing orders); the Corollary follows by cases over the 24 orders (`IsPerm` excludes the other 232 quadruples).
- Not needed: verify.py's division pruning and brute-force cross-check, NOTE.md's facts (F1), (F2).

## Build

One module at a time, plain `lean` (no lake process, no Mathlib), oleans to .scratch/orderslean/olean (build output, deleted at the end; `bash Riemann/Deep/Orders/Lean/audit/verify.sh` regenerates the data files, checks they are byte-identical, and rebuilds all 26 modules in about 7 min):

    d=Riemann/Deep/Orders/Lean   # per module m, in the order of the table
    LEAN_PATH=.scratch/orderslean/olean bash .scratch/memcap.sh 290 1800 /usr/bin/time -l \
      ~/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean -j1 -o .scratch/orderslean/olean/$d/$m.olean $d/$m.lean

Clean rebuild 12 Nov 17:24–17:31 UTC: the generators reproduced every .lean file byte for byte; 26 modules rc = 0, 433 s in total; peak memory footprint per module (/usr/bin/time -l, the top MEM measure; memcap's 3 s sampling missed short runs), maximum 258 MB (D_dxpm). Check.lean recompiled 17:33 with the sanity examples: 8 s, 86 MB.

| module | s | MB | module | s | MB |
|---|---|---|---|---|---|
| Defs | 8 | 95 | D_dxpm | 29 | 258 |
| Ext | 8 | 97 | D_dpxm | 22 | 227 |
| Bulk | 24 | 254 | D_dpmx | 22 | 228 |
| Patch | 24 | 187 | Search | 6 | 93 |
| Induct | 8 | 86 | Fail | 7 | 96 |
| D_pxmd | 15 | 205 | Small1 | 6 | 126 |
| D_pxdm | 21 | 210 | Small2 | 13 | 215 |
| D_pdxm | 25 | 231 | Small3 | 7 | 213 |
| D_pdmx | 33 | 235 | Small4 | 8 | 96 |
| D_xpmd | 25 | 199 | Small5 | 7 | 107 |
| D_xpdm | 28 | 205 | Small6 | 11 | 234 |
| D_xdpm | 25 | 243 | Main | 18 | 138 |
| D_xdmp | 25 | 247 | Check | 8 | 84 |
## Static scan and axioms

    rg -w 'sorry|admit|native_decide|skipKernelTC' Riemann/Deep/Orders/Lean   # no match (rc 1)
    rg '^\s*axiom\s' Riemann/Deep/Orders/Lean                                # no match (rc 1)

`decide +kernel` (the kernel evaluates a Bool to true; no extra axiom): 694 uses = 132 patches + 381 base sequences + 30 sporadic witnesses + 138 small-N non-existence cases + 13 dead ends. Check.lean: 82 #print axioms lines, all 82 standard (55 × [propext, Classical.choice, Quot.sound], 27 × [propext, Quot.sound]; no sorryAx, no Lean.ofReduceBool):

    'Orders.check_iff' depends on axioms: [propext, Quot.sound]
    'Orders.perfect_extend' depends on axioms: [propext, Classical.choice, Quot.sound]
    'Orders.bulk_steps' depends on axioms: [propext, Classical.choice, Quot.sound]
    'Orders.bulk_vals' depends on axioms: [propext, Classical.choice, Quot.sound]
    'Orders.patch_sound' depends on axioms: [propext, Classical.choice, Quot.sound]
    'Orders.works_of' depends on axioms: [propext, Classical.choice, Quot.sound]
    'Orders.reach_sound' depends on axioms: [propext, Quot.sound]
    'Orders.noPerfect_sound' depends on axioms: [propext, Quot.sound]
    'Orders.bound_of_reach' depends on axioms: [propext, Classical.choice, Quot.sound]
    'Orders.corollary' depends on axioms: [propext, Classical.choice, Quot.sound]

and the same for: D_pxmd.works, D_pxdm.works, D_pdxm.works, D_pdmx.works, D_xpmd.works, D_xpdm.works, D_xdpm.works, D_xdmp.works, D_dxpm.works, D_dpxm.works, D_dpmx.works, fail_pmxd, fail_pmdx, fail_xmpd, fail_xmdp, fail_mpxd, fail_mpdx, fail_mxpd, fail_mxdp, fail_mdpx, fail_mdxp, fail_dxmp, fail_dmpx, fail_dmxp, small_pxmd, small_pxdm, small_pmxd, small_pmdx, small_pdxm, small_pdmx, small_xpmd, small_xpdm, small_xmpd, small_xmdp, small_xdpm, small_xdmp, small_mpxd, small_mpdx, small_mxpd, small_mxdp, small_mdpx, small_mdxp, small_dpxm, small_dpmx, small_dxpm, small_dxmp, small_dmpx, small_dmxp, thm1_pxmd, thm1_pxdm, thm1_pmxd, thm1_pmdx, thm1_pdxm, thm1_pdmx, thm1_xpmd, thm1_xpdm, thm1_xmpd, thm1_xmdp, thm1_xdpm, thm1_xdmp, thm1_mpxd, thm1_mpdx, thm1_mxpd, thm1_mxdp, thm1_mdpx, thm1_mdxp, thm1_dpxm, thm1_dpmx, thm1_dxpm, thm1_dxmp, thm1_dmpx, thm1_dmxp.

Sanity examples in Check.lean (`decide`): the paper's Table (v2 l. 110–119, N = 9) is accepted by `check`; the same moves with 2 − 8 for 8 − 2, under × + − ÷, or one move short are rejected.
## Review

One `reviewer` subagent (12 Nov 17:32 UTC), given the paper's rule lines (v2 l. 70–86, 101, 702), NOTE.md §1 (rules, table, Corollary) and the Lean definitions and main statements verbatim: **FAITHFUL**, no issues. Checked: schedule (`σ.get (i % 4)` gives σ1..σ4 for moves m = i + 1, no shift); IsPerm scope; operands from the prior collection, reuse and a = b allowed; − and ÷ ordered and exact, positivity via 1 ≤ v; naturals lose no perfect sequence; start {2}, newness, coverage; `ms.length + 1 = N` (N = 0 impossible, N = 1 fails, N = 2 exactly for ÷-first orders); all 24 Theorem 1 statements match NOTE.md; the Corollary's ∃ N0 ∀ N ≥ N0 is "works for all large N". Optional remark: the counts (11 working orders; finite sets ∅, {2}, {2, 4, 5}) follow from thm1_* but are not separate lemmas. Prompt and answer: Riemann/Deep/Orders/Lean/audit/review_prompt.txt, review_out.md.

## Entry for FORMALIZATIONS.md

- Operation orders (side result; open question of arXiv 2610.10095): Lean 4 core, sorry-free, standard axioms — the exact perfect-N set of each of the 24 orders of + × − ÷ and the Corollary (exactly 11 orders work); patches, base and small cases kernel-checked by proved-sound checkers → [LEAN.md](Riemann/Deep/Orders/LEAN.md)

## Housekeeping

- The repo-root .lake (Mathlib cache of other agents) was left as it is: no lake command was run and nothing imports Mathlib. The orchestrator restores the placeholder.
- Deleted at the end: the TeX downloads (/tmp/ops.tgz, /tmp/ops_v1.tgz, extracted sources), the test files, the oleans (.scratch/orderslean/olean; verify.sh rebuilds them), and the `__pycache__` that importing verify.py had created in Riemann/Deep/Orders/ (the generators now set `sys.dont_write_bytecode`).
- Logs: Riemann/Deep/Orders/Lean/audit/peaks.log (every compile of the session), verify.log, logs/log_<module>.txt, err_<module>.txt, STATE.md.

