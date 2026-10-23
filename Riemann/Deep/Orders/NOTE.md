# Which orders of + × − ÷ give perfect sequences

Answer to the open question of arXiv 2610.10095 (TeX l. 921–925). Files: this note,
REPORT.md, PATCHES.md (all patches, readable), patches.json (same data), verify.py
(all checks, about 6 s).

## 1. Result

Rules (TeX l. 110–129, 151–153). Start with the set {2}. An order σ = σ1σ2σ3σ4 is a
permutation of + × − ÷; move m applies σ_k with k ≡ m (mod 4) to two numbers already
made (possibly equal); a difference must be positive, a quotient a positive integer; the
result joins the set. A sequence is perfect for N if after N − 1 moves every number
1..N is made and every move made a new number in [1, N]. N = 1 is never perfect (1 is
not made by zero moves).

**Theorem 1 (classification).** The numbers N with a perfect sequence are exactly:

| order | perfect N | order | perfect N |
|---|---|---|---|
| + × − ÷ (the paper) | N ≥ 9 | ÷ + × − | 2, 3, 4, N ≥ 9 |
| + × ÷ − | N ≥ 10 | ÷ + − × | 2, 3, 4, N ≥ 12 |
| + ÷ × − | N ≥ 9 | ÷ × + − | 2, 4, 5, 6, 8, 9, 10, N ≥ 12 |
| + ÷ − × | 4, 6, 10, N ≥ 12 | ÷ × − + | 2, 4, 5 only |
| × + − ÷ | N ≥ 9 | ÷ − + ×, ÷ − × + | 2 only |
| × + ÷ − | N ≥ 8 | − first (6 orders) | none |
| × ÷ + − | 4, 5, 6, 8, 9, 10, N ≥ 12 | + − × ÷, + − ÷ × | none |
| × ÷ − + | 4, 5, 6, N ≥ 12 | × − + ÷, × − ÷ + | none |

**Corollary.** An order works (perfect sequences for all large N) if and only if
subtraction is third or fourth and the order is not ÷ × − +. So exactly 11 of the 24
orders work. For the other 13, perfect N form a finite set: ∅, {2} or {2, 4, 5}.

Proof plan: §2 the 13 failing orders (hand proofs); §3–§5 a construction for every
N ≥ N0 for the 11 working orders (one uniform bulk lemma, 132 small patches checked
generically in q, induction N = 3P + c); §6 the finitely many small N (exact search).

## 2. The thirteen failing orders (Proposition 1)

Two facts used throughout. (F1) A quotient a/b is a or at most a/2, so a division never
makes a number above half the largest number made. (F2) From a complete interval [1, P],
a difference or a quotient is never new.

1. σ1 = − (6 orders): move 1 can only be 2 − 2 = 0, not positive. No legal move, so no
   N ≥ 2 is perfect (and N = 1 never is).
2. σ1σ2 = +− or ×− (4 orders): move 1 makes 4 (2 + 2 or 2 × 2; nothing else is new).
   Move 2 is a subtraction in {2, 4}: the only positive difference is 4 − 2 = 2, old. So
   N ≥ 3 is impossible, and N = 2 needs move 1 to make 1, impossible.
3. σ1σ2 = ÷− (2 orders): move 1 makes 2/2 = 1; move 2 gives only 2 − 1 = 1, old. So
   only N = 2 (perfect by 2/2 = 1).
4. ÷ × − +: every legal start is forced: 2/2 = 1, then × in {1, 2} gives only 4 new, then
   − in {1, 2, 4} gives only 3 new, then + gives some s ∈ {5, 6, 7, 8}. Move 5 is a
   division in {1, 2, 3, 4, s}: by (F1) a new quotient would be at most s/2 ≤ 4, but
   1..4 are made. So no fifth move makes a new number; a perfect sequence for N needs
   N − 1 moves, hence N ≤ 5. N = 2, 4, 5 are perfect (stop after 1, 3, 4 moves with
   s = 5); N = 3 is not (move 2 is ×, it cannot make 3).

verify.py part F re-checks these dead ends by enumerating all legal prefixes.

## 3. Choice of P and the rounds

For a working order σ let r ≡ 1 + (index of × in σ, counted from 0) (mod 4). If P ≡ r
(mod 4), move P is a multiplication. After a perfect sequence for P the set is exactly
[1, P] and the next move is ×. (By (F2) the move after a complete interval must be
+ or ×; we always arrange ×.)
From move P on, the moves come in rounds of four starting with ×: × followed by the
other three operations in the cyclic order of σ (the "round type").

| σ | r | P | round type | offsets c used | N0 | N1 |
|---|---|---|---|---|---|---|
| + × − ÷ | 2 | even | × − ÷ + | 0..11 | 9 | 41 |
| + × ÷ − | 2 | even | × ÷ − + | 0..11 | 10 | 41 |
| + ÷ × − | 3 | odd | × − + ÷ | 0..11 | 9 | 44 |
| + ÷ − × | 0 | even | × + ÷ − | 0..11 | 12 | 47 |
| × + − ÷ | 1 | odd | × + − ÷ | 0..11 | 9 | 38 |
| × + ÷ − | 1 | odd | × + ÷ − | 0..11 | 8 | 38 |
| × ÷ + − | 1 | odd | × ÷ + − | 1..12 | 12 | 51 |
| × ÷ − + | 1 | odd | × ÷ − + | 1..12 | 12 | 51 |
| ÷ × + − | 2 | even | × + − ÷ | 0..11 | 12 | 53 |
| ÷ + × − | 3 | odd | × − ÷ + | 0..11 | 9 | 44 |
| ÷ + − × | 0 | even | × ÷ + − | 0..11 | 12 | 47 |

## 4. The extension [1, P] → [1, 3P + c]

Write P = 4q + r, H = ⌊3P/2⌋ (= 3P/2 or (3P − 1)/2), T = 3P, N = T + c, fix d ∈ {0, 1}
and let J = P/2 − 1 − d (P even) or (P − 3)/2 − d (P odd). If J ≥ 2, the intervals

  [1, P], L = [P+1, P+d], X = [P+1+d, H−3], M = [H−2, H+e+3d], Y = [H+e+1+3d, T−8], U = [T−7, T+c]

(e = 1 for P even, e = 3 for P odd) partition [1, N]. The N − P new moves are: the head
(rounds 0, 1, 2: twelve moves), the bulk (rounds j = 3..J), the tail (c + 4d moves, plus
2 if P is odd). Head and tail together make the pool L ∪ M ∪ U; the bulk makes X ∪ Y.

**Lemma 2 (bulk).** Suppose the set before round 3 contains [1, P] and n3, where
n_j = 3P − 2j (P even) or 3P − 1 − 2j (P odd). Let round j make, in round-type order,

  × : 3·(P − j) = C_j,   + : C_j + 1,   − : C_j − 1,   ÷ : n_j / 2 = H − j.

Then every bulk move is legal and the bulk makes each number of X ∪ Y exactly once.

*Proof.* × uses 3 and P − j ∈ [1, P] (j ≤ J < P/2). ± use C_j, made first in the same
round, and 1. The triples {C_j − 1, C_j, C_j + 1} are consecutive (C_{j−1} = C_j + 3)
and tile [C_J − 1, C_3 + 1] = Y; the quotients H − j (j = 3..J) are X; n_j is even.
The numerator n_j was made earlier: P even, j = 3t, 3t+1, 3t+2: n_j = C_{2t},
C_{2t+1} + 1, C_{2t+1} − 1, made in round k(j) = 2t, 2t+1, 2t+1; P odd: n_j =
C_{2t} − 1, C_{2t+1}, C_{2t+2} + 1, round 2t, 2t+1, 2t+2. For j ≥ 4, 3 ≤ k(j) ≤ j − 1,
so n_j is a bulk number of an earlier round; for j = 3 it is n3, given. ∎

**Patches.** All head and tail numbers are written as a·q + b with integers a, b. A
patch lists, for each head and tail move, the result t (a pool number) and two operands
u, w with t = u ∘ w as an identity in q (for ÷, w is a constant k ≥ 2 and u = k·t; for
×, u is a constant). Each operand has a witness: an earlier patch result, or a value
in [1, P], or (tail only) a value in X or Y. A witness "lo ≤ v ≤ hi" between two linear
functions of q holds for all q ≥ some explicit q0 (lexicographic comparison of (a, b));
q* is the maximum of these q0 and of the q for which J ≥ 2.

**Lemma 3 (patch).** If a patch satisfies: its operations follow the schedule from move
P on, its results are exactly the pool, n3 is a head result, and every operand has a
witness, then for every q ≥ q* a perfect sequence for P = 4q + r extends to one for
3P + c (prefix, head, bulk, tail).

*Proof.* For q ≥ q* the six intervals partition [1, N] (J ≥ 2), so pool, bulk and
[1, P] are disjoint and cover [1, N]: every move makes a new number of [1, N], each
once, in N − P moves. Operands: a head operand is in [1, P] (made by the prefix) or an
earlier head result; a tail operand may also be in X ∪ Y, made by the bulk, which is
legal by Lemma 2 (n3 is a head result). The identities give the right results; every
division is exact because u = k·t. ∎

All 132 patches (11 orders × 12 offsets) are listed in PATCHES.md and patches.json:
d ∈ {0, 1}, 12 to 30 moves, q* ≤ 2. verify.py part C checks every condition of Lemma 3
symbolically and recomputes q*; part D also builds and checks the full sequences for
q* ≤ q ≤ q* + 8. They were found by a depth-first search over pool assignments with
symbolic witnesses (bounded constants 1..4, sometimes 1..6 or 1..12).

**Example (the paper's order, c = 3).** P = 4q + 2, H = 6q + 3, T = 12q + 6,
N = 12q + 9; rounds × − ÷ +. The paper found no patch for c = 3 and used c = 15
(23 moves, two rounds deleted); this one has 15 moves, q* = 1, so it extends every
perfect P ≡ 2 (mod 4) (P ≥ 10) to 3P + 3.

| round | × | − | ÷ | + |
|---|---|---|---|---|
| 0 | 4·(3q+2) = T+2 | (T+2)−1 = T+1 | (T+2)/2 = H+1 | (T+2)+1 = T+3 |
| 1 | 3·(4q) = T−6 | (T+1)−1 = T | T/2 = H | (T−6)+2 = T−4 |
| 2 | 4·(3q+1) = T−2 | (T)−1 = T−1 | (T−2)/2 = H−1 | (T−6)+1 = T−5 |
| 3..2q | bulk, Lemma 2 (needs n3 = T−6, made in round 1) | | | |
| tail | 3·(4q+1) = T−3 | (T−6)−1 = T−7 | (T−4)/2 = H−2 | |

Pool: U = [T−7, T+3], M = [H−2, H+1]; bulk: X = [P+1, H−3], Y = [H+2, T−8].

## 5. Theorem 2: every N ≥ N0 is perfect (11 working orders)

*Proof.* Strong induction on N. For N0 ≤ N < N1 (table in §3) an exact search gives a
perfect sequence (verify.py part A re-finds and checks each one). Let N ≥ N1. The 12
offsets c of σ meet every residue class mod 12 once in 3r + c, so there is exactly one
c with N ≡ 3r + c (mod 12); put P = (N − c)/3 ≡ r (mod 4), q = (P − r)/4. N1 is the
largest, over c, of the least N in the class of c with P ≥ max(4q*_c + r, N0); hence
q ≥ q*_c and N0 ≤ P < N. P is perfect by induction, and Lemma 3 gives N. ∎

For σ = + × − ÷ this reproves the paper's theorem with a different patch set. The
construction is explicit: verify.py part E builds and checks it for every N0 ≤ N ≤ 400
and N = 1000, 4321, 10000, 12345, 30001; a one-off run checked every N ≤ 2500 and 15
random N < 200000 per order.

## 6. Small N (exact search)

verify.py part A runs an exhaustive depth-first search over the new number made at
each move, for all 24 orders and N ≤ 60, with memo of dead sets and one valid pruning
(each remaining division must make a missing number ≤ N/2, by (F1)). For N ≤ 14 a plain
brute force without memo or pruning agrees on all 336 cases. This gives the finite
parts of Theorem 1, e.g. N = 7 and N = 11 have no perfect sequence for × ÷ + − (search
trees of 9 and 61 nodes). Independently of the construction, direct search found
perfect sequences for every 12 ≤ N ≤ 160 and all 11 working orders.

## 7. Remarks

- The scout table was wrong in three places: N = 11 is not perfect for + ÷ − ×,
  × ÷ + − and ÷ × + − (N0 = 12, not 10, 8, 8), and it omitted the sporadic small N.
- Whether an order fails is decided in its first five moves. Among working orders the
  cyclic class only fixes the round type; the rest of the argument is uniform.
- verify.py: about 6 s, peak memory about 100 MB, no dependencies beyond Python 3.
- Open: counting perfect sequences; other starting numbers; removing an operation.
