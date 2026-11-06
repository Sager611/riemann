RECORD: from the Lean-checked half-plane (B = 0.8733819666132076, E = 1/(1−B) = 7.8978), for every Q ≥ 3 there is a set P of at most 2.1067 log Q + 58 primes, all ≤ (820·max(log Q, log 9))^E, such that for every 2 ≤ q ≤ Q the subset products of the primes of P not dividing q are all of (ℤ/qℤ)^× — the first covering of every modulus by O(log Q) small primes with numbers, = (1/log 2 + 0.664) log Q + 58 (per-Q table: Q = 10¹²: 77 primes ≤ 10^33.33; Q = 10¹⁰⁰: 498 primes ≤ 10^39.91).

## Claim
Assume (H): ζ(s) and every Dirichlet L(s, χ) are zero-free on Re s > B (a theorem: Bebel's Lean-kernel-checked Re s > B₂ + 2.2·10⁻²⁵⁴,
B₂ = 0.87338196661320752 < B). Say P *covers* q if {∏_{p∈F} p mod q : F ⊆ P, every p ∈ F coprime to q} = (ℤ/qℤ)^×.

**Theorem 1 (every modulus, every Q).** For every Q ≥ 3 there is a set P of primes, all in (e^{−7}y, y] with y = (c'·max(log Q, log 9))^E,
that covers every 2 ≤ q ≤ Q, with
| c' | #P ≤ | = (1/log 2 + η) log Q + c |
|---|---|---|
| 820 | **2.1067 log Q + 58** (k = 4); 2.2237 log Q + 44 (k = 3); 2.4624 log Q + 31 (k = 2); 3.1832 log Q + 15 (k = 1) | η = 0.6640, 0.7810, 1.0197, 1.7405 |
| 820 | 2.1067 log Q + 15 for Q ≥ 10¹⁰⁰ | η = 0.6640 |
| 820000 | **1.9492 log Q + 98** (k = 4); 2.0694 log Q + 74 (k = 3) | η = 0.5065 |
In particular P covers every q ≤ Q coprime to ∏P (Hughes's form). P is only shown to exist: N independent draws from the prime
law ν of the Method work with probability ≥ 1 − U ≥ 0.45 (k = 4 row), but no explicit description is known and checking a draw for
all q ≤ Q is infeasible for large Q.

**Theorem 2 (table, least certified count, c' per row).**
| Q | log₂ Q | #P ≤ | primes ≤ | c' |
|---|---|---|---|---|
| 10³ | 10.0 | 27 (explicit: the first 15 primes, ≤ 47, suffice) | 10^29.34 | 750 |
| 10⁶ | 19.9 | 44 | 10^31.33 | 670 |
| 10¹² | 39.9 | 77 | 10^33.33 | 600 |
| 10³⁰ | 99.7 | 164 | 10^36.05 | 530 |
| 10¹⁰⁰ | 332.2 | 498 (470 with primes ≤ 10^63.60) | 10^39.91 | 490 |
| 10¹⁰⁰⁰ | 3321.9 | 4858 (4510 with primes ≤ 10^71.28) | 10^47.59 | 460 |
Any covering set needs ≥ log₂ φ(q) primes, i.e. about log₂ Q.

**Proposition 3 (small Q, exact).** The first 17 primes 2, …, 59 cover every q ≤ 2000 (first 15: q ≤ 1000; 13: q ≤ 300; 11: q ≤ 100);
the 16 primes 2003, …, 2113 cover every q ≤ 2000.

**Not reached:** the brief's (1/log 2 + η) with η small. The method gives slope (1 + 1/k)/log(1/ρ_g) + 2/log(1/β) → 5/(4 log 2) = 1.803
for k = 4 as c' → ∞; η → 0 would need the k-th moment bound for every k (DERIVATION §9).

## Before
- Klurman–Shparlinski–Teräväinen (arXiv 2609.25393, Thm 1.1): (log Q)^{1+ε} primes ≤ (log Q)^{C_*(A+1)} covering all but O(Q(log Q)^{−A})
  moduli; C_* not computed; the exceptional set "is unavoidable with present technology" (TeX l.129). Their GRH sketch (l.101–105)
  covers every q with all primes ≤ (log q)^B, i.e. polylogarithmically many primes.
- Hughes (arXiv 2610.04678): (1/log 2 + o(1)) log Q primes for all but O_A(Q(log Q)^{−A}) moduli (Thm main; one exceptional set comes from
  random failures and survives any zero-free region); 2 log Q/log 2 + O_A(log Q/log log Q) primes for every q ∉ 𝓔_an coprime to ∏P
  (Thm structured), Q sufficiently large, "We do not track effectivity". Under (H) 𝓔_an = ∅ (every character is good, his (ZF) with
  σ < 1 − B), which gives an all-moduli statement with 2/log 2 + o(1) and no numbers; nobody wrote it down.
- No all-moduli statement with O(log Q) primes and numbers existed; no all-moduli leading constant below 2/log 2 existed.

## Method
DERIVATION.md has every step and constant. In short:
- Draw N primes from the law ν(p) ∝ log p·w(p/y) on (e^{−7}y, y] (explprog's triangle kernel); for a modulus q a draw dividing q or
  repeating an earlier draw is a no-op. Lemma D: for every q ≤ Q, every nonprincipal ψ mod q and every step,
  Re E[ψ(g_j) | past] ≤ δ, from the explicit formula with all zeros in b ≤ β ≤ B (Lemmas EF, P, Z, Zζ, K1 of explprog) plus the
  no-op and prime-power masses; δ ≤ 10⁻³ for all Q ≥ 3 at c' = 820 (664 certified cells and a monotone tail). This is where the
  half-plane empties Hughes's exceptional set: no character needs repair.
- Lemma C, collision stage: the k-th moment of Hughes's collision count C_r = Σ_{χ≠1}|μ̂_r(χ)|² is a sum over k-tuples of characters of
  (π + δ(1 − π))^r, π = the probability of a lazy ±1 relation among the tuple; tuples are grouped by their closed relation set
  (enumerated exactly: 1, 4, 48, 2340 sets for k = 1..4) with kernel counts H^{k−ρ}∏|G[d_i]| (Smith form). Pr(hole ≥ αH) ≤ E[C_r^k]/t^k.
  For k = 2, 3 a hand classification gives the same classes. Hole stage: Hughes's contraction E|T ∩ gT| ≤ (δ + (1 − δ)α)|T| while |T| < αH.
- Union over 3 ≤ q ≤ Q with explicit sums of φ(q)^m∏|G[d_i]| (|G[d]| ≤ 2d^{ω(q)}, Σ τ_D ≤ Q(1 + log Q)^{D−1}); the k-th moment makes
  the collision stage cost (1 + 1/k) log₂ Q instead of 2 log₂ Q (k = 1, Hughes's lemma).
- Theorem 1: r = ⌈a₁ log Q + a₀⌉, s = ⌈b₁ log Q + b₀⌉ and a certified sup over all Q ≥ 3 of the failure sum (< 0.55 for k = 4);
  every printed constant is a ball upper end rounded up.

## Certificate
One command (about 50 s, 35 MB peak, one thread; python-flint Arb at 200 bits; every inequality a certified ball sign, else abort):

    cd Riemann/Deep/UnitCover && uv run --no-project --with python-flint --with mpmath python3 cert.py all

Output: logs/cert_output.txt (sections core, delta, cover, table, small). It imports Riemann/Deep/ExplProg/cert.py unchanged.

## Checks
- classcheck.py (exact integers, 610 s): for every tuple of nonprincipal characters of 118 (k = 2), 66 (k = 3), 29 (k = 4) small groups
  (ℤ/qℤ)^×, its relation set is one of the enumerated closed sets, π matches, the class count is ≤ the Smith-form kernel, the
  torsion bounds |G[2]| ≤ τ(q), |G[d]| ≤ 2d^{ω(q)} hold, and the hand classes for k = 2, 3 dominate (logs/classcheck_output.txt).
- Exact closure for every q ≤ 2000 with the printed sets (Prop. 3), incl. 3, 4, 8, 30, 210, 1024, 1155 = 3·5·7·11, 1365, primes 1993–1999.
- Monte Carlo (mc.py, not proof): random primes from (10⁶, 2·10⁶] cover all q ≤ 100, 300, 1000 with median N = 10, 12, 14; the
  certified N(10³) = 27 is about twice that (logs/mc_output.txt).
- Lemma D: cell bounds agree with direct evaluation (δ = 0.00095 at log Q = log 9, 0.00059 at log Q = 230); the tail margin is 34%.
- enumcheck.py (56 s): brute force over all symmetric subsets of the cube reproduces the 4 and 48 closed sets (k = 2, 3); 32403 random
  and structured 4-tuples in 12 groups with element orders up to 840 meet 691 distinct relation sets, all enumerated
  (logs/enumcheck_output.txt).
- Three verdict-only reviewer passes (6 Nov, inline statements): (1) moment lemma: steps sound, but δ ≥ 0 must be assumed (else β < 0
  and the over-count can fail) and the finite-horizon bound on δ must include the repeat mass — both hold here (Lemma D has
  δ > 0 and prices N − 1 earlier draws); the Smith-form count and closedness were confirmed. (2) Lemma D: proved, y-monotonicity
  confirmed; its remark that the cell bound needs no L-monotonicity is now explicit (Lemma Z at L₁ covers every conductor ≤ e^{L₁};
  den is certified at the cell corner). (3) sums over moduli, torsion bounds, the linear choice of r, s and the sup over Q: SOUND
  (budget N − 1 ≤ 2.10668 L + 57 ≤ 4L + 200 confirmed). DERIVATION §2, §3 updated.
- Self-review against the brief: moduli with many small primes (no-ops cost ≤ 2‖w‖_∞ log q/M), q ≤ 4 (H = 2), q prime, 2-powers,
  repeated draws, imprimitive characters; no GRH, no zero-density estimate, no unrefereed zero count; all roundings upward.

## Novelty
Searched 6 Nov 2026 (NOVELTY.md): KST and Hughes TeX (both v1), Frei–Sofos–Widmer 2610.18491 v2 TeX (no covering), arXiv math.NT and
math.CO listings Oct–Nov 2026 (two regexes; only Hughes, FSW and half-plane papers), Zenodo (11 phrase queries; Tolbert 23581412 is
about Regev relation lattices for almost all inputs: different, cited), HAL (4 queries, nothing), GXL evidence API (2 claims: only
KST/Hughes with exceptional sets), our repo (only the harvest entry). Nothing states an all-moduli covering with numbers or an
all-moduli constant below 2/log 2.

## Caveat
- Conditional only on (H), which is Lean-kernel-checked (tracker admission) but not refereed. Not refereed here; three short reviewer
  passes only (see Checks).
- Existence only; the primes are huge for practical Q (10^26 to 10^48 in Theorems 1–2; with c' = 820 they exceed Q for Q below 10^38.4): the content is asymptotic.
  For small Q exact sets are far better (17 primes ≤ 59 for all q ≤ 2000).
- The constants are far from the truth: observed ≈ log₂ Q + 4 primes; certified slope ≤ 1.4603 log₂ Q (c' = 820). Losses:
  union over ≈ Q moduli (k ≤ 4), the hole stage (2/log(1/δ)), and the BMOR zeros below height 5/7 (72% of c_∞ = 0.4533).
- The k = 4 classification is computer-assisted (exact enumeration), verified exhaustively on small groups but not by hand.
- Routine-level parts: Theorem 1 with k = 1 is Hughes's argument with explicit constants; the gain to k = 4 is a new but elementary
  refinement.

RH relevance: N11

## Files
- REPORT.md (this file), DERIVATION.md (all steps and constants), NOVELTY.md (searches), cert.py (Arb certificate: core, delta,
  cover, table, small), classcheck.py (exhaustive check of the tuple classes on small groups), enumcheck.py (independent check of the
  closed-set enumeration), mc.py (Monte Carlo sanity), logs/ (cert_output.txt, classcheck_output.txt, enumcheck_output.txt,
  mc_output.txt).
