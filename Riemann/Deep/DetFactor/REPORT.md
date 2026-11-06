KNOWN: OpenAI (QRH preprint, 30 Sep 2026: deterministic polynomial-time square roots mod p; "Deterministic Polynomial Factorization over Prime Fields", 4 Oct 2026: deterministic polynomial-time factorization over F_p in every degree, from a uniform Hecke strip) and Hexagon 2610.00178 (11 Oct 2026: roots of X^n − a, Rónyai's even and 3-divisible splitting, a small non-cube in F_p[ω], all without ERH). Our additions are routine: OpenAI's factoring theorem needs only Cureton–Thorner's Lean-checked strip, its auxiliary primes are ≤ (1.96 q² log(12q²p))^{32}, its exponent 10^{12} becomes 500000, and the constants of the small non-residues become numbers.

# DetFactor (node N11): deterministic finite-field algorithms without ERH

Agent detfactor, 6 Nov 2026, 11:26–12:40 UTC. Full arguments: [DERIVATION.md](DERIVATION.md); searches: [NOVELTY.md](NOVELTY.md).

## Claim (this round; routine, candidates for COROLLARIES)
1. **Numeric auxiliary primes (Thm B).** For primes q and p with p ∤ 6q there is a prime ℓ ≡ 1 (mod 12q), ℓ ∉ {2, 3, p, q},
   with p not a q-th power mod ℓ and ℓ ≤ (0.49(log 3 + D(q, p)))^{32} ≤ (1.96 q² log(12q²p))^{32},
   D = qφ(12q) log(12q) + φ(12q)(q log q + (q − 1) log p). Proof: our numeric Chebotarev theorem (Cureton–Thorner only)
   for the class Gal(L/ℚ(μ_{12q})) \ {1} of L = ℚ(μ_{12q}, p^{1/q}), with Δ_L ≤ e^D.
2. **OpenAI's factoring on a Lean-checked input (Prop. C).** The preprint's only analytic input (Section 10: a strip for
   ζ_K, ζ_M with K = ℚ(μ_{12q}), M = K(p^{1/q})) is supplied by Cureton–Thorner's Lean-checked ζ_K ≠ 0 on Re s > 31/32,
   so the companion's unrefereed Theorem 1.2 (width 10^{−6}) is not needed; with Thm B its bit bound
   O(B^{500000} + B^{100}E^{20}) becomes O(((n + 1)⌈log₂ p⌉)^{500000}) (stated: exponent 10^{12}).
3. **Explicit exponents in F_p (Thm A).** The least r-th power non-residue mod p is a prime ≤ (log p)^{7.8978} for every
   prime r | p − 1; deterministic square roots in O((log p)^{8.8978} log log p) bit operations, r-th roots in
   O((log p)^{10.8978} + r log³ p).
4. **F_{p^k}.** Lemma N (x is an r-th power in F_{Q^k} iff its norm is one in F_Q, r | Q − 1): for r | p − 1, r ∤ k the
   least prime r-th non-residue of F_{p^k} is ≤ (log p)^{7.8978}. Generators of the r-Sylow subgroup of any explicit
   field k ∋ ζ_r: OpenAI's table and Lemma ff:unitroot (eigenvector branch) with Thm B, plus the one case they leave
   out (r = 2, p ≡ 3 mod 4, v_2([k:F_p]) = 1: x = a + bi with a² + b² = ℓ ≡ 1 mod 4, (ℓ/p) = −1,
   ℓ ≤ (0.49(log 48 + 2 log p))^{32}).
5. **Small non-residues of F_{p²} (Thm E).** K imaginary quadratic, p inert, r | p² − 1, units r-th powers: a non-r-th
   power α ∈ O_K with N(α) ≤ (0.49(log 3 + h r (log|Δ_K| + 2 log p)))^{32}; e.g. a non-cube a + bω of F_p[ω] with
   a² − ab + b² ≤ (1.96 log 3p)^{32} (Hexagon: C(log 3p)^{32}, C not computed).
6. **Classical route (task 3).** With 4., the GRH subroutine of Rónyai 1988, Evdokimov 1994 and Ivanyos–Karpinski–Saxena
   2009 (primitive r-th non-residues and roots of x^r − a over the working field; IKS: "needed only" for this) is supplied
   without GRH and without OpenAI's odd-degree splitter, so Evdokimov's (n^{log n} log p)^{O(1)} holds unconditionally,
   taking IKS's description as faithful. Superseded by OpenAI's claimed polynomial bound.

Not achieved: small non-residues in a fixed polynomial basis of F_{p^k} for k ≥ 3 — blocking step: the residue character
of O_K/pO_K is a Hecke character only on the part trivial on global units, and for k ≥ 3 the unit group is infinite, so
small prime ideals from the zero-free region need not have small generators.

## Before
- Under ERH: Ankeny/Bach non-residues; Rónyai 1988 (poly(n^r, log q) for a prime r | deg f), Huang 1991 (abelian Galois
  group), Rónyai 1992, Evdokimov 1994 ((n^{log n} log q)^{O(1)}), Ivanyos–Karpinski–Saxena 2009; Shoup/Buchmann–Shoup
  non-residues in F_{p^k}. Unconditional: Shoup p^{1/2+o(1)}, IKRS 2012 (factor or automorphism), Schoof (fixed a).
- After the 2026 half-planes: square roots (OpenAI 30 Sep), AMM roots, Rónyai r = 2, 3 and a non-cube in F_p[ω]
  (Hexagon 2610.00178), complete factorization in polynomial time (OpenAI 4 Oct, analytic input δ = 10^{−6}, auxiliary
  primes ≤ c₀B^{2·10^7}, c₀ unspecified), F_{p^k} consequences (Zenodo 23205034). Our OAINum audit (R255, 5 Nov) found no
  flaw in the factoring preprint.

## Method (per algorithm: ERH step, source line, replacement; full table in DERIVATION §5)
- Tonelli–Shanks / Cipolla: quadratic non-residue → Thm A(a) (OpenAI QRH paper.tex l. 320–362: C(log p)^{32}).
- Adleman–Manders–Miller: r-th non-residue mod p → Thm A(a) (Hexagon 2610.00178 l. 10430–10470).
- Rónyai 1988, Evdokimov 1994, IKS 2009: "GRH … needed only to find primitive r-th nonresidues in a finite field F_q"
  (IKS 0804.1974 l. 94–97; IKRS 0811.3165 l. 150–167), in the field k ⊇ F_p with s-th non-residues for s ≤ m (IKS
  l. 465–466) → OpenAI 01-table + 02 Lemma ff:unitroot (l. 106–) + Thm B + our q = 2 case; r = 2, 3 already in Hexagon
  2610.00178 (l. 10605, 10735).
- Huang 1991, Rónyai 1992: same mechanism (IKS l. 100–105, 485); originals not read (paywalled).
- BIMS 2017: weak-GRH Dirichlet strip (Conj. 6.3, l. 732–735) ⇒ n(p, r) = O(log^{4/(1−2ε)} p); Bebel's B gives exponent
  15.8 there; Thm A gives 7.8978 with constant 1.
- OpenAI 2026: Section 10 only (10-analytic l. 9–18, 39–48, 91–98); rg of sections 02–09 finds no other analytic input
  → Cureton–Thorner + Thm B (Prop. C).

## Checks (all single-threaded, < 60 MB)
- logs/nr_output.txt (nr.c, 5 s): all primes p ≤ 10⁷: least quadratic, cubic, quintic non-residues ≤ 53, 29, 19;
  largest log n/log((log p)^{7.8978}) = 0.187 (p ≥ 17).
- logs/aux_output.txt (aux.py, 4 s): least auxiliary primes for all p ≤ 2·10⁵, q ≤ 13 are ≤ 2113 (≤ 20 q² log p); bound
  slack ≥ 10^{48}; adversarial p up to 10^{37.8} far below the bound.
- logs/table_output.txt (table.py, < 1 s): OpenAI's table for q = 2, 3 at 10 primes up to 10^{18} + 9, with negative
  controls; F_{p²} example p = 10⁹ + 7 (ℓ_3 = 181, ω_3 of order 9).
- logs/primary_output.txt (primary.py, 12 min): r-Sylow generators in 70 random fields k = F_p[t]/(g), r ≤ 7,
  [k:F_p] ≤ 14, every branch (eigenvector 25, auxiliary prime 39, q = 2 special 3, trivial 3): exact order every time.
- logs/eisen_output.txt (eisen.py, 15 s): Thm E for F_p[ω] and F_p[i], p ≤ 3·10⁴, r ≤ 7: least non-residue norms ≤ 73;
  slack ≥ 10^{15.7}.
- Self-review: every ERH/GRH/hypothesis hit in the OpenAI factoring TeX read (DERIVATION §3); line numbers rechecked by rg;
  one bug found and fixed in primary.py (r = 2, e = 0 took the q = 2 special branch).
- Reviewer (verdict only, two passes, proofs inline): Thm B and the q = 2 case SOUND; Thm E SOUND after restricting its
  non-square example to odd p (fixed); Prop. C CLOSED once the dependency evidence (paper's own statements, rg, OAINum) was
  supplied.

## Novelty
Pre-empted as above (NOVELTY.md: OpenAI preprints and CONTENTS.md, Hexagon API all 716 works, Zenodo phrase queries,
arXiv math.NT/cs.CC/cs.DS/cs.SC/math.AG listings Oct–Nov 2026, Baldwin–Sivek, FSW, 2610.37179, Pathak, HAL, our files).
Nowhere found: the Cureton–Thorner substitution in OpenAI's factoring theorem, numeric auxiliary-prime or F_{p²}
non-residue constants, the resulting exponent, or the q = 2 special case. All routine.

## Caveat
OpenAI's factoring preprint is unrefereed (our OAINum audit: finite checks, no flaw). Cureton–Thorner's theorem is
Lean-checked, its paper unrefereed; our Chebotarev theorem is unrefereed (reviewer: correct). Huang, Rónyai, Evdokimov
and Buchmann–Shoup were read only through IKS, IKRS and OpenAI's history section, so claim 6 rests on IKS's own
statement of their GRH use. The constants are far from the truth (true auxiliary primes ≤ 20 q² log p in our range).

RH relevance: N11

## Files
REPORT.md, DERIVATION.md, NOVELTY.md, nr.c, aux.py, table.py, primary.py, eisen.py, novelty/ (search scripts), logs/ (nr_output.txt,
aux_output.txt, table_output.txt, primary_output.txt, eisen_output.txt). Run: `cc -O2 -o /tmp/nr nr.c -lm && /tmp/nr 10000000`;
`uv run --no-project --with sympy python3 aux.py`; `python3 table.py`; `python3 primary.py`; `python3 eisen.py`.
Suggested COROLLARIES line: "OpenAI's deterministic polynomial-time factorization over F_p (4 Oct 2026) rests on
Cureton–Thorner's Lean-checked 31/32 alone; auxiliary primes ≤ (1.96 q² log(12q²p))^{32}; bit cost
O(((n + 1) log p)^{500000}) (stated 10^{12}); small non-cube a + bω of F_{p²} with a² − ab + b² ≤ (1.96 log 3p)^{32}.
Routine; DetFactor."
