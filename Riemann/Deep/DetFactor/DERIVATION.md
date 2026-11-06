# DetFactor: derivations (node N11)

Agent detfactor, 6 Nov 2026. Line numbers refer to the TeX sources named in NOVELTY.md (downloaded, read, deleted).

## 0. Inputs

- **(H)** Bebel (Lean-kernel-checked, QuasiRiemannTracker, 1 Nov 2026): ζ and every Dirichlet L-function are zero-free on
  Re s > B = 0.8733819666132076. Through our explicit Linnik theorem ([ExplProg](../ExplProg/REPORT.md), Theorem 3,
  [page](../../../discoveries/linnik-numeric.md)): every nonprincipal character χ mod m ≥ 3 has χ(n) ∉ {0, 1} for a
  prime n ≤ (log m)^{7.8978}.
- **(CT)** Cureton–Thorner, Hexagon 2610.00200v2, Theorem thm:dedekind (l. 171–172): for every number field K,
  ζ_K(s) ≠ 0 for Re s > 31/32, s ≠ 1; formalized in Lean (l. 1569–1576, github.com/AxiomMath/DedekindZetaZeroFree;
  our R399 checked scope and absence of sorry).
- **(Cheb)** Our numeric Chebotarev theorem ([ChebExp](../ChebExp/REPORT.md), Theorem A,
  [page](../../../discoveries/chebotarev-numeric.md)), which uses only (CT): for every Galois L/K, L ≠ ℚ, and every
  conjugacy class C of Gal(L/K) there is a degree-one prime 𝔭 of K, unramified in L, with Frob_𝔭 = C and
  N𝔭 ≤ (0.49 log 3Δ_L)^{32}.

## 1. Square roots and r-th roots in F_p (task 1): KNOWN in polynomial form; explicit exponent routine

**Theorem A.** Let p be an odd prime and r a prime dividing p − 1.
(a) The least r-th power non-residue n_r(p) of F_p is a prime ≤ (log p)^{7.8978}.
(b) Square roots: given p and a, decide by Euler's criterion; if a is a nonzero square, test n = 2, 3, … by the Jacobi
symbol until (n/p) = −1, then run Tonelli–Shanks with this n. Deterministic, O((log p)^{8.8978} log log p) bit
operations with schoolbook arithmetic (search ≤ (log p)^{7.8978} candidates at O(log p · log n) each: p mod n, then a
Jacobi symbol of small numbers; Tonelli–Shanks O(log² p) multiplications, O(log⁴ p)). The stopping rule needs no constant.
(c) r-th roots: test n = 2, 3, … by n^{(p−1)/r} ≠ 1 (≤ (log p)^{7.8978} tests, O(log³ p) each), then
Adleman–Manders–Miller with this non-residue: s = v_r(p − 1) rounds, each one exponentiation and a search over r
exponents, O((r + log p) log³ p). Total O((log p)^{10.8978} + r log³ p) bit operations, polynomial in r and log p.

*Proof.* (a) F_p^× is cyclic, so a character χ of exact order r exists; its kernel is the subgroup of r-th powers.
(H)/ExplProg Thm 3 with m = p gives a prime n ≤ (log p)^{7.8978} with χ(n) ∉ {0, 1}: n ≢ 0 and n is not an r-th power.
The least non-residue is prime, because a product of r-th powers is an r-th power. (b), (c): standard (Tonelli–Shanks,
AMM 1977); deterministic once a non-residue is supplied; costs as stated. ∎

Known in polynomial form: OpenAI QRH preprint (30 Sep 2026), Corollary cor:nonresidues-square-roots (paper.tex
l. 320–362: n(p) ≤ C(log p)^{32}, C unspecified, via BIMS's weak-GRH theorem with ε = 7/16; deterministic square roots by
Tonelli–Shanks); Hexagon 2610.00178 (11 Oct), Section sec:roots (l. 10430–10470): all roots of X^n − a in F_p in
deterministic time polynomial in n and log p. New here only: the exponent 7.8978 with constant 1 in (a) (already
ExplProg Thm 3) and the resulting explicit running-time exponents. Routine.

*Check* (logs/nr_output.txt, nr.c, 5 s): for all primes p ≤ 10⁷: max n_2 = 53 (p = 9,257,329), max n_3 = 29
(p = 4,775,569, over p ≡ 1 mod 3), max n_5 = 19; the largest log n_r(p)/log((log p)^{7.8978}) for p ≥ 17 is 0.187.

## 2. Auxiliary primes with numbers (the only analytic input of OpenAI's factoring theorem)

OpenAI, *Deterministic Polynomial Factorization over Prime Fields* (4 Oct 2026), eq. (auxiliary)
(00-introduction.tex l. 49–57): for primes q ≤ n, ℓ is an *auxiliary prime for (p, q)* if ℓ ∉ {2, 3, p, q},
ℓ ≡ 1 (mod 12q) and p^{(ℓ−1)/q} ≢ 1 (mod ℓ).

**Theorem B.** Let q and p be primes with p ∤ 6q. Put m = 12q and
D(q, p) = q φ(m) log m + φ(m) (q log q + (q − 1) log p).
Then the least auxiliary prime ℓ_q(p) satisfies
ℓ_q(p) ≤ (0.49 (log 3 + D(q, p)))^{32} ≤ (1.96 q² log(12 q² p))^{32}.
Only (CT) is used (through (Cheb)). With Liu's Hecke half-plane B_2 (Hexagon 2611.00082, unrefereed) chebexp's part (B)
gives (0.73 (log 3 + D))^{7.8978} instead.

*Proof.* Let F = ℚ(μ_m), α = p^{1/q}, L = F(α): the splitting field of (X^m − 1)(X^q − p) over ℚ, so L/ℚ is Galois.
(i) [L : F] = q: p ∤ m, so p is unramified in F; at a prime of F above p, X^q − p is Eisenstein, hence irreducible.
(ii) N = Gal(L/F) ≅ μ_q by σ ↦ σ(α)/α. For τ ∈ Gal(L/ℚ) with τ(ζ_m) = ζ_m^a one computes τ σ_ζ τ^{−1} = σ_{ζ^a}; since
(ℤ/m)^× → (ℤ/q)^× is onto, C = N \ {1} is one conjugacy class of Gal(L/ℚ) (for q = 2, N is central and C = {σ}).
(iii) Each of 2, 3, q, p ramifies in L (2 in ℚ(i), 3 in ℚ(√−3), q in ℚ(ζ_q) or q = 2, p totally in ℚ(α)).
(iv) (Cheb) with K = ℚ gives a prime ℓ unramified in L with Frob_ℓ ∈ C and ℓ ≤ (0.49 log 3Δ_L)^{32}. Unramified:
ℓ ∉ {2, 3, q, p}. Frob_ℓ|_F = 1: ℓ splits completely in F, so ℓ ≡ 1 (mod 12q). Frob_ℓ ≠ 1: a prime 𝔩 | ℓ of F
(residue field F_ℓ) does not split completely in L. As 𝔩 ∤ disc(X^q − p) = ±q^q p^{q−1}, the order O_F[α] is
𝔩-maximal and Dedekind–Kummer applies: X^q − p is not a product of distinct linear factors mod 𝔩. Since μ_q ⊂ F_ℓ, a
single root would split it completely; so p is not a q-th power mod ℓ, i.e. p^{(ℓ−1)/q} ≢ 1 (mod ℓ).
(v) Δ_L = Δ_F^{q} · N_{F/ℚ}(𝔡_{L/F}), 𝔡_{L/F} divides disc(X^q − p)·O_F, so N(𝔡_{L/F}) ≤ (q^q p^{q−1})^{φ(m)}; and
Δ_F = m^{φ(m)} / ∏_{ℓ|m} ℓ^{φ(m)/(ℓ−1)} ≤ m^{φ(m)}. Hence log Δ_L ≤ D(q, p); (0.49 log 3Δ)^{32} increases with Δ.
(vi) Simplification: φ(12q) = 4(q − 1) for q ≥ 5 gives D ≤ 4q(q − 1) log(12q²p) and
0.49 (log 3 + D) ≤ 1.96 q² log(12q²p); q = 2 (φ(24) = 8) and q = 3 (φ(36) = 12) directly (D = 8 log(2304p),
resp. 36 log 108 + 24 log p). ∎

*Check* (logs/aux_output.txt, aux.py, 4 s): all primes 5 ≤ p ≤ 2·10⁵ and q ∈ {2, 3, 5, 7, 11, 13}: the least
auxiliary prime is ≤ 2113 and ≤ 20 q² log p; the bound exceeds it by a factor ≥ 10^{48}. Adversarial primes
p ≡ 1 mod (product of the first k primes ≡ 1 mod 12q), up to 10^{37.8}: ℓ_q(p) grows like the k-th such prime
(≤ 2521), far below the bound.

## 3. OpenAI's factoring theorem on a Lean-checked analytic input (routine input substitution)

**Proposition C.** In OpenAI's factoring preprint the analytic theorem (Theorem analytic:hecke, 10-analytic.tex
l. 9–18: Hecke L-functions over cyclotomic fields ⊇ μ_12 zero-free on Re s > 1 − 10^{−6}; proof in the companion
*Primitive roots for every admissible integer base*, Theorem 1.2, not Lean-checked) is used only through Lemma
analytic:explicit (l. 39–48) (dependency audit: the scope check below), which needs a zero-free strip for the Dedekind zeta functions of K = ℚ(μ_{12q}) and
M = K(p^{1/q}), inside Proposition analytic:auxiliary (l. 91–98). (CT) gives that strip with δ = 1/32 for every number
field, and Theorem B gives the auxiliary primes outright. Hence Theorem thm:factorization holds with (CT) in place of
the companion's theorem. Inserting Theorem B into Proposition complexity:bound (09-complexity.tex l. 19–30: bit cost
O(B^{500000} + B^{100} E^{20}), B = 20 + (n + 1)(L + 1), L = ⌈log₂ p⌉, E = largest auxiliary prime): auxiliary primes
are needed only when p > B^{200000}; then log n ≤ L/200000, so 1.96 q² log(12q²p) ≤ 1.96 n² (2.485 + 2 log n + 0.6932 L)
≤ n²(L + 1)² < B², E < B^{64}, B^{100}E^{20} < B^{1380}, and the total is O(B^{500000}) = O(((n + 1)⌈log₂ p⌉)^{500000})
bit operations, instead of the stated O(((n + 1)L)^{10^{12}}) (their derivation: E ≤ c₀B^{2·10^7}, c₀ unspecified,
giving B^{4·10^8+100}).

*Scope check.* rg -i 'ERH|GRH|Riemann hyp|hypothesis|Weil|Chebotarev|Linnik|Burgess|density|zero-free|L-function|analytic'
over sections 02–09: the only hits are 'characteristic hypotheses' (02 l. 174: p > h²; 06 l. 92; 07 l. 182, 209:
p > (N + 4g + 1)²), all elementary. GRH appears only in the introduction (l. 37, 98–106) and in Section 10. The paper says
the same: "The reduction in Theorem thm:reduction does not use that theorem" (00-introduction l. 79–80) and "all
dependence on that paper is isolated in this input" (10-analytic l. 6–7); our OAINum audit read all eleven sections.
So Section 10 is the only analytic input. (Second reviewer pass on Prop. C with this evidence: CLOSED.) The algebraic sections 02–09 (≈ 140 KB: cyclic curves Y^q = F(X), divisor-class division,
norm equations via maximal orders) were audited in our OAINum round ([report](../OAINum/REPORT.md), R255, 5 Nov: all
sections read, finitely checkable lemmas retested, no flaw; complexity bookkeeping not run); not re-audited here.

*Executable check of the table section* (01-table.tex; table.py, logs/table_output.txt; OAINum's test_table.py covered
p ≤ 103): for 10 primes p (7 to
10^{18} + 9, both classes mod 3): U_q (Gauss periods of ℓ_q) is a field of degree q exactly when p is not a q-th power
mod ℓ (η_0^p ≠ η_0, η_0^{p³} = η_0; split in the control case); for q = 2, c = v² is a non-square and equals ℓ_2
(quadratic reciprocity); for q = 3, ω_3 lies in K_3 = F_p(ζ_3) and has order exactly 3^{v_3(|K_3| − 1)}, including
K_3 = F_{p²} for p ≡ 2 (mod 3) (e.g. p = 10⁹ + 7: ℓ_3 = 181, ω_3 of order 9).

## 4. Non-residues in F_{p^k} (task 2)

**Lemma N.** For a finite field F_Q, k ≥ 1 and a prime r | Q − 1, x ∈ F_{Q^k}^× is an r-th power iff
N_{F_{Q^k}/F_Q}(x) is an r-th power in F_Q. *Proof:* x^{(Q^k−1)/r} = (x^{(Q^k−1)/(Q−1)})^{(Q−1)/r} = N(x)^{(Q−1)/r}. ∎

**Corollary N1.** If r | p − 1 and r ∤ k, the least prime that is not an r-th power in F_{p^k} is ≤ (log p)^{7.8978}
(N(n) = n^k for n ∈ F_p, and n^k is a non-residue iff n is, as r ∤ k; then Theorem A(a)). If r | p − 1 and r | k, or
r ∤ p − 1, no element of F_p is an r-th non-residue of F_{p^k}; with d = ord_r(p), Lemma N reduces the problem to
F_{p^d} = F_p(ζ_r) when v_r(k/d) = 0.

**Algorithmic form: KNOWN.** OpenAI's table (01-table.tex l. 45–193) builds K_q = F_p(ζ_q) (one call to factor Φ_q,
degree q − 1) and a generator ω_q of the q-primary subgroup of K_q^× — a q-th power non-residue of F_p(ζ_q) — in
O(ℓ_q^{10} B^{40}) bit operations plus that call (09 l. 244–268); with Theorem B the search bound is explicit. For any
explicit field k ⊇ K_q (q odd), OpenAI's Lemma ff:unitroot (02-finite-fields.tex l. 106–) gets the q-primary generator of
every component by its eigenvector branch (solve v^{Q^{d/q}} = ζ_q v when q | d; else ω_q, by v_q(Q^d − 1) = v_q(Q − 1) +
v_q(d)) and then extracts promised q-th roots in k[X]/H. Factoring over explicitly given F_{p^k} reduces to F_p
(Kriebel, Zenodo 23205034).

*Check of the composite-field construction, including q = 2* (primary.py, logs/primary_output.txt, 70 random fields
k = F_p[t]/(g), p from 1009 to 1,000,033, r ∈ {2, 3, 5, 7}, [k : F_p] ≤ 14, ζ_r given): the output has order exactly
r^{v_r(|k|−1)} in all 70 (25 by the eigenvector branch, 39 by the auxiliary-prime base step with U_r ⊗ k, 3 trivial,
3 by the case below). The one case outside OpenAI's lemmas: r = 2, p ≡ 3 (mod 4), v_2([k : F_p]) = 1, where
v_2(p² − 1) = 1 + v_2(p + 1) exceeds what the eigenvector step gives. Then: a nonzero u ∈ F_{p²} with u^p = −u (a kernel vector), i = u/c with
c² = −u² in F_p (so i² = −1), the least prime ℓ ≡ 1 (mod 4) with (ℓ/p) = −1 written as ℓ = a² + b², and x = a + bi has
norm ℓ, a non-square of F_p, so x is a non-square of F_{p²} (Lemma N); x^{(|k|−1)/2^t}, t = v_2(|k| − 1), generates the 2-Sylow subgroup ([k : F_{p²}] is odd, so x stays a non-square of k).
By our Chebotarev theorem for L = ℚ(i, √p) (Δ_L = 16p², class {σ: σ(i) = i, σ(√p) = −√p}),
ℓ ≤ (0.49 (log 48 + 2 log p))^{32}. (Routine.)

**Theorem E (small non-residues of F_{p²} in the bases {1, ω} and {1, i}; routine).** Let K be imaginary quadratic with
class number h and discriminant Δ_K, p a prime inert in K (O_K/p ≅ F_{p²}), r a prime dividing p² − 1, and suppose every
unit of K is an r-th power in O_K/p. Then some α ∈ O_K \ pO_K is not an r-th power mod p and
N(α) ≤ (0.49 (log 3 + h r (log|Δ_K| + 2 log p)))^{32}; for h = 1, N(α) ≤ (0.49 (log 3 + r log|Δ_K| + 2(r − 1) log p))^{32}.
In particular for odd p ≡ 2 (mod 3) and F_{p²} = F_p[ω]/(ω² + ω + 1): a non-cube a + bω with a² − ab + b² ≤ (1.96 log 3p)^{32}
(if ω itself is a cube, i.e. 9 | p² − 1; otherwise ω is the non-cube), so |a|, |b| ≤ √2 (1.96 log 3p)^{16}; a non-square
with norm ≤ (0.49 (3 log 3 + 2 log p))^{32}. For p ≡ 3 (mod 4) and F_p[i]: norm ≤ (0.49 (log 3 + r log 4 + 2(r − 1) log p))^{32}.

*Proof.* Cl_p = ray class group of K mod pO_K; S = image of (O_K/p)^×, so Cl_p/S = Cl_K and S ≅ (O_K/p)^×/im(O_K^×).
ψ(x) = x^{(p²−1)/r} is trivial on the units, so it induces ψ̄ on S onto μ_r; U = ker ψ̄. The class field L of Cl_p/U has
[L : K] = hr, conductor dividing p; its characters trivial on S/U (h of them) are unramified, the others have conductor
norm ≤ p², so (conductor–discriminant) log Δ_L ≤ hr log|Δ_K| + 2h(r − 1) log p. Take σ ∈ Gal(L/K) ≅ Cl_p/U a nontrivial
element of S/U ({σ} is a conjugacy class: abelian). Our Chebotarev theorem (base K) gives a degree-one prime 𝔮 of K,
unramified in L (so 𝔮 ∤ p), with Artin symbol σ and N𝔮 ≤ (0.49 log 3Δ_L)^{32}. Its class lies in S, so 𝔮 = (α) is
principal, and the class of (α) is the image of α mod p; it is ≠ 1 modulo U, so ψ(α mod p) ≠ 1. For {1, ω},
a² − ab + b² ≥ (a² + b²)/2. Units of ℚ(√−3): ±ω^j, all r-th powers unless r = 3 and 9 ∤ p² − 1; units of ℚ(i): always. ∎

Precursor: Hexagon 2610.00178, Lemma lem:eisenstein-cubic-nonresidue (l. 10745–10815): r = 3, K = ℚ(√−3),
x = C (log 3p)^{32} with C effective but not computed; Theorem E gives C = 1.96^{32} and every r.
*Check* (eisen.py, logs/eisen_output.txt, 15 s): all p ≤ 3·10⁴ in the class, r ∈ {2, 3, 5, 7}: least non-residue norms
≤ 73 ({1, ω}) and ≤ 61 ({1, i}); the bound exceeds them by a factor ≥ 10^{15.7}.

**Small elements in a fixed polynomial basis for k ≥ 3: NO RESULT.** Blocking step: for F_{p^k} = F_p[t]/(g) and
K = ℚ(θ) of degree k with p inert, the r-th power residue character of O_K/pO_K extends to a ray class (Hecke) character
only on the part trivial on the image of the global units, and (CT)/(Cheb) produce small prime *ideals* 𝔮 with
χ(𝔮) ≠ 1, not small *elements* a_0 + a_1θ + …: a generator of 𝔮 is small only up to units, and for k ≥ 3 the unit group
is infinite (its covering radius, hence the generator size, is not polylogarithmic in general). For imaginary quadratic
K (k = 2) the unit group is finite and Theorem E goes through. The dual (Kummer) form — small primes ℓ whose Gauss
periods give the non-residue — is Theorem B, which is what OpenAI's table uses.

## 5. ERH steps, algorithm by algorithm (task 3)

| Algorithm | ERH step (faithful source, line) | Replaced by | Status |
|---|---|---|---|
| Tonelli–Shanks, Cipolla | a quadratic non-residue mod p (Ankeny/Bach) | Thm A(a) | KNOWN: OpenAI QRH Cor. (l. 320–362) |
| Adleman–Manders–Miller 1977 | r-th power non-residue mod p | Thm A(a) | KNOWN: Hexagon 2610.00178 sec:roots (l. 10430–10470) |
| Huang 1985/1991 (abelian Galois group) | primitive r-th non-residues in finite fields (IKS 0804.1974 l. 94–97 for the common mechanism; Huang named as special case of Rónyai 1992, l. 100–105) | OpenAI table + Thm B | subsumed by OpenAI's theorem; original not read (paywalled) |
| Rónyai 1988 (prime r ∣ deg f; time poly(n^r, log q)) | primitive r-th non-residues in a field k ⊇ F_q(ζ_r), used to solve x^r = a (IKS l. 94–97, 111–121; IKRS 0811.3165 l. 150–165) | r = 2, 3: Hexagon 2610.00178 Thms even-splitting (l. 10605), three-divisible-splitting (l. 10735); all r: OpenAI table ω_r + Thm B | KNOWN (r = 2, 3); subsumed (all r) |
| Rónyai 1992 (Galois groups) | decomposing algebras with explicit automorphisms under GRH (IKS l. 485: "Theorem 2.3 in [ro4]") — same non-residue mechanism | as above | subsumed |
| Evdokimov 1994 ((n^{log n} log q)^{O(1)}) | as Rónyai 1988, in algebras of dimension n^{log n} (IKS l. 124–128; IKRS l. 165–167) | as above | unconditional by the classical route below; subsumed by OpenAI (claimed polynomial) |
| Ivanyos–Karpinski–Saxena 2009 (m-schemes) | "needed only to find primitive r-th nonresidues in a finite field F_q" (IKS l. 94–97); the field k ⊇ F_p with s-th non-residues for all s ≤ m, "poly(log p, m^m) time under GRH" (l. 465–466) | OpenAI table + Lemma ff:unitroot + Thm B (+ q = 2 case) | unconditional (classical route below) |
| BIMS 2017 (r-th roots) | Weak GRH, Conj. 6.3 (l. 732–735): a Dirichlet strip of half-width ε; n(p, r) = O(log^{4/(1−2ε)} p) (Thm thm-npr, l. 784–787) | (H): ε = B − ½ gives exponent 15.8 in their analysis; Thm A(a) gives 7.8978 with constant 1 | KNOWN (OpenAI used ε = 7/16, A = 32) |
| OpenAI 2026 (poly time, every degree) | none; uniform Hecke strip of width 10^{−6} (companion Thm 1.2) only in Section 10 | (CT) + Thm B (Prop. C) | claimed, unrefereed; analytic input now Lean-checked |

*Classical route without OpenAI's splitter.* The subroutine that IKS name as the only use of GRH — primitive r-th
non-residues of the working field k ⊇ F_p(ζ_r) and roots of x^r − a in k-algebras (IKS l. 94–97 say this is the only
use of GRH in their paper, which includes the algebra decompositions of l. 485, 548, 600 via Rónyai 1992's Theorem 2.3;
the field k: l. 465–466) — is supplied without GRH by OpenAI's
elementary table (01-table) and Lemma ff:unitroot (02 l. 106–), with auxiliary primes from Theorem B (Cureton–Thorner
only), plus the q = 2 exceptional case above. The table's one call to factor Φ_q (degree q − 1) is made with the classical
algorithm itself, by induction on q using only entries for primes < q (the same induction as OpenAI's Prop.
table:construction), and the working field k is the compositum of the K_s so obtained. These parts do not use OpenAI's
odd-degree splitter (sections 04–07); our OAINum audit and primary.py checked them. Hence, taking IKS's description of the GRH use as faithful, Rónyai's
poly(n^r, log p) splitting (prime r | n) and Evdokimov's (n^{log n} log p)^{O(1)} complete factorization hold
unconditionally. Inputs with p at most the square of the working dimension (ff:unitroot needs p > h²) go to Berlekamp's
deterministic algorithm, polynomial in p and hence within the same bound. Routine; superseded by OpenAI's claimed
polynomial bound.

**Strongest unconditional deterministic factoring theorem** (claimed, unrefereed): OpenAI's — complete factorization
over F_p in O(((n + 1)⌈log₂ p⌉)^{500000}) bit operations with (CT) as the only analytic input (Prop. C); over explicitly
given F_{p^k} via Zenodo 23205034's reductions.
