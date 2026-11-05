NO IDEA: homological stability supplies only RH-free topology (stable ranges, stable homology, Betti bounds); it becomes arithmetic only through Frobenius (Grothendieck–Lefschetz on moduli of curves) for the main term and Deligne's Weil II on every unstable degree for the error, with Weil's RH for each curve built in, and its only shadow over Q is the CFKRS/DGH moment recipe, which L(1/2, χ_8d)^k obeys numerically (k ≤ 4, d ≤ 10^6) and which a Davenport–Heilbronn function obeys too, so it carries no zero information.

# Exploration explore9: homological stability, Hurwitz spaces and the scanning map → RH tree (2026-11-05)

Field: homological stability of configuration and Hurwitz spaces and its consequences over F_q(t). Brief: .scratch/orca/explore9-brief.md. Derivations, TeX locations and checks: [derivations.md](derivations.md). Notation: n = deg d = 2g + 1; θ = uniform stable range; Θ̄_d = q^{−1/2}·Frobenius of y² = d(t), a conjugacy class in USp(2g); DH = Davenport–Heilbronn.

## Survey (one page)
- **EVW** (0912.0325, Annals 2016; II: 1212.0923). Stability for dihedral Hurwitz spaces. Cohen–Lenstra densities over F_q(t) are reached as q → ∞. The point count |X_n(F_q)| = Σ_j (−1)^j Tr(Frob | H^j_c) (Grothendieck–Lefschetz). With dim H^i ≤ C^i and Deligne's bound q^{n−i/2} on H^{2n−i}_c, |#X_n(F_q)/q^n − 1| ≤ 1/(√q/C − 1) for q > C² (EVW §1.6, §8). The main term counts the Frobenius-fixed components (EVW II, lifting invariant).
- **Landesman–Levy** (2410.22210, 2410.22222, 2503.03861, 2510.02068). Integral stability for every finite rack. The stable homology of each component equals that of Conf_n (after inverting finitely many primes). Consequences, for fixed large q: Cohen–Lenstra moments, Cohen–Lenstra–Martinet, Malle (leading constant: Santens 2512.12838), Bhargava's S_d counts, Selmer moments. Point counts use Grothendieck–Lefschetz, Sun's Deligne bounds for stacks and an exponential cohomology bound (Lemma "weak-point-bound").
- **Bergström–Diaconu–Petersen–Westerland** (2302.07664).
  - Thm A (topology, via the scanning map and Koszul duality): a plethystic formula for the stable homology of β_n with coefficients in S^λ(V), V the integral Burau representation.
  - Thm B: the stable homology is pure Tate of weight −2i + k. It follows from a topological "indecomposables" theorem plus the Galois action on degree-0 classes; Weil II is not used.
  - The stable Frobenius traces reproduce exactly the CFKRS arithmetic factor R = Π_P (1 + (1 + |P|^{−1})^{−1} Σ_{j>0} h_{2j}).
  - Thm C (= 11.3.19): moments = Q_r(2g+1) + O(4^{g(r+1)} q^{−θ(2g+1)/2}). The error uses Weil II on the unstable degrees and Fuks cells (dim H_k ≤ binom(n−1, k) dim V).
- **Miller–Patzt–Petersen–Randal-Williams** (2402.00354): θ(n) = (n − 35)/34 holds uniformly in λ. Hence CFKRS holds for every moment over F_q(t) once q > 4^{34(r+1)}, already explicit (BDPW intro). Conjectured slope 1/4.
- **Ellenberg–Landesman** (2310.16286): Hurwitz spaces over punctured bases of any genus; BKLPR heuristics for Selmer groups of quadratic twists, with an error that tends to 0 as q grows.
- **Wang** (2402.01214): the ratios conjecture over F_q(t) for large fixed q, at distance ≥ q^{−δ} from the line ("cancellation beyond GRH" on average over the family).
- **Secondary terms.**
  - Diaconu–Goldfeld–Hoffstein conjecture D^{1/2+1/(2n)} terms for r ≥ 4.
  - r = 3: the D^{3/4} term is proved over F_q(t) (Diaconu 1801.00486) and over Q (Diaconu–Whitehead 1804.00690, residue −0.0034).
  - BDPW §13: the first term T_1 is the stable homology; T_n for n ≥ 2 are conjecturally secondary (higher) stability; tr_λ ~ T_2 unconditionally for λ = (3^g).
  - Jung 2608.22717, 2609.12764: Florea's q^{2g/3} first-moment term has period-3 coefficients.
- **2025–26 follow-ups** add counts, not zeros. Ellenberg–Shusterman 2601.01242 average over conductors. Sawin–Wood 2602.21032 study unramified extensions and give a number-field conjecture with non-existence results. Lalín–Lee–Oliver–Pozdnyakov 2608.01337 identify the murmuration density with a stable-homology correction (already in R296).
- **Tree memory.** R296's "lift to twisted local systems for uniform family moment bounds (idea)" is exactly BDPW + MPPRW (2023–24): known. R320 already noted "Hurwitz-space Cohen–Lenstra carries no zero information".

## (i) Where Frobenius, Grothendieck–Lefschetz and Deligne enter

| step (BDPW / EVW / LL) | input | RH-free? | counterpart over Q |
|---|---|---|---|
| 1. homological stability, uniform θ (RW–Wahl, MPPRW, EVW, LL) | topology of β_n, Hurwitz spaces over C | yes | none: no space whose points are the quadratic fields |
| 2. stable homology (scanning, Koszul duality; LL: = Conf_n) | topology | yes | shadow only: the CFKRS arithmetic factor, an Euler product (as in Farb–Wolfson–Wood "homological = arithmetic densities") |
| 3. Galois action on stable classes (Thm B; EVW II components) | Frobenius on étale homology, Galois-equivariant E_2 operad | uses Frobenius, not Weil II | none (R35: no Frobenius over Q) |
| 4. Grothendieck–Lefschetz on H_{g,1,0} (Behrend for stacks) | Frobenius | — | none; family averages must be done by Poisson summation / approximate functional equations instead |
| 5. Howe duality: q^{−n} Σ_d L(1/2,χ_d)^r = Σ_{λ ⊆ (r^g)} tr_λ(g) dim V_{λ†} | algebra on the eigenvalues of Θ̄_d | algebra yes, but Θ̄_d ∈ USp(2g) is Weil's RH for each curve | a unitary "Frobenius" with the zeros as eigenvalues = Hilbert–Pólya (N13): circular |
| 6. error on degrees k > θ: Weil II (Frobenius eigenvalues of absolute value ≤ q^{−k/2}) × Fuks-cell Betti bound; V_λ pure of weight equal to the size of λ | Deligne [Del80] + Weil's RH for each fiber | no | none; the analytic stand-ins (Pólya–Vinogradov, Heath-Brown's quadratic large sieve) give the full polynomial only for r ≤ 3 unconditionally (r = 4: Shen under GRH; Shen–Stucky unconditionally with four main terms) |
| 7. q > 4^{34(r+1)} | Betti growth 4^{g(r+1)} against the saving q^{−θ/2} | — | Q is the fixed-q regime: no large parameter |

The RH-free steps (1, 2, the algebra of 5) never act on zeros. Their only number-field counterpart is the moment recipe.

## (ii) Numeric test over Q: does the stable-homology main term give the right lower-order terms?

The family is χ_8d with d odd and squarefree, d ≤ 10^6 (405,286 discriminants).
- L(1/2) comes from the approximate functional equation (lvals.py, 198 s). It matches mpmath.dirichlet to 5·10^{−13}.
- The full CFKRS polynomial Q_k(L), a polynomial in L = log(8d/π), comes from an FFT Taylor extraction (cfkrs.py). Three checks pass:
  - the direct recipe sum at small shifts agrees (k = 2: 3·10^{−7}, k = 3: 5·10^{−6} relative);
  - k = 1: the leading coefficient is A(0)/2 with A(0) = ½·(6/5)·0.70444220 (carefree constant);
  - the Diaconu–Whitehead residue comes out at −0.003420 (theirs: ≈ −0.0034).
- Q_1 = 0.21133266 L + 0.41727854.
- Q_2 = 0.0028577887 L³ + 0.050189885 L² + 0.16557888 L + 0.041810766.
- Q_3 = 1.1183240e−6 L⁶ + 7.4002891e−5 L⁵ + 1.6752752e−3 L⁴ + 0.015011480 L³ + 0.048693395 L² + 0.033843654 L + 0.018334688.
- Q_4 (degree 10, leading coefficient 6.9159159e−12, i.e. a_4 = 3.346·10^{−5}): all coefficients are in Q4.json. Over F_q(t) this is BDPW + MPPRW's theorem; over Q the full polynomial is only conjectural (GRH: Shen; four main terms unconditionally: Shen–Stucky).

| k | data/recipe − 1, sharp, X = 10^5 | same, X = 10^6 | smoothed W(d/X), X = 10^6 | (data − recipe)/√X, X = 10^6 | leading term only / recipe − 1, X = 10^6 |
|---|---|---|---|---|---|
| 1 | +1.1·10^{−5} | −2.9·10^{−5} | −1.2·10^{−5} | −0.04 | −12.6% |
| 2 | +1.7·10^{−4} | −1.1·10^{−4} | −4.9·10^{−5} | −0.88 | −61.1% |
| 3 | +6.9·10^{−4} | −2.9·10^{−4} | −4.2·10^{−5} | −18.7 (smoothed −0.96) | −94.9% |
| 4 (conjectural over Q) | +1.7·10^{−3} | −5.7·10^{−4} | +1.2·10^{−4} | −358 (smoothed +27.6) | −99.87% |

Readings:
- The lower-order terms, which over F_q(t) are the stable classes with λ ≠ 0, carry most of the moment at this height: the leading (random-matrix) term alone misses 13%, 61%, 95% and 99.9% for k = 1, 2, 3, 4.
- With them the relative agreement is 10^{−5} to 6·10^{−4} at X = 10^6, k = 4 included (over Q a conjecture). The residual over √X grows with k only because the k-th powers are large.
- The stable-homology secondary term (Diaconu–Goldfeld–Hoffstein's X^{3/4}, residue −0.00342, Ŵ(3/4) = 0.325) predicts −35 at X = 10^6, against a smoothed residual of −962. It becomes visible only once 0.00111·X^{3/4} exceeds the residual ≈ √X, i.e. X ≈ 10^{11}–10^{12}. Alderson–Rubinstein (1110.0253, L(1/2, χ_d) for −5·10^{10} < d < 1.3·10^{10}) reported evidence for it at the edge of that range.
- Verdict: everything matches. For k ≤ 3 this is a known theorem (k = 1 Jutila and Young, k = 2 Jutila and Sono, k = 3 Soundararajan, Young and Diaconu–Whitehead). For k = 4 it is the CFKRS conjecture, which Alderson–Rubinstein already tested numerically to |d| ≈ 5·10^{10}. There is no new information and no zero information. Output: compare_out.txt.

## Transfers (both adversarial tests)

**T1 (N9). Stable homology = CFKRS main term, so the recipe over Q.**
- Test (ii) passes.
- Frobenius test: over F_q(t) the main term needs steps 3 and 4. Over Q the same polynomial comes from the approximate functional equation and the orthogonality of χ_d(n). The topological content of Thm A collapses to the Euler product of the arithmetic factor; nothing beyond CFKRS survives.
- DH test: fails. The DH function f = ((1−iκ)/2)L(s,ψ) + ((1+iκ)/2)L(s,ψ̄) (ψ mod 5) vanishes at 0.808517182457 + 85.699348485378 i (|f| = 6.8·10^{−13}), yet ∫_0^T |f(1/2+it)|² dt follows its own recipe ((1+κ²)/2)(4/5)T(log(5T/2π) + 2γ − 1 + ½ log 5) + 0.395675·T: −1.2%, −0.30%, +0.06%, +0.32%, +0.04% at T = 100, 200, 400, 800, 1200 (dh_test_out.txt). The recipe for L(s,ψ) alone (zeros on the line as far as computed) is off by −0.75%, +0.24%, +0.45% at T = 400, 800, 1200: no better than DH. Moment main terms see variances, not zeros (NOGO R8.2).
- Verdict: dead (known, blind to DH).

**T2 (N9). Unstable and secondary stability = the DGH secondary terms D^{1/2+1/(2n)}.**
- Frobenius test: over F_q(t) unstable classes are controlled only by Weil II. Over Q the terms are poles of multiple Dirichlet series, continued by functional equations alone; this is RH-free, and the continuation reaches only Re w > 2/3 for r = 3.
- DH test: fails. DH-type twist families (L(s, ψ^{±}χ_d) combinations) have multiple Dirichlet series with the same group of functional equations.
- Numerics: the D^{3/4} term is invisible at d ≤ 10^6 (predicted −35 against a residual of −962) and below about 10^{11}.
- A structural remark (standard, made explicit here). Over Q, in the known continuations, the zeros of ζ enter the quadratic-family moments only through ζ(2w)-type factors (the squarefree sieve): Σ_a μ(a) a^{−2w} has poles at w = ρ/2, giving terms of size X^{β/2} ≤ X^{1/2}. Over F_q(t) they are absent, because 1/ζ_{F_q[t]}(2w) = 1 − q^{1−2w} is entire (a genus-0 base has no zeros). Seeing them would need moment error terms below √X, which is at least as hard as the squarefree count with error O(X^{1/4+ε}) (that count already implies RH).
- Verdict: dead (known structure, chance ≈ 0).

**T3 (N11/N13). Uniform stability for all λ ⇒ every moment, large-support n-level densities, ratios (Wang).**
- Frobenius test: this needs step 5, so the average over the family of s_⟨λ⟩(Θ̄_d) for a unitary Θ̄_d. That exists over F_q only because of Weil's RH for each member. Over Q it would be a spectral (Hilbert–Pólya) realisation of the zeros of each L(s,χ_d).
- DH test: the mechanism uses what DH lacks (unitarity), but that is RH for the members itself.
- Wang's ratios likewise need 1/L(s,χ_d) to be holomorphic for Re s > 1/2.
- Verdict: circular.

**(iii) Explicit q-threshold over F_q(t).**
- Moot: BDPW + MPPRW already give q > 4^{34(r+1)} (r = 3: 7.6·10^{81}).
- Bounding Σ_{k>θ} binom(2g, k) q^{−k/2} by binom(2g, θ) q^{−θ/2} instead of 2^{2g} q^{−θ/2}, using BDPW's own Fuks-cell lemma, gives q > 8.3·10³·4^{34r} (r = 3: 2.1·10^{65}).
- The conjectured slope 1/4 would give q > 90·256^r.
- This is bookkeeping on their lemma: routine, insignificant, not recorded.

| transfer | RH relevance | chance | gain | verdict |
|---|---|---|---|---|
| T1 stable homology → recipe | low (N9, values not zeros) | high (already true) | none | dead: known, DH-blind |
| T2 secondary stability → DGH terms; ζ-zeros via squarefree sieve | low–medium (zeros of ζ appear) | ≈ 0 (needs errors below √X) | small | dead: known structure |
| T3 uniform stability → all moments, ratios, n-level densities | medium (N11/N13) | 0 (needs RH per member) | — | circular |
| (iii) q-threshold | none (F_q(t)) | high | insignificant | not recorded |

## Side observation (routine; for COROLLARIES triage, N10)
- In the known continuations, zeros of ζ reach the quadratic-family moments over Q only through ζ(2w)-type factors (the squarefree sieve, and the normalising factor of the multiple Dirichlet series), as terms X^{β/2} ≤ X^{1/2}.
- The simplest instance is the squarefree count. Evelyn–Linfoot's argument with M(x) ≪ x^{B+ε} from the Lean-checked half-plane gives Q(X) = 6X/π² + O(X^{1/(3−B)+ε}) = O(X^{0.47023}) unconditionally. Before 2026 the best was Walfisz's X^{1/2} exp(−c(log X)^{3/5}(log log X)^{−1/5}).
- k-free numbers: O(X^{1/(k+1−B)+ε}), e.g. k = 3: X^{0.3198}.
- No explicit gain: with 1.84·10^6 x^{0.9382} it beats √X only beyond 10^{240}.
- Status: not reviewed; novelty checked by an arXiv API query only. Derivation in [derivations.md](derivations.md) §5.

## Tree note (for the orchestrator)
- N12 records to add (function-field analogues): all CFKRS moments of quadratic L-functions over F_q(t) for q > 4^{34(r+1)} (BDPW 2302.07664 + MPPRW 2402.00354); Cohen–Lenstra moments over F_q(t) for q large (Landesman–Levy 2410.22210, 2503.03861).

## Suggested NOGO line
- explore9 (Orca exploration: homological stability / Hurwitz spaces / scanning map → N9, N11, N12) NO IDEA: the topology (uniform stable range, stable homology, Fuks-cell Betti bounds) is RH-free but becomes arithmetic only through Grothendieck–Lefschetz on moduli of curves and Weil II on every unstable degree, with Weil's RH for each curve built into the unitary Frobenius; its only shadow over Q is the CFKRS/DGH recipe (full polynomials match L(1/2, χ_8d)^k, k ≤ 4, d ≤ 10^6, to ≤ 6·10⁻⁴; the leading terms alone miss 13% to 99.9%), which a Davenport–Heilbronn function obeys too; the q-threshold is already explicit (q > 4^{34(r+1)}, BDPW + MPPRW), and R296's twisted-local-system idea is that work. Filter: an error term bounded by Frobenius weights on unstable cohomology has no number-field counterpart. [ExploreHomStab](Riemann/Deep/ExploreHomStab/REPORT.md)

## Files
- lvals.py, check_lvals.py: L(1/2, χ_8d) table by the approximate functional equation, and the mpmath check.
- cfkrs.py, check_cfkrs.py: CFKRS polynomials by FFT Taylor extraction, and the direct recipe-sum check. Outputs: Q1.json … Q4.json.
- compare.py → compare_out.txt: data against the recipe, sharp and smoothed, plus the Diaconu–Whitehead term.
- dh_test.py → dh_test_out.txt: DH zero check and mean square against the recipe.
- derivations.md: TeX locations of every Weil II and Grothendieck–Lefschetz use, the recipe derivation, the checks and the threshold arithmetic.
- The 405,286-value table (.scratch/explore9/L1e6.npz) is regenerated by "lvals.py 1000000" in about 200 s and is not kept.

RH relevance: none. Shadows only: N9 (moment recipe, values not zeros), N11/N13 (needs RH per member), N12 (function-field side; the q-threshold is already explicit).
