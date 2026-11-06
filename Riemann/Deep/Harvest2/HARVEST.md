CANDIDATES: 12 (best: arXiv 2609.25393 + 2610.04678, node N11, every modulus covered by O(log Q) small primes once the Lean-checked half-plane removes the exceptional set)

Scope: arXiv monthly listings math.NT/CA/CV/SP for 2026-07, 2026-08, 2026-09 (3415 entries; 1815 math.NT with cross-lists); 2173 abstracts read through the API; 191 papers selected by node keywords and read in TeX (`.scratch/harvest2/dl_scan.py`: the brief's phrase list plus 14 similar phrases; 164 had a match, 89 a strong one: 227 strong lines, 990 "sufficiently large"/"absolute constant" lines; hits kept in `.scratch/harvest2/hits/`). Zenodo API, quoted phrases only, records 2026-07-01 to 2026-11-06: 25 phrase queries, 662 records, 35 on-tree records whose full descriptions were read (all unrefereed). Each hit checked against TREE, COROLLARIES, DISCOVERIES, NOGO, REFUTATIONS, FORMALIZATIONS, discoveries/, Riemann/Tree and Riemann/Deep. Harvest agent harvest2, 2026-11-06, 05:11–06:35 UTC.

Main finding: the Jul–Sep core-node papers (N1–N9) are already consumed by our rounds (e.g. R126, R129, R332, R345, R386); what remains open there is asymptotic or conditional. The live harvest is N11: a cluster of October Zenodo notes that state consequences of the 7/8 half-plane as "conditional on QRH" with constants never computed. QRH is Lean-checked (TREE N7), so these are unconditional theorems without numbers. Our explprog and chebexp machinery supplies the numbers.

## Ranked candidates

| # | id (source line) | node | statement in the source | what is missing | inputs (explicit where?) | h | gain | risk |
|---|---|---|---|---|---|---|---|---|
| 1 | arXiv 2609.25393 Thm 1.1 l.109–122, optimality l.129, Rem. l.251–253 (Klurman–Shparlinski–Teräväinen, 21 Sep); arXiv 2610.04678 Thm l.71–79, l.87–95 (Hughes, 3 Oct) | N11 | ≤ (log Q)^{1+ε} primes (KST) or (1/log 2 + o(1)) log Q primes (Hughes), all ≤ (log Q)^{C_A}, whose subset products give (ℤ/qℤ)^× for all q ≤ Q outside an exceptional set; KST l.129: the exceptional set "is unavoidable with present technology" without Vinogradov's least non-residue conjecture; C_*, C_A not computed | an all-moduli statement and any numbers | Lean-checked Re s > B_2 + 2.2·10⁻²⁵⁴ for every Dirichlet L (Bebel; TREE N7); taking σ < 1 − B in Hughes's lemma (ZF) (l.262, any σ ≤ 1/7 allowed), no character is bad (l.314), so his 𝓔_an (l.323–324) is empty and no repair is needed; explprog's per-character smoothed bound and least non-residue ≤ (log q)^{7.8978} (Riemann/Deep/ExplProg) | 4 | for every Q ≥ 3 one explicit set of ≤ (1/log 2 + η) log Q + c primes ≤ (c' log Q)^{7.90} that covers every q ≤ Q (subset products of its primes not dividing q); first all-moduli form, optimal count | medium: Hughes or KST v2 may add a QRH remark; Zenodo 23608018 (generation with arbitrary exponents only) and 23581412 (22 Oct, Regev relation lattices, almost all inputs) are near but different |
| 2 | Zenodo 23341593 (13 Oct), 23394711 (15 Oct), 23423890 (16 Oct); Zenodo 23623384 (24 Oct, Erdős 985) | N11 | every odd prime p has a prime primitive root g*(p) ≤ C(log p)^{24} (Martin 1997 Prop. 3 + 7/8 half-plane), C absolute, never computed; 23623384: g*(p) < p for p ≥ 1.4·10³³ plus a finite check | any numeric polylogarithmic (or p^{o(1)}) bound for g(p) or g*(p) valid for all p; unconditional explicit bounds are of size √p (Grosswald-type) | Martin, Acta Arith. 80 (1997) Prop. 3 (shifted sieve); explprog per-character bound; Bebel's B; the plain Vinogradov sieve gives only p^{O(1/log log p)}, useless in any computable range | 6+ | first numeric g*(p) ≤ (c log p)^{3/(1−B)} = (c log p)^{23.69} for every prime p: Vinogradov's least primitive root conjecture with numbers | high: explicit shifted sieve may overrun 6 h; crowded field (OpenAI "Primitive roots for every admissible integer base", 4 Oct, not read; Johnston–Naik 2610.08569) |
| 3 | Zenodo 23608373 (24 Oct) | N11 | P(n) ≤ h(n) ≤ max{P(n), ⌈C(log n)^8⌉} for large n (h(n) least h with gcd(2ⁿ−1,…,hⁿ−1) = 1, P(n) largest prime with p−1 \| n), "QRH-conditional resolution of Erdős Problem 770(iii)"; C and the threshold not computed | numbers, and the unconditional statement | explprog: primes ≤ (log q)^{7.8978} generate (ℤ/qℤ)^× (from its least non-residue for every character); the note's counting argument; Bebel's B | 4 | Erdős 770(iii) unconditional with explicit C and n₀: h(n) = P(n) whenever P(n) > (C log n)^{7.90} | medium: the implication was "previously recorded by Plasma AI" (Oct); erdosproblems.com/770 (read 6 Nov) is unchanged since Sep 2025, with one proof claim and six comments |
| 4 | arXiv 2609.07918 Thm local-pair l.229–239 (B. Wang); Zenodo 23027645 (29 Sep) | N8 | 23027645 derives every-T proportions of simple-or-real zeros in (T, T+T^θ] for θ > 0.2554312 from Wang's short-interval pair-correlation formula (W), which it calls an unrefereed "transposition" of Baluyot–Goldston–Suriajaya–Turnage-Butterbaugh | a verdict on (W) | Wang's TeX l.229–790; BGST (refereed); 23027645's unconditional inequality S(Z) ≥ x + (\|Z\|−x)²/N_m | 5 | audit: if (W) holds, the first every-T explicit proportion in (T, T+T^θ], 0.2554 < θ ≤ 7/22, becomes unconditional; if not, a REFUTATIONS entry for 2609.07918 | low pre-emption; medium value |
| 5 | Zenodo 23607492 (24 Oct) | N11 | exponent E(D) of the class group of ℚ(√D) satisfies E(D) ≥ (1/8 + o(1)) log\|D\|/log log\|D\| (Iwasawa 1969: does E(D) → ∞?), via Frei–Sofos–Widmer's small split prime; o(1) and constants unspecified | an explicit bound for every D | chebexp bound B for quadratic L (Riemann/Deep/ChebExp/REPORT.md: least split prime ≤ (0.62 log 3\|D\|)^{7.898}; Lean-checked here, since ζ_L = ζ·L(s,χ_D) and Bebel covers both at Re s > B_2 + 2.2·10⁻²⁵⁴; recheck 0.62 at that abscissa); 𝔭^E = (α) ⇒ p^E ≥ \|D\|/4 | 1 | first numeric E(D) ≥ log(\|D\|/4)/(7.898 log(0.62 log 3\|D\|)) for every fundamental D < −4 | low; routine (COROLLARIES) |
| 6 | Zenodo 23607811, 23608018 (24 Oct) | N11 | least prime quadratic residue r(p) ≪ (log p)^8; a prime ℓ ∤ Q with prescribed Legendre symbols modulo p_1, …, p_k, ℓ ≪ (2^k log Q)^8; constants absolute, not computed | numbers | chebexp (split prime in ℚ(√p*)); explprog smoothed sums with weights Π(1 + ε_iχ_i(n)) | 2 | first numeric r(p) ≤ (0.62 log 3p)^{7.898} and pattern bound | low; routine |
| 7 | Zenodo 23621868 (24 Oct) | N11 | every reduced class modulo a sufficiently large q contains a product of three primes ≤ q (QRH, mean-square argument); threshold not computed | explicit q₀, then all q by computation | explprog per-character bound; the note's orthogonality + mean square | 5 | every reduced class mod q contains p₁p₂p₃ with p_i ≤ q, for every q ≥ q₀ explicit (Erdős–Odlyzko–Sárközy-type) | medium: q₀ may be too large to finish by computation |
| 8 | Zenodo 21326653 (12 Jul) | N12 | distinct primitive characters of conductor ≤ Q differ in their zero ordinates below a polylogarithmic height, under GRH; unconditionally only polynomial height | unconditional polylog height | Bebel's half-plane (off-line zeros cost x^B in the explicit formula); the note's band-limited test function | 3 | unconditional distinction below (c log Q)^κ, κ ≈ 2/(1−B) | low; routine |
| 9 | arXiv 2609.24051 l.71–73, l.1133 (prime plus a square) | N10 | E(X) = O(X^{4/5}) for integers that are not a prime plus a square, constant ineffective (Siegel–Walfisz); "effectively only O(X^{9/10})" | effective X^{4/5} | no exceptional character under the half-plane; explprog (explicit Siegel–Walfisz with power saving) | 2 | effective E(X) ≪ X^{4/5}; no numeric constant without explicit circle method | low; routine |
| 10 | arXiv 2607.04283 abstract, l.86–88 (Acta Arith. 2008, posted 5 Jul) | N10 | π(x+N) − π(x) ≤ 2N/(log N + 3.53) "when N is large enough" | numeric N₀ | the paper's local-model sieve; explicit large sieve (Montgomery–Vaughan) | 6+ | first explicit N₀ | medium–high; off the RH path |
| 11 | arXiv 2609.00537 l.368 | N9 | explicit approximate functional equations need x, y ∈ ℤ + 1/2; "it should be possible to remove, a posteriori, the condition … retaining constants very close" | the removal, with constants | the paper's Table of ε₀ maxima and its §§ on the B-process | 4 | explicit AFE for all x, y | low gain: R386 found its B-process does not move our sub-Weyl certificate |
| 12 | arXiv 2609.06167 Thm l.113–122 | N11 | κ(q) ≥ 3/8 − ε (proportion of primitive χ mod q with L(1/2,χ) ≠ 0) for every sufficiently large q ≢ 2 (mod 4) | explicit q₀(ε) | explicit mollified first and second moments for all moduli (none available) | 6+ | first explicit non-vanishing proportion for all large q | high: not feasible in 6 h |

## Mini-briefs (top 5)

**1. Boolean covering for every modulus (N11).** Goal: for every Q ≥ 3, numeric c, c' and a set 𝒫 of ≤ (1/log 2 + η) log Q + c primes, all ≤ (c' log Q)^{1/(1−B)}, such that for every 2 ≤ q ≤ Q the subset products of {p ∈ 𝒫 : p ∤ q} are all of (ℤ/qℤ)^× (and Hughes's form: 𝒫 covers every q ≤ Q coprime to Π𝒫).
Inputs: Hughes 2610.04678 (Erdős–Rényi doubling Lemma cover, l.104–113; structured theorem l.87–95; 𝓔_an defined by zeros near 1); KST 2609.25393 (GRH sketch l.101–105, obstruction l.129); explprog DERIVATION.md and cert.py (|Σ χ(n)Λ(n)k(n/y)| ≤ c y^B log q per primitive χ); Bebel's B.
Method: per character, the weighted mean of χ(p) over primes ≤ y is ≤ δ once y ≥ (c log q/δ)^{1/(1−B)}; then Hughes's doubling with an empty exceptional set, or a simpler Hoeffding + union bound over the ≤ Q² characters (≈ 20 log Q primes) as a fallback; imprimitive χ and p | q handled explicitly.
Certificate: Arb for y(Q, δ) and the probability budget; brute-force check of the covering for all q ≤ 2000 with a printed prime set.
Verify: reviewer on the doubling lemma's constants and on moduli with many small prime factors (q = primorials, 2^k); novelty: KST and Hughes v2, Zenodo "cover the unit group"/"subset products", Oct–Nov listings 'cover|unit group|Regev|subset'.

**2. Least prime primitive root with numbers (N11).** Goal: numeric c with g*(p) ≤ (c log p)^{3/(1−B)} for every prime p (3/(1−B) = 23.69).
Inputs: G. Martin, "The least prime primitive root and the shifted sieve", Acta Arith. 80 (1997), Prop. 3 and its sieve lemmas; explprog per-character bound; Zenodo 23341593/23394711/23423890 (the QRH deduction) and 23623384 (explicit g*(p) < p under QRH, a structural cross-check).
Method: re-derive Prop. 3 with every constant tracked (shifted sieve over the prime factors of p − 1, character sums over primes ≤ x via the explicit formula with Re ρ ≤ B); small p by direct computation.
Certificate: Arb chain; sanity: exact g*(p) for p ≤ 10⁷ against the bound, and p with p − 1 highly composite.
Verify: reviewer on each sieve inequality; stop rule: if the explicit shifted sieve is not closed by 5 h, write up the partial constants and the blocking lemma.
Novelty: OpenAI "Primitive roots for every admissible integer base" (4 Oct), Johnston–Naik 2610.08569, Fan–Pollack, Zenodo 'primitive root', listings 'primitive root'.

**3. Erdős 770(iii) with numbers (N11).** Goal: explicit C, n₀ with P(n) ≤ h(n) ≤ max{P(n), ⌈(C log n)^{7.90}⌉} for n ≥ n₀; hence h(n) = P(n) whenever P(n) > (C log n)^{7.90}.
Inputs: Zenodo 23608373 (its subgroup estimate for 𝔽_p^× and the counting of products of distinct small primes); explprog least non-residue ⇒ generation by primes ≤ (log q)^{7.8978}; Bebel's B.
Method: substitute the explicit generation lemma, carry constants through the counting step; compute h(n) exactly for n ≤ 10⁴ (gcd of aⁿ − 1).
Certificate: Arb for the inequality chain; exact table h(n), P(n).
Verify: reviewer on the large-p branch of the counting argument; novelty: erdosproblems.com/770 (read only), Plasma AI's October record, Zenodo 'Erdős Problem 770'.

**4. Audit of Wang's short-interval pair correlation (N8).** Goal: verdict on 2609.07918 Thm local-pair (l.229–239): for H = T^θ, supp g ⊂ [−λ, λ], λ < θ, the zero-pair sum equals (HL/2π)(g(0) + ∫\|α\|g) + O(H + T^λ L²).
Inputs: Wang TeX l.229–790 (lemmas norm, localization, coefficients); BGST's published unconditional Montgomery theorem; Zenodo 23027645 (unconditional multiset inequality and the corollaries it draws from (W)).
Method: line-by-line check of the transfer to (T, T+H] (mean value of Dirichlet polynomials over a short interval, uniformity of Gonek's lemma, the weight w); numeric stress test of the formula with LMFDB zeros near T = 10⁶–10⁸ for λ close to θ.
Certificate: if sound, Arb values of c*(C_θ) at θ = 0.26, 0.3, 7/22; if not, the failing line and a numeric witness.
Verify: reviewer verdict; known: Wang v2, citing papers, 23027645 versions.

**5. Iwasawa's class-group exponent with numbers (N11).** Goal: E(D) ≥ log(\|D\|/4)/(7.898·log(0.62 log 3\|D\|)) for every fundamental discriminant D < −4.
Inputs: Riemann/Deep/ChebExp (bound B for quadratic L/ℚ, C = {1}); the 3-line argument 𝔭^E = (α), α ∉ ℤ ⇒ p^E ≥ \|D\|/4 (α ∈ ℤ would force 𝔭 = 𝔭̄); Zenodo 23607492.
Method: one page; same lemma gives row 6 (r(p) via ℚ(√p*)).
Certificate: PARI table of E(D) for \|D\| ≤ 10⁶ against the bound.
Verify: reviewer checks that bound B for quadratic fields uses only ζ and L(s,χ_D) (Bebel, Lean-checked), not Liu's Dedekind B_2; novelty: FSW 2610.18491, Zenodo 'class group exponent', listings 'exponent'.

## Dropped hits
- Already ours: 2607.16795 (gap k < 10¹⁸r³ in NOGO, R332), 2608.08682 (hollandk), 2608.19160 (R332; it solves Open Problem 4.14 itself), 2608.24827 and 2609.20367 (Weil rounds), 2609.37164 (its Λ certificate at τ₀ = 0.161 is far above our 0.10756), 2609.23222, 2608.01498, 2609.26787, 2609.33794, 2609.27715, 2607.07566, 2609.02882, 2608.13637, 2609.33043; 2609.24167 has no source (404) and is superseded by 2610.08965.
- Superseded by the Lean-checked half-plane: 2609.22624 (its 13 % with real part in [1/16, 15/16] is 100 %, since every zero has 1 − B ≤ Re ρ ≤ B); 2608.31060 (its conditional S(x,t) ≪ x/(log x)^{1/100} follows from μ(σ) = 0 for σ > B); 2608.07399, 2609.30207 (consequences of a Siegel zero, now vacuous); 2609.35950 (Landau-type real-zero bound; 2610.21511 is stronger); 2607.22930, 2607.15574, 2609.07371, 2608.23500, 2608.07198, 2607.27282 (ineffective only through Siegel–Walfisz: effective via explprog, numbers need full circle-method work).
- Conditional on GRH or other hypotheses: 2607.06433, 2607.14515, 2609.14471, 2609.16048, 2609.23390, 2609.02392, 2608.05961, 2608.15773, 2609.18569, 2609.21195, 2608.21573, 2607.04316, 2608.27130, Zenodo 21585399, 21581938.
- Asymptotic only (no constant to compute): 2609.15329, 2609.07918 Thm main, 2608.16034, 2609.27808, 2609.11619, 2609.01101, 2608.09906, 2609.18938, 2609.18939, 2609.25985, 2607.03654.
- Input missing or off the tree: 2609.11070 l.1477 (C(A) needs an explicit Carleson–Hunt constant); 2609.28526 (already effective); 2609.14295 (partition sequences); 2609.17806, 2609.16893, 2607.10053, 2607.04366; Zenodo 23527836 (Eisenstein series, already explicit); Zenodo 23870366 (CCM: its author notes the missing step is equivalent to RH).

## Flags for our files (novelty)
- discoveries/littlewood-half-plane.md says "before explicit only 1/log q". Zenodo 23400116 (15 Oct) and 23459046 (18 Oct, Lean 4) prove L(1,χ) ≥ 0.0621/log log q (earlier 0.0617) for real primitive χ, q ≥ 10¹⁰⁰, from Re s > 7/8. Ours (0.11, every primitive χ, every q) is still stronger, but the "before" line and Littlewood78/REPORT.md l.147 ("Zenodo: nothing relevant") need correcting.
- Riemann/Deep/ExplProg/NOVELTY.md l.5 ("nonresidue quasi-Riemann: nothing relevant") misses Zenodo 23442654 (17 Oct), 23442990 (18 Oct) and 23702717 (27 Oct): n(p) ≪ (log p)^8 under QRH, constants not computed. explprog's numeric claim stands; cite them.
- The imqclass page omits Zenodo 23547020 (21 Oct): all fields with h ≤ 3000, conditional on Kimura's real-zero bound (Lean-checked by its author, unrefereed), between 2610.21511 (2000) and ours (5000).
- Harvest1 scanned 2610.04678 (Hughes) without noting that the half-plane empties its exceptional set (candidate 1 here).

## Papers scanned
2607.00282 hit · 2607.00332 none · 2607.01184 hit · 2607.01956 hit · 2607.02828 none · 2607.03654 hit
2607.04283 hit · 2607.04316 hit · 2607.04366 hit · 2607.04476 hit · 2607.04632 hit · 2607.06022 none
2607.06433 hit · 2607.07566 hit · 2607.07911 hit · 2607.08318 hit · 2607.08342 none · 2607.09056 hit
2607.09110 hit · 2607.09721 hit · 2607.10053 hit · 2607.11068 hit · 2607.12084 hit · 2607.13490 none
2607.13518 none · 2607.14515 hit · 2607.14612 hit · 2607.14981 hit · 2607.15137 hit · 2607.15311 hit
2607.15574 hit · 2607.16032 hit · 2607.16695 hit · 2607.16795 hit · 2607.17249 hit · 2607.17418 hit
2607.17731 none · 2607.20041 hit · 2607.21259 hit · 2607.21532 hit · 2607.21883 none · 2607.22930 hit
2607.23150 hit · 2607.24311 none · 2607.25002 hit · 2607.25088 hit · 2607.25403 hit · 2607.27131 hit
2607.27282 hit · 2607.28091 hit · 2607.28931 hit · 2608.01337 hit · 2608.01498 hit · 2608.02381 hit
2608.04249 hit · 2608.04391 hit · 2608.04923 hit · 2608.04948 hit · 2608.04982 hit · 2608.05934 hit
2608.05961 hit · 2608.07198 hit · 2608.07399 hit · 2608.07804 none · 2608.08149 hit · 2608.08682 hit
2608.09488 hit · 2608.09906 hit · 2608.10378 hit · 2608.12257 hit · 2608.12315 hit · 2608.12400 hit
2608.13266 none · 2608.13299 hit · 2608.13509 hit · 2608.13637 none · 2608.13687 hit · 2608.14777 hit
2608.14852 hit · 2608.15063 hit · 2608.15773 hit · 2608.16034 hit · 2608.16679 hit · 2608.17045 hit
2608.17199 hit · 2608.18197 hit · 2608.19160 hit · 2608.19262 hit · 2608.20431 none · 2608.21573 none
2608.22717 hit · 2608.23500 hit · 2608.24827 hit · 2608.24878 hit · 2608.24930 hit · 2608.25721 none
2608.25812 hit · 2608.27130 hit · 2608.27323 hit · 2608.27356 hit · 2608.29558 none · 2608.29938 hit
2608.30061 hit · 2608.30920 none · 2608.31060 hit · 2608.31126 hit · 2609.00537 hit · 2609.01101 hit
2609.01358 hit · 2609.02327 hit · 2609.02392 hit · 2609.02478 hit · 2609.02713 hit · 2609.02882 hit
2609.03961 hit · 2609.04092 hit · 2609.04155 hit · 2609.04993 hit · 2609.05950 hit · 2609.06167 hit
2609.07080 hit · 2609.07371 hit · 2609.07918 hit · 2609.09942 hit · 2609.10328 none · 2609.10446 hit
2609.11070 hit · 2609.11119 hit · 2609.11619 hit · 2609.11855 hit · 2609.12764 none · 2609.13095 hit
2609.13894 hit · 2609.14295 hit · 2609.14714 hit · 2609.14909 hit · 2609.14915 none · 2609.15166 hit
2609.15329 hit · 2609.15862 hit · 2609.15937 hit · 2609.16181 hit · 2609.16360 hit · 2609.16678 hit
2609.16807 hit · 2609.16893 hit · 2609.17478 hit · 2609.17806 hit · 2609.17875 hit · 2609.18220 hit
2609.18414 hit · 2609.18569 hit · 2609.18938 hit · 2609.18939 hit · 2609.20367 hit · 2609.20505 hit
2609.20522 hit · 2609.22624 hit · 2609.22777 hit · 2609.23222 hit · 2609.23390 hit · 2609.24051 hit
2609.24167 nosrc · 2609.24347 hit · 2609.25171 hit · 2609.25393 hit · 2609.25879 none · 2609.25885 hit
2609.25985 hit · 2609.26786 hit · 2609.26787 none · 2609.27434 hit · 2609.27715 hit · 2609.27808 hit
2609.28200 hit · 2609.28253 hit · 2609.28526 hit · 2609.29898 none · 2609.30207 hit · 2609.32053 hit
2609.32673 hit · 2609.33043 none · 2609.33278 none · 2609.33692 hit · 2609.33794 hit · 2609.35496 hit
2609.35950 hit · 2609.36200 none · 2609.37164 hit · 2609.37466 hit · 2609.39706 hit
zen:23027645 hit · zen:23341593 hit · zen:23394711 hit · zen:23423890 hit · zen:23623384 hit · zen:23608373 hit · zen:23607811 hit · zen:23608018 hit · zen:23607492 hit · zen:23621868 hit · zen:23643108 hit · zen:23645699 hit · zen:21326653 hit · zen:23400116 hit · zen:23459046 hit · zen:23442654 hit · zen:23442990 hit · zen:23702717 hit · zen:23405184 hit · zen:23547020 hit · zen:23461828 hit · zen:23527836 hit · zen:23870366 hit · zen:23855214 hit · zen:23007659 hit · zen:22898988 hit · zen:21581938 hit · zen:21585399 hit · zen:23396432 hit · zen:23436577 hit · zen:23881837 hit · zen:23608576 hit · zen:23462601 hit · zen:23605077 hit · zen:23793085 hit
("hit" = at least one phrase match, strong or weak; "nosrc" = no TeX source (HTTP 404); "zen:" = Zenodo record, unrefereed, description read.)

RH relevance: N11, N8, N12, N10, N9
