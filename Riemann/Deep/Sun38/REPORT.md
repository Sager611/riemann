COUNTEREXAMPLE: Sun's conjecture fails exactly for (q, a0) ∈ {(11, 2)} among q ≤ 100000. Novelty: KNOWN — anonymous Zenodo preprint 10.5281/zenodo.23307915 (2026-10-12) has the same proof and claims (11, 2) is the only exception for every q > 7. We found it independently on 2026-10-22 (Riemann/Deep/SunPrim), ten days later.

# Sun38: circular permutations with primitive values, plus three more Sun conjectures

Orca agent sun38, 2026-11-05. arXiv 1309.1679 v9 (2 Mar 2020) is the journal paper (Nanjing Univ. J. Math. Biquarterly 36 (2019) 134–155). Its TeX labels Conj3.3, Conj3.4, Conj3.5 and Conj3.10 are printed as Conjectures 3.2, 3.3, 3.4 and 3.8.

## 1. Conj3.10 (printed 3.8): evidence
- **Statement** (v6 of 24 Sep 2013 to v9, and the journal copy on Sun's homepage, 196a.pdf of 8 May 2020): for q > 7 and any a0 ∈ F_q, some circular permutation (a_1, …, a_{q−1}) of F_q^* has every a0 + a_i a_{i+1} primitive. In graph terms, G_{q,a0} on F_q^* (x ~ y iff x ≠ y and a0 + xy is primitive) is Hamiltonian. v5 (19 Sep 2013) had only a0 = ±1, and that case never fails here.
- **The exception.** The primitive roots mod 11 are 2, 6, 7, 8, so in G_{11,2} we have x ~ y iff xy ∈ {4, 5, 6}. Then N(2) = N(9) = {3, 8}. A Hamilton cycle must use both edges at a degree-2 vertex, so it would contain the closed 4-cycle 2–3–9–8–2 on only 4 of the 10 vertices. Exhaustive search agrees (5 nodes).
- **Method.** Write x = g^i with g primitive and n = q − 1. Then x ~ y iff i + j ∈ T, where T = {s : a0 + g^s primitive}, a sum graph on Z/n.
  *Lemma (two products, n even).* Suppose t ∈ T is odd, t + d ∈ T and gcd(d, n) = 2. Then x_{2m} = g^{md}, x_{2m+1} = g^{t−md} (0 ≤ m < n/2) is a Hamilton cycle. The even positions run through the squares and the odd positions through the non-squares. Adjacent products alternate g^t and g^{t+d}; the closing product is g^{t+d} because (n/2)d ≡ 0. Luo–She's cycle for a0 = 0 is the case t = −1, d = 2.
  Variants (proofs in the gen.py header): 'u' (n even, t odd, d a unit, t + d ∈ T and t + d + n/2 ∈ T); 'a' (n odd, t, t + e, t + 2e ∈ T with e a unit); 'z' (a0 = 0 and n odd: Luo–She's cycle). If none applies, a cycle is found by search: a rotation heuristic, then exhaustive search pruned by forced edges and connectivity.
- **Results.** Range 7 < q ≤ 100,000: 9,695 prime powers, every a0, 457,028,131 pairs (q, a0).

| range | how | pairs | outcome |
|---|---|---|---|
| q ≤ 2000 (328 fields) | gen.py; every certificate stored (data/certs.txt.gz) | 290,830 | 288,762 two-product, 2,036 'a', 6 'u', 4 'z', 21 explicit cycles, 1 exception |
| primes 2000 < q ≤ 100,000 (9,289) | ext.c: two-product certificate for every a0, each re-checked by modular exponentiation | 454,119,487 | all Hamiltonian |
| prime powers 2000 < q ≤ 100,000 (78) | gen.py; re-derived separately by verify.py --recompute-pp | 2,617,814 | all Hamiltonian |

- **Where search was needed.** Explicit cycles were needed only at q = 9 (4 values of a0), 11 (3), 13 (8), 25 (4), 31 (a0 = 3) and 43 (a0 = 5). Every case left open by the earlier scan is Hamiltonian. (31, 3) has an explicit cycle. The other seven have two-product certificates (t, d): (31, 11) → (11, 14); (37, 13) → (1, 22); (37, 25) → (5, 2); (37, 34) → (3, 14); (43, 10) → (19, 8); (43, 12) → (15, 8); (43, 38) → (23, 8).
- **Side observation** (checked, not proved). For odd q with 43 < q ≤ 100,000, every a0 admits a permutation for Sun's conjecture whose adjacent products take only two values. Two-product certificates fail only at q ∈ {9, 11, 13, 25, 31, 43}.
- **Checks.**
  - `python3 verify.py` (stdlib only, 90 s) rebuilds every field from the stored modulus and generator. The generator must have order q − 1, which also proves the modulus irreducible.
  - It then checks every stored certificate and writes out and edge-checks the full permutation for all 82,631 pairs with q ≤ 1000. It also re-proves (11, 2) by exhaustive search, checks the extension summaries, and recomputes one random extension prime.
  - `--recompute-pp` re-derives all 78 extension prime powers (191 s; done, OK). `--recompute LO HI` re-derives any range.
  - Run times: ext.c 730 s single-threaded in 1.8 MB; gen.py for q ≤ 2000 takes 6 s.

## 2. Novelty of (11, 2): KNOWN
- **Zenodo 10.5281/zenodo.23307915**, "Sun's circular primitive-permutation conjecture: one exception and a complete classification" (Anonymous, 2026-10-12, v1.1.0-candidate, unrefereed; code at github.com/ipitchford/sun-circular-primitive-permutation). Its Prop. 3 is our 4-cycle proof. Its Thm 2 says G_{q,a0} is Hamiltonian iff (q, a0) ≠ (11, 2), for all q > 7: a sieve leaves 510 fields up to 92,821, which are closed by heuristic search, and the cycles are not kept. Our run covers every field up to 100,000 with a different construction and stored or recheckable certificates, so it confirms their finite part on its own.
- **Luo–She, arXiv 2610.06289** (only v1, 5 Oct 2026) proves the conjecture for q > 18,888,871 and, for a0 = 0, for all q > 4. It mentions no small exceptions. Its "she" reference is Y.-F. She, arXiv 2609.19114 (16 Sep 2026), which treats prime circles and the Sun–Hou sum version (TeX Conj3.7(i), printed 3.5(i)), not 3.8.
- **Other sources.** Sun made no correction in v6–v9 or in the journal copy, and his later paper 1811.10503 does not restate the conjecture. OEIS A229232 states only a0 = ±1; its revision 34 (27 Oct 2026) added terms only. Semantic Scholar lists 11 papers citing 1309.1679; only Luo–She discusses 3.8, and the newest (2610.18431, Ai–Hu–Yan) is on reciprocal sums.

## 3. Extra scope: three more conjectures of 1309.1679 (own code: verify.py part B, tables.py)
**(a) Conj3.3 (printed 3.2, 2013-09-02).** If n ∤ |G| (or n is even and the Sylow 2-subgroup is cyclic), some ordering with b_1 = a_1 has the n − 1 differences b_i − b_{i+1} pairwise distinct.
- Smallest counterexample: |G| = 6, n = 4; nothing fails for |G| ≤ 5. Take G = Z/6, (a_1, …, a_4) = (0, 1, 3, 5) and b_1 = 0. Then b_2, b_3, b_4 run through {1, 3, 5}. The differences b_2 − b_3 and b_3 − b_4 lie in {2, 4} and must differ, so they sum to 0. But their sum is b_2 − b_4 ≠ 0.
- General form: a_1 together with a coset of a subgroup H with |H| ≥ 3 and Σ H = 0. The n − 2 inner differences would be all of H ∖ {0}, which sums to 0, yet their sum is b_2 − b_n ≠ 0. All 6 failing pairs (A, a_1) in Z/6 are of this kind (tables.txt). It also fails in Z/2×Z/4 and (Z/2)³ (n = 5).
- The text is the same from v1 (6 Sep 2013, then Conj. 2.3) through v9 and the journal copy.

**(b) Conj3.4(ii) (printed 3.3(ii) in v9, 2013-09-03).** If |G| is odd, n > 2 and n ∤ |G|, some cyclic order has the n differences a_i − a_{i+1} pairwise distinct.
- Smallest counterexample: Z/5, n = 3; all ten 3-subsets fail. For A = {0, 1, 2} the two cyclic orders give differences (−1, −1, 2) and (−2, 1, 1).
- General form: every (|G| − 2)-subset of an odd-order G fails. The n differences are distinct, nonzero and sum to 0, so they would be all of G ∖ {0, m}, but that set sums to −m ≠ 0.
- Older wordings fail as well:

| version | wording of (ii) | smallest counterexample |
|---|---|---|
| v1–v3 | "n < |G|" (with 2 ∤ n or n ∤ |G|) | Z/4, {0, 1, 2} |
| v4–v8 | "3 < n < |G|" | Z/6, {0, 2, 3, 4}; also every 5-subset of Z/7 |
| v9 TeX | "n ∤ |G|" | Z/5, {0, 1, 2} |
| journal copy | "n = |G| − 1" | none: this is the known theorem that odd-order abelian groups are R-sequenceable |

- So the printed version is correct. Sun changed it in 2020 without saying why.

**(c) Conj3.5(i) (printed 3.4(i), 2013-09-20).** If |A| = n > 3 and G is finite with 3 ∤ |G|, some cyclic order has the n sums a_i + 2a_{i+1} pairwise distinct.
- Smallest counterexample: |G| = 7, n = 5; nothing fails for |G| ≤ 5. Every 5-subset of Z/7 fails. The map x ↦ αx + β sends each sum s to αs + 3β, and the affine group of Z/7 is 2-transitive, so all 21 five-subsets behave like {0, 1, 2, 3, 4}. Each of its 24 cyclic orders repeats a sum (tables.txt).
- It also fails in Z/8 (n = 4, 5, 6: 8, 8 and 4 subsets, e.g. {0, 1, 2, 5}) and in Z/2×Z/4 (n = 5).
- v6 (24 Sep 2013) had no condition on G; v7 (31 Oct 2013) added 3 ∤ |G|. The text is unchanged since then, journal copy included.

**Later sources.** None found. The 20 OEIS entries citing 1309.1679 carry no note; the only one changed since 22 Oct is A229232, with new terms. The citing papers listed above say nothing about these failures. The 251 GitHub repositories of the Zenodo author cover only 3.8 and Sun–Hou.

**Status in this project.** All three were already found on 22 Oct 2026 (Riemann/Deep/SunPrim, NOTE.md and NOVELTY.md). (a) and (c) are already in DISCOVERIES.md ("Seven published conjectures of Z.-W. Sun fail…"). (b) is excluded there because the printed version is corrected. Today's new items are only the extension of 3.8 to Q = 100,000 and the version table for (b).

## 4. Files (Riemann/Deep/Sun38/)
- verify.py: checker (stdlib only).
- gen.py: Python generator (certificates, search, exact search).
- ext.c: C sweep over primes.
- tables.py → tables.txt: case tables for (a) and (c).
- data/certs.txt.gz: every certificate and cycle, 7 < q ≤ 2000.
- data/exceptions.json: the exception with its proof.
- data/meta.json: stored and extension ranges.
- data/ext_primes.txt.gz: per-prime counts from ext.c.
- data/ext_pp_stats.json: per-field counts for the prime powers.
