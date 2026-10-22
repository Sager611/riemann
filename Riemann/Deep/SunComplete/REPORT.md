DISCOVERY: Theorem A (Sun–Hou Conjecture 3.5(i), Nanjing Univ. J. Math. Biquarterly 36 (2019)). For every prime power q > 7 there is a numbering a_1, …, a_q of F_q in which all cyclic adjacent sums a_i + a_{i+1} are primitive. Theorem B (Sun Conjecture 3.8, corrected; this classification was posted first by an anonymous Zenodo preprint, 10.5281/zenodo.23307915, so only Theorem A is new). For every q > 7 and every a₀ ∈ F_q with (q, a₀) ≠ (11, 2) there is a circular permutation of F_q^* in which every a₀ + a_i a_{i+1} is primitive; (11, 2) has none.

# REPORT — suncomplete (2026-10-22)
- **Status.** Both conjectures are now proved for all q > 7. Before this, both were known only for q > 18,888,871 (She, arXiv:2609.19114; Luo–She, arXiv:2610.06289), plus A for primes p > 61 and B for a₀ = 0.
- **Graph bound.** Chvátal–Erdős κ ≥ α with two sharper bounds: α via e(I, I) ≤ |I| (Hoffman-type factor n/(Ñ + E), no splitting of I), and κ from the exact cut inequality Ñ²ab ≤ E²(a+c)(b+c). Together they reduce She's 3·W√q to about qE/(q+E).
- **Sieve.** Cohen's prime sieve gives E = (W(k)Δ − 1)√q + 1, via Gauss sums for A and Jacobi sums for B.
- **Analytic coverage.** Uniformly, every q ≥ 30,140 (Lemma 6, by cases on ω(q−1)). Checked per q, every q < 30,140 except 79 values for A and 83 for B, the largest being 6,091. So the analytic threshold drops from 18,888,871 to 6,091.
- **Exceptions for A.** She's three-consecutive-primitive-roots construction for 62 primes; explicit cycles for q = 9, 13, 16, 25, 27, 29, 49, 61, 121, 169, 361, 841, 1681.
- **Exceptions for B.** All 59,012 pairs (q, a₀ ≠ 0) are certified. Two new explicit constructions do almost all the work: alternating products u, ur (u a non-square, r a generator of the squares) for 58,962 pairs, and geometric products u, us, us² (q even) for 22. The other 27 pairs have explicit cycles. Exhaustive search confirms (11, 2) has no Hamilton cycle.
- **a₀ = 0.** Covered by Luo–She's explicit construction, re-proved in NOTE §6.
- **Cross-check.** Independent of the analytic part, compiled C code checked every q < 30,140 directly: 46.7M pairs (q, a₀), plus all 74 non-prime fields for A. Results are stored as checksums.
- **Sharpness.** A fails for q = 3, 5, 7 and holds for q = 4. B fails for q = 7 and every a₀ ≠ 0.
- **Rigour.** All inequalities use exact `Fraction` arithmetic, with √q replaced by rational upper bounds. Field arithmetic is honest: F_p[x]/(f), with x verified to have order q−1. A brute-force DFT on 16 fields confirmed the character-sum bound (ratio at most 0.84). For q ≤ 81, exact α and κ (networkx) agree with the proved bounds.
- **Novelty.** arXiv still has only v1 of both papers. An arXiv search on 2026-10-22 found no completion.
- **Files** (`Riemann/Deep/SunComplete/`, under 1 MB):
  - `NOTE.md` — full proofs.
  - `verify.py` — about 5–10 s; checks Lemma 6, re-derives E_A and E_B, checks all certificates and rebuilds 541 (B)-cycles; `--cross` reruns the C checks.
  - `certA.json`, `certB.json.gz` — certificates.
  - `gen_certs.py`, `gen_table.py` — generators.
  - `crossA.c`, `crossB.c`, `crossA.out`, `crossB.out` — the C cross-check.
- **Numbering.** Following the printed paper, the conjectures are 3.5(i) and 3.8; the arXiv v9 TeX labels them Conj3.7 and Conj3.10.
- **Open.** Still open: a character-sum proof that the alternating construction (Lemma 9) alone works for all odd q > 43. It holds for every q < 30,140 tested.

## Review (2026-10-22)
> INCOMPLETE — **FLAWED as written:** `NOTE.md:267–275`, Lemma 10 omits **u ≠ 0**.
> Significance: A’s proof checks out; B needs only this local correction, already enforced for every actual certificate.
> `REPORT.md:1`: “DISCOVERY” overstates B’s novelty; the earlier Zenodo preprint already claims that classification.
> Counterexample to Lemma 10: q=8, u=0, a₀=s primitive satisfies all stated hypotheses, but t=0; the proposed cycle cannot enumerate F₈*.
> The independence/connectivity bounds, inequality directions, and Chvátal–Erdős application are sound for these simple graphs with n≥3.
> The Gauss/Jacobi sieve estimates and Lemma 6’s finite table plus infinite-tail argument check out.
> `verify.py` ran successfully: all checks passed in 6.6 seconds, including exhaustive exceptional-pair coverage and extension-field validation.
> Sampling reconstructed B-cycles is not a gap: every construction’s algebraic conditions are checked, and the nonzero-u construction is proved.
> TODO @author: add **u ∈ F_q*** to Lemma 10; qualify B’s novelty in `REPORT.md:1`.
- Fixed: Lemma 10 now assumes u ∈ F_q^* (every certificate already had u ≠ 0); line 1 now qualifies the novelty of Theorem B.
