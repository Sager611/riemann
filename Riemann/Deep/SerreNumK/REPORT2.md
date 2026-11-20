RESULT: Momose Type 2 isogeny primes are bounded without GRH over every number field k, and so Larson–Vaintrob's isogeny corollary holds without GRH. (1) No k-rational p-isogeny of Momose Type 2 has log p ≥ y*(k) ≤ max(3000n_M, 2 log Δ_M + n_M log 4), M the Galois closure of k. (2) Certified: y*(k) = 162.0–218.0 on Banwait–Derickx's table fields of degree 2–10 (p < 10^70.4 … 10^94.7, against their GRH bounds 5.65·10^10 … 1.24·10^14), and y* = 198.75 for ℚ(2^{1/3}). (3) If k contains no Hilbert class field of an imaginary quadratic field, every k-isogeny prime divides MMIB(k)·TypeOneBound(k) or is < e^{y*(k)}. So IsogPrimeDeg(k) is finite with an explicit bound (LV14 Cor. 'isogeny'). (4) LV14 Thm 1 alternative (2) is void for log ℓ ≥ y*(K). (5) Allen–Genao 2511.23381's GRH statement (l.77) becomes unconditional with the enlarged exceptional set S_F ∪ {p : log p < y*(F)} in place of S_F; case (1) of LV is excluded by its own conclusion, a CM curve over F with CM field inside F. Genao 2405.05507 is only partly settled: one of its steps rests on [Gen24], which was not checked.

# SerreNumK follow-up: Lemma 2 for every number field; Larson–Vaintrob without GRH (node N11)

Agent serrenumk, 20 Nov 2026, 00:18–00:32 UTC (one process at a time, ≤ 35 MB). Full proofs: `DERIVATION2.md`.

## Claim
Write y = log p, M = Galois closure of k, n_M = [M:ℚ], Δ_M = |disc M|.
- **Lemma 2(k) (every number field k; no hypothesis).** No E/k has a k-rational p-isogeny of Momose Type 2 (λ^{12} = χ_p^6)
  with y ≥ y*(k). Lemma 2′: y*(k) ≤ max(3000n_M, 2 log Δ_M + n_M log 4).
- **Theorem 3 (LV14 Cor. 'isogeny' without GRH).** Let k contain no Hilbert class field of an imaginary quadratic field.
  If some E/k has a k-rational p-isogeny, then p | MMIB(k)·TypeOneBound(k) or log p < y*(k). Hence IsogPrimeDeg(k) is
  finite, with an explicit bound, unconditionally.
- **Corollary 4 (LV14 Thm 1).** For ℓ ∉ S_K with log ℓ ≥ y*(K), alternative (2) ("GRH fails for K[√−ℓ]", ψ^{12} = χ_ℓ^6)
  cannot occur. This is Harvest4 item 5's goal.
- **Corollary 5 (Allen–Genao 2511.23381, l.77, without GRH).** Let F have no rational CM (no imaginary quadratic subfield,
  their l.75). Let E/F and p satisfy the hypotheses of l.77 (F(E[p])/F abelian, p ≥ 17, p ∤ Δ_F or p > (6[F:ℚ])! + 1).
  Then p ∈ S_F ∪ {p : log p < y*(F)}: an enlarged exceptional set, not necessarily LV's S_F.
  - LV's case (1), verbatim (1110.0255 l.192–196): "There exists a CM elliptic curve $E'$, which is defined over $K$ and
    whose CM-field is contained in $K$". That needs an imaginary quadratic field inside F, so it is impossible.
  - Case (2) is void by Corollary 4.

| k (LMFDB) | d | n_M | log Δ_M | **y*(k)** | Type 2 primes < | BD's GRH bound | Harvest4 route (float) |
|---|---|---|---|---|---|---|---|
| 2.2.5.1 | 2 | 2 | 1.609 | 162.00 | 10^70.36 | 5.65·10^10 | 10^139.5 |
| 3.3.49.1 | 3 | 3 | 3.892 | 172.00 | 10^74.70 | 4.09·10^11 | 10^151.0 |
| 4.0.125.1 | 4 | 4 | 4.828 | 183.50 | 10^79.70 | 1.46·10^12 | 10^160.8 |
| 5.5.14641.1 | 5 | 5 | 9.592 | 192.25 | 10^83.50 | 4.75·10^12 | 10^168.3 |
| 6.6.300125.1 | 6 | 6 | 12.612 | 198.75 | 10^86.32 | 1.12·10^13 | 10^174.4 |
| 7.7.594823321.1 | 7 | 7 | 20.204 | 204.50 | 10^88.82 | 2.65·10^13 | 10^179.6 |
| 8.0.64000000.2 | 8 | 8 | 17.974 | 209.50 | 10^90.99 | 4.16·10^13 | 10^183.8 |
| 9.9.16983563041.1 | 9 | 9 | 23.556 | 214.00 | 10^92.94 | 7.60·10^13 | 10^187.7 |
| 10.10.572981288913.1 | 10 | 10 | 27.074 | 218.00 | 10^94.68 | 1.24·10^14 | 10^191.2 |
| ℚ(2^{1/3}), non-Galois | 3 | 6 | 10.463 | 198.75 | 10^86.32 | — | 10^151.1 |
| x³ − x² − 2x − 20, non-Galois | 3 | 6 | 20.016 | 199.00 | 10^86.43 | — | 10^151.4 |

The first nine rows are BD's 'tab:output_summary_type_2': the smallest Galois fields of class number one without an
imaginary quadratic HCF. |Δ_M| for the cubics comes from PARI nfdisc(nfsplitting(·)) (logs/galclosure.txt). The last
column is floating point and only a comparison, see Method.

## Before
- **Banwait–Derickx** (ANT 2025; arXiv 2203.06009v3). They make Momose's classification explicit and algorithmic for every
  k: MMIB(k) and TypeOneBound(k) are unconditional, and the passage to non-Galois k is fixed. Their Type 2 bound is
  p ≤ (8d log(12p) + 16 log Δ_k + 10d + 6)⁴ under GRH, "the only result in our work which requires the Generalised
  Riemann Hypothesis" (TeX l.1457).
- **Larson–Vaintrob** (Compositio 2014; arXiv 1110.0255).
  - Thm 1 is unconditional, but keeps alternative (2) "GRH fails for K[√−ℓ]".
  - Cor. 'isogeny' (finiteness of IsogPrimeDeg(K) for K without an imaginary quadratic HCF) is stated "Under GRH".
  - BD restate it as their Thm 1.1: "Assume GRH ... finite iff ..." (TeX l.210–214).
- **Unconditional before:** semistable curves (BD; Type 2 cannot occur), ℚ (Mazur), specific fields. Over the real quadratic
  fields with D < 50, our REPORT.md Lemma 2 gave y* = 162.

## Method
- **Lemma 2(k).** This is REPORT.md's argument with M in place of the quadratic field.
  - chebexp over ℚ on F = M(√−p), trivial class: q₁. Here log Δ_F ≤ 2 log Δ_M + n_M log 4p, from 𝔡_{F/M} | 4p.
  - chebexp on F(√q₁), trivial class: q₂ ≠ q₁. Here log Δ ≤ 2 log Δ_F + 2n_M log 4q₁.
  - Both q_i split completely in M, so they have residue degree 1 in k and split in ℚ(√−p).
  - Both are < p/4 once y ≥ y*(k). They cannot both be primitive cube roots of 1 mod p.
  - So BD's unconditional Condition CC (Prop. 'cond:CC') fails.
- **Lemma 2′.** A closed-form tail, valid for every n_M ≥ 1, using G(R) ≤ R for R ≥ 171.
- **Theorem 3.** BD's case split:
  - not of Momose Type 1, 2 or 3: p | MMIB(k), including signature Type 2 that is not Momose Type 2 (Cor. 'type_2_not_momose');
  - Momose Type 1: p | TypeOneBound(k);
  - Momose Type 2: Lemma 2(k);
  - Momose Type 3: needs an imaginary quadratic HCF inside k.
- **Corollary 4.** LV case (2) forces reducibility over F_ℓ (LV l.205–206). Then λ^{12} = χ_ℓ^6 for the isogeny character
  (ψ ∈ {λ, χ_ℓλ^{−1}}), which is Momose Type 2.
- **The Harvest4 item 5 route, compared.** It puts chebexp into LV's own step (Cor. 'avo' on K(√±ℓ, √±3), then
  ℓ ≤ (1 + √N v)⁴).
  - That gives exponent 64 instead of 32: y ≈ 321–440 against our 162–218.
  - It needs the claim "a split prime is not over 3", which fails when 3 ramifies in K with even index, e.g.
    ℚ(√3)(√−3)/ℚ(√3) is unramified above 3.
  - CC has no condition at 3 and needs only q < p/4.

## Certificate (one command)
```
cd Riemann/Deep/SerreNumK && uv run --no-project --with python-flint python3 cert2.py > logs/cert2_output.txt
```
- Ends with `ALL CERTIFIED (7.1 s)`, one core, max RSS 35 MB (logs/cert2_time.txt).
- Per field it runs a chain from y = 100 to Y₀ = max(3000n_M, ⌈2 log Δ_M + n_M log 4⌉ + 1): width ¼ below 2000,
  width 4 above.
  - Each interval [y_a, y_b] checks log q_i + log 4 < y_a, with R_F and R₂ taken at y_b.
  - y* is the right end of the last failing interval.
  - Above Y₀ it checks Lemma 2′'s inequalities at Y₀ in Arb.
- The envelope constants are the same as in cert.py (Arb balls from decimal strings).

## Checks
- **Consistency.** For 2.2.5.1, y* = 162.00 matches REPORT.md's real-quadratic Lemma 2 (same R_F and R₂). y* grows slowly:
  log y* rises with log n_M, since the envelope bottleneck is G = 151.17 or 0.108R. Every y* lies far above BD's GRH
  bounds (log 5.65·10^10 = 24.8, log 1.24·10^14 = 32.4), as it must.
- **Inputs re-read in the TeX.**
  - BD: Prop. 'cond:CC' and its proof (l.1427–1450). It uses only Momose's Lemmas 3–5, with no condition at 2 or 3 and
    any residue degree f.
  - BD: Def. Momose Type 2 (l.858), Thms 'main-brief' and 'type_1_brief' ("Let k be a number field", l.251–263),
    Cor. 'type_2_not_momose' (l.887–910), the signature Type 3 → Momose Type 3 step (l.840–846), and the GRH remark (l.1457).
  - LV: Thm 1 (l.186–212), the case (2) proof (l.2655–2690), and the GRH remarks (l.2721, 3044–3068).
  - Allen–Genao: every GRH hit (l.50, 75, 77, 276, 283).
  - Genao: every GRH and 'Thm_LV' hit (l.124–182, 314–341, 507) and the references inside §477–580.
- **Galois closures.** PARI (galclosure.gp): ℚ(2^{1/3}) has n_M = 6 and Δ_M = −34992; x³ − x² − 2x − 20 has n_M = 6 and
  Δ_M = −493114979, with Δ_k = −108 and −2891.

## Review
One `reviewer` subagent saw S1–S5 pasted with exact statements (logs/review2.txt). Lemma 2(k), Lemma 2′, Theorem 3,
Corollary 4 and the Allen–Genao consequence are all CORRECT. Overall: ROUTINE relative to the inputs ("a valid
elementary two-prime argument followed by direct classification consequences").

## Novelty
Searched 2026-11-20 (logs/novelty2_*.txt).
- **arXiv listings** (7 categories, 2026-10 and 2026-11), `lst.py 'isogen|Momose|Larson|Vaintrob'`: 7 title hits, none
  on isogeny primes over number fields, Momose types or Larson–Vaintrob.
- **Zenodo.** "isogeny primes": 0. "Larson-Vaintrob": 0.
- **arXiv API.**
  - all:"isogeny primes": 0.
  - Momose AND "Type 2": one unrelated hit.
  - Larson AND Vaintrob AND unconditional: 0.
  - "isogenies of prime degree" AND unconditional: only BD (2203.06009) and 1705.10401 (torsion).
- **Our files** (rg 'Larson|Vaintrob|Momose|isogeny prime|IsogPrimeDeg|Banwait'). Only Harvest4 item 5 (planned, not done),
  this directory, and the COROLLARIES/REFUTATIONS lines for REPORT.md.
- **Assessment.** The finiteness of IsogPrimeDeg(k) without GRH, for every k without an imaginary quadratic HCF, was
  stated only under GRH by LV and by BD (2025). Here it is unconditional, with explicit bounds. The proof is a routine
  combination of chebexp with BD's unconditional Condition CC. The two-prime and cube-root step is the only new
  ingredient. Candidate for COROLLARIES (the orchestrator decides).

## Caveat
- **Inputs.**
  - chebexp: Cureton–Thorner's Lean-checked 31/32, paper unrefereed; chebexp audited.
  - BD's unconditional parts are cited, not re-run: CC, MMIB(k), TypeOneBound(k), the non-Galois Lemmas 1–2 (Props.
    'momose_1_non_galois_md' and 'momose_lemma_2_non_galois'), and the signature → Momose Type 2 and 3 steps.
  - LV Thm 1 as stated.
- **"Effectively computable".** Theorem 3 gives an explicit finite superset of IsogPrimeDeg(k): divisors of BD's integers,
  plus primes < e^{y*(k)}. Deciding which p below 10^70–10^95 really occur, via k-points on X₀(p), is not addressed. BD's
  exact determinations (their Thm 1.9 etc.) remain GRH-conditional for Type 2 primes between their GRH bound and e^{y*}.
- **Size.** e^{y*} is 10^70–10^95, against 10^11–10^14 under GRH. Lemma 2′'s closed form (3000n_M) is crude. The
  certified y* is the usable value.
- **Genao 2405.05507** is not claimed beyond its l.507 step. Its GRH theorems also rest on Thm 'BaseChangeImprovement' /
  [Gen24, Thm 2], not checked.
- **Allen–Genao.** Corollary 5 bounds p by S_F ∪ {p : log p < y*(F)}, which is larger than the GRH statement's S_F.
  The exclusion of LV case (1) uses LV's full conclusion (a CM curve over F with CM field inside F), not merely "F contains
  a CM field". (Repair after the orchestrator's reviewer, 20 Nov.)
- **S_K.** Corollary 4 keeps LV's finite set S_K. Its explicit size is GRH-conditional in LV (Thm 'effectivethm'),
  unconditional per their Remark 'besteff'. Theorem 3 does not need it.

RH relevance: N11

## Files
- `REPORT2.md` (this file), `DERIVATION2.md` (§1 Lemma 2(k) and 2′, §2 Theorem 3 and Corollary 4, §3 Allen–Genao,
  Genao and LV's S_K, §4 the Harvest4 route).
- `cert2.py`, the Arb certificate.
- `galclosure.gp`, PARI Galois-closure discriminants.
- `logs/`: cert2_output.txt, cert2_time.txt, galclosure.txt, novelty2_lst.txt, novelty2_zenodo_lv.txt, novelty2_api.txt,
  review2.txt.
