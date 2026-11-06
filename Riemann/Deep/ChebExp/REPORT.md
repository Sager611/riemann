RECORD: first numeric polylogarithmic least prime in every Chebotarev class. For every Galois L/K (L ≠ ℚ) and every class C there is an unramified degree-one prime 𝔭 of K with Frob_𝔭 = C and N𝔭 ≤ (0.484·log 3Δ_L)^32, from Cureton–Thorner's Lean-checked half-plane Re s > 31/32. That constant is optimal to 0.07 %, since ℚ(√−3) forces ≥ 0.48365; it drops to 0.094 once log Δ_L ≥ 10⁴. Conditional on Liu's B_2: N𝔭 ≤ (0.62·log 3Δ_L)^{7.8978} (forced ≥ 0.58228; 0.344 once log Δ_L ≥ 10⁴). Both are certified in Arb.

# ChebExp: explicit least prime in Chebotarev classes from a zero-free half-plane (node N11)

Agent chebexp, 2026-11-06, 04:12–06:45 UTC. Brief: `.scratch/orca/chebexp-brief.md`. Full proof: `DERIVATION.md`.

## Claim
L/K is Galois with L ≠ ℚ, G = Gal(L/K), C ⊂ G a conjugacy class, Δ_L = |disc L|. Then there is a prime 𝔭 of K,
unramified in L, with N𝔭 = p prime, Frob_𝔭 = C and:

| input (only ζ_L of the top field is used) | bound, every L ≠ ℚ | for log Δ_L ≥ 10 | ≥ 100 | ≥ 1000 | ≥ 10⁴ |
|---|---|---|---|---|---|
| **A.** ζ_L ≠ 0 on Re s > 31/32 (Cureton–Thorner, Lean-checked) | N𝔭 ≤ (0.484·log 3Δ_L)^32 | 0.29 | 0.151 | 0.108 | 0.094 |
| **B.** ζ_L ≠ 0 on Re s > B_2 = 0.8733819666 (Liu, unrefereed) | N𝔭 ≤ (0.62·log 3Δ_L)^{1/(1−B_2)}, 1/(1−B_2) = 7.89777 < 7.90 | 0.595 | 0.426 | 0.366 | 0.344 |

Each table entry replaces the constant 0.484 (respectively 0.62) in that range.

- **Optimality at small fields.** For ℚ(√−3)/ℚ and C = {1}, the least split prime is 7. That forces c ≥ 7^{1/32}/log 9
  = 0.48365 and c' ≥ 7^{1−B_2}/log 9 = 0.58228 in any bound of this form. A is optimal to 0.07 %; B is within 6.5 %.
- **Limit of the method.** Floating point, not certified: as Δ_L → ∞ the constant tends to 2.596·δ, where δ = 1/32
  for A and δ = 1 − B_2 for B. That is 0.0811 (A) and 0.3287 (B). The table reaches 3.0δ (A) and 2.7δ (B).

## Before
- **Unconditional numeric bounds were all powers of Δ_L.**
  - Kadiri–Wong (arXiv 2105.14181): N𝔭 ≤ Δ_L^{310} for every L ≠ ℚ. This was the best fully numeric bound.
  - Ahn–Kwon: Δ_L^{12577}.
  - Asymptotic, with explicit exponent only: Zaman (1508.00287) ≪ Δ_L^{40}; Kadiri–Ng–Wong (1902.08640) ≪ Δ_L^{16};
    Thorner–Zaman (1604.01750).
  - Bound A lies below Δ_L^{16} for every Δ_L ≥ 3: it is 7.2 at Δ_L = 3, and its log grows more slowly.
- **After the 2026 half-planes, the polylog shape was known without constants.**
  - Frei–Sofos–Widmer, arXiv 2610.18491v2 Cor. 1.2: c(δ)(log Δ_L)^{1/δ}(log log Δ_L)^{2/δ}, with c(δ) "effectively
    computable", never computed. Only v1 and v2 exist.
  - Hu, Hexagon 2610.00216 (9.10): ≪ (log eΔ_L + n_L)^8, constant unspecified.
  - Our Riemann/Deep/ArtinHalfPlane/NOTE.md §3.3: (C*·log 3Δ_L)^{1/δ}, C* unspecified.
  - COROLLARIES (R317): Chebotarev with error O_L(x^{B_2}log²x), not explicit.
- **Under GRH:** Bach–Sorenson 1996, (4 log Δ_L + 2.5n_L + 5)².

## Method
- **Deuring reduction.** Fix g ∈ C, H = ⟨g⟩, E = L^H. The prime sum over E is weighted by
  θ = |H|^{−1}Σ_{χ∈Ĥ} χ̄(g)χ(·) ∈ [0, 1] (exact also at ramified primes). Each χ is a primitive Hecke character of E,
  and ζ_L = ∏_χ L(s, χ).
- **Kernel.** k(t) = t^{−s*}B_3(log 1/t), with B_3 the quadratic B-spline of knot width h. Its Mellin transform is
  Φ(s) = ((1 − e^{−h(s−s*)})/(s − s*))³.
- **Explicit formula.** The contour goes to −∞, so the trivial zeros appear as an explicit, summed series.
- **The zero sum needs no zero counting.**
  - The Hadamard identity for ζ_L at one real point σ_0 > 1 gives
    Σ_ρ Re 1/(σ_0 − ρ) ≤ ½log Δ_L + 1/σ_0 + 1/(σ_0 − 1) − n_L g(σ_0). The Γ-factor term is a saving.
  - The ratio |x^{s−1+δ}Φ(s)| / Re(1/(σ_0 − s)) is subharmonic on the strip δ ≤ Re s ≤ 1 − δ. So its sup is attained
    on the two boundary lines: a 1-D sup, certified by branch and bound.
- **Bad terms.**
  - Ramified primes, using v_𝔭(d_{L/K}) ≥ |G|/2 and a Riemann-sum bound per prime:
    ≤ 2 log Δ_L·min(‖k‖ log x, ‖k‖(U + log 2), Φ(0) + ‖k‖ log 2)/log 2.
  - Prime powers and higher-degree primes: ≤ n_L·P_k(x). P_k is summed exactly for x ≤ 2·10¹². Above that, a
    partial-summation bound with θ(y) < 1.01624y is used.
- **Uniformity in n_L.** By Minkowski, n_L ≤ log Δ_L/log(π/2). The inequality is affine in n_L, so both ends are
  checked.
- **Covering all fields.**
  - Analytic part: chains of intervals of log Δ_L, each with its own kernel and a fixed
    x = (c(Λ_a + log 3))^{1/δ}, plus a monotone tail.
  - Small fields directly. For A: Δ_L ≤ 12, where Minkowski leaves 8 quadratic fields. For B: Δ_L ≤ 257, where
    Minkowski gives degree ≤ 4 and PARI nflist lists all 200 fields; 563 (L, K, C) cases with their least primes.

## Certificate (one command)
```
cd Riemann/Deep/ChebExp && uv run --no-project --with python-flint python3 cert.py chains.json > logs/cert_output.txt
```
- Ends with `ALL CERTIFIED (104 s)`; one core, peak RSS 38 MB.
- Covers 10 chains (405 kernels) and both direct ranges (the B range is read from logs/direct_range.txt, made by
  `gp -q direct.gp`). The smallest margin is LHS/Φ(1) ≥ 0.0030.
- Every constant is an Arb enclosure used in the safe direction. Φ(1) and g(σ_0) enter as lower bounds; every error
  constant as an upper bound.
- `design.py` (floating point) produced `chains.json`; it is not part of the proof.

## Checks
- **Least primes.** `sanity.py`, logs/sanity_output.txt. All are below both bounds.
  - Every class of ℚ(ζ_q)/ℚ, 3 ≤ q ≤ 2000 (1,013,926 classes). The largest empirical constant is 0.4837, at q = 3,
    p = 7: the ℚ(√−3) extremal case.
  - The S_3 closures of all 508 cubic fields with |d| ≤ 3000 (PARI nflist, `s3_least.gp`), 1524 classes.
    The largest empirical constant is 0.108.
- **Normalisation.** For ℚ(ζ_q)/ℚ, q = 5, 7, 11, 13, x = 10⁷: |H|S(x)/(xΦ(1)) = 1.00001, 1.00002, 0.99904, 1.00134.
  The errors lie below both the proven-shape bound (δ = 1/32) and the GRH-shape bound (δ = 1/2).
- **Hadamard identity for ζ at σ = 2:** 2000 zeros plus the tail give 0.0690665; the identity gives 0.0690662.
- **Independent re-evaluation.** `xcheck.py` (mpmath/numpy) reproduces the first-interval margins of both uniform
  chains (logs/xcheck_output.txt) and confirms the maximum principle on a dense grid.
- **Self-review.** DERIVATION.md §7 goes through the Deuring step, the primes ramified in L/E and in L/K, the absence
  of truncations, uniformity in n_L, rounding directions, coverage, and the forced lower bound.
- **Bugs found and fixed during the run.**
  - Arb cancellation in the spline overestimated ‖k‖ about 5×. F is now evaluated piecewise.
  - The crude prime-power bound ‖k‖(ψ − θ) failed for wide kernels. It was replaced by the partial-summation bound.
  - The float log 13 is 2.6·10⁻¹⁸ above the true value. The chain now starts at 2.5649, asserted < log 13 in Arb.

## Novelty
Searched 2026-11-06. Not found anywhere:
- **arXiv listings.** math.NT, Oct–Nov 2026, `lst.py 'Chebotarev|prime ideal|least prime|Frobenius'`: 9 title hits.
  Only 2610.36464 and 2610.42022 are on Chebotarev, and both are conditional (Mumford–Tate plus standard analytic
  conjectures; motivic RH).
- **arXiv API.** "least prime ideal" AND Chebotarev; Chebotarev AND "least prime" AND explicit; au:Thorner AND
  au:Zaman. All hits are power bounds: 1508.00287, 1902.08640, 1604.01750, 1803.02823, 2105.14181.
- **Frei–Sofos–Widmer.** 2610.18491: only v1 and v2.
- **Hexagon.** All 716 works via the API, titles and abstracts. "Chebotarev" appears only in 2610.00216 (Hu:
  "an effective Chebotarev theorem", non-explicit), 2610.00451 and 2610.00274 (unrelated).
- **Zenodo.** "Chebotarev": 145 records, none on least primes.
- **HAL.** "Chebotarev" 2025–2026: 10 records, none relevant.
- **Our files.** DISCOVERIES, COROLLARIES, FORMALIZATIONS, REFUTATIONS, ArtinHalfPlane, ExplProg and the NOVELTY files
  give only the non-explicit forms above. ExplProg is the Dirichlet/Linnik analogue in the q-aspect.

Assessment: the qualitative statement is known (FSW Cor. 1.2 with an extra log log factor; Hu (9.10); our NOTE §3.3).
What is new:
- the first numeric constants for the polylogarithmic bound;
- a uniform constant for A that is optimal to 0.07 %;
- certified constants for large discriminants.

## Caveat
- **Inputs.** A rests on Cureton–Thorner's 31/32 (Lean-checked, paper unrefereed). B rests on Liu's B_2 (unrefereed;
  our audit found it sound given his 2610.26835v2).
- **Classical facts cited, not re-proved.** The qualitative contour estimates (zero counts in unit windows, L'/L at
  good heights), the Hadamard identity, Artin reciprocity with Hecke conductors, Minkowski, and Rosser–Schoenfeld's
  θ(y) < 1.01624y. None of the qualitative estimates contributes a constant.
- **PARI.** B's direct range trusts the completeness of PARI nflist for degree ≤ 4 and |disc| ≤ 257. The field counts
  match the classical tables.
- **Review.** Two quick passes by a `reviewer` subagent, each on pasted excerpts (§3.2–4.3; §3.1, §5 and the
  cert code), found only presentation gaps, now fixed (DERIVATION §7, end). A full independent review is still
  pending.
- **Not optimal for large fields.** The kernel family (tilted quadratic B-splines, one Hadamard point) is restricted,
  so the large-discriminant constants are not optimal.
- **Size of the numbers.** Exponent 32 makes A huge for real fields: Δ_L = 10⁶ gives about 10^{27.5}. B gives about
  10^{7.6} there.

RH relevance: N11

## Files
- `REPORT.md` (this file), `DERIVATION.md` (every step and constant).
- `cert.py` (certificate), `chains.json` (intervals and kernels), `design.py` (how the chains were found).
- `direct.gp` (B's direct range), `sanity.py`, `s3_least.gp`, `xcheck.py`.
- `logs/`: cert_output.txt, cert_time.txt, direct_range.txt, sanity_output.txt, s3_least.txt, xcheck_output.txt.
