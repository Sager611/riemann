NO FLAW FOUND: every step of Banks, arXiv 2610.31402v1, checked against its TeX, Masser 1981 Thm 1 and Stark 1974 (both read from the GDZ scans), with exact and numerical tests. No false step; the proofs use no unproved input.
# Audit: W. D. Banks, "Landau–Siegel zeros over a base field and effective Brauer–Siegel" (22 Oct 2026)
TeX lines: thm:quadratic 249, cor:siegel 271, thm:bs 288, thm:main 358, thm:precise 380, lem:generator 598, lem:bias 724,
lem:stark 837, thm:masser 1085, prop:mu 1100, lem:frob 1260, lem:pigeon 1426, prop:main 1548. Details: NOTE.md; scripts check_*.
1. The "OpenAI argument" is OpenAI, "Uniform exclusion of Landau–Siegel zeros" (1 Oct 2026, github.com/openai/math preprints),
   Thm 1: an absolute c > 0 with (1−β) log q ≥ c. Banks never uses this theorem. He re-runs Kerr's form of it (2610.21099:
   pigeonhole + Masser), which we already audited as correct (Riemann/Deep/Siegel). OAI, Green, Barreto, Kominers: background and remarks only.
2. The zero estimate is Masser, Invent. Math. 63 (1981), Thm 1, p. 82: "For real N ≥ 0, D ≥ 1 let P be a non-zero polynomial
   of total degree at most D such that P(γ) = 0 for all γ in Γ(N). Then we have D ≥ (N/n)^μ." Banks's thm:masser states the same thing.
   Every hypothesis holds: the ω_t^S form a Z-basis (1 ∈ S); total degree ≤ weighted degree; the box is Γ(N−1); the bound
   depends only on (n, μ), so it is uniform in F.
3. μ(Γ) = n_K/|S| (prop:mu): we re-derived the Galois-averaging proof and it is correct. About 72k subspaces of Q(√d, ξ_5) and Q(√d, ξ_7): 0 violations.
4. The pigeonhole chain (lem:frob through lem:pigeon) is correct. GP check of θ^{Np} ≡ σθ and of N(Q_{p,σ}) = Np^{|G|/|C_p|}
   (S_3 over Q, and C_3 over Q(√−3)): 2820 checks, 0 failures.
5. lem:count, prop:main and thm:precise are correct. Exact lattice sums match I, (1+(s−1)/H)I and I/Np, and Kerr's
   384 / 6144 / 1536 are reproduced. A toy run of inequality (pigeon) in Q(√5, √−3) gives lhs/rhs ≤ 0.30 (< 1, consistent with Masser).
6. lem:generator is correct. GP on 750 quadratic extensions F of Q(√10), Q(√−5) and a cubic field: max|ιτ|/d_F^{1/n_k} ≤ 0.40.
7. thm:quadratic and cor:siegel are correct. The arithmetic is exact in 302,952 (ε, n_k, |G|) cases. lem:bias is Stark's (8)–(10)
   with an absolute constant. The finitely many small fields are handled effectively (Hermite, mean value theorem, Rademacher).
8. thm:bs is correct. Stark's Lemma 3, Theorem 3 and Lemma 4 are quoted faithfully (scan pp. 139–145).
Effectivity (task 2): it is gained in thm:precise, an unconditional upper bound for the log-mass of inert primes, built from
Frobenius congruences, a pigeonhole and Masser's uniform zero estimate. β enters only in lem:bias. No second character and
no Deuring–Heilbronn argument is needed. For k = Q: (1−β) log q ≥ c := min(1/(16 c_1 C), finite check), C = max(32, 8c_1),
for every real primitive χ. This c is effective, but nobody has computed it.
Caveats (not flaws): c(k) is not explicit. ℓ = 3 would already suffice for cor:siegel. Kominers' bound and its base-field
form (rem:compare) are claims without proof, in remarks only. Line 891, "ζ_k(1) ≠ 0", should read "(s−1)ζ_k(s) ≠ 0 at s = 1".
For k = Q the result is known (OAI, Kerr).
Action: no change to REFUTATIONS.md or DISCOVERIES.md (the paper is correct, and the result is not ours).
