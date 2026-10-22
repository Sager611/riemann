# Working notes: audit of arXiv 2609.22336 (Deng), 2026-10-21

Sources: the v2 TeX e-print (main, intro, local, end: 1281 lines, all read) and v1 (main.tex, 1130 lines, read for differences). For the cited results, Cullinan–Hajir's preprint PDF (faculty.bard.edu/cullinan/papers/legendre.pdf): Thm 1.9(b) says Gal(L_{2n+δ}) ≃ S₂≀S_n for half-degree n ≤ 60 (Magma), and Cor 3.4(b) says L_m is irreducible for p−4 ≤ m ≤ p+3. All downloads were deleted after use.

## v1 to v2
Same reduction and lower bound. v1's large-prime estimate was weaker: log𝓡 < k√(2n)log(2n) + 2m log2·h_{2k−1}. It needed n ≥ 10⁵ and Groth's 2025 thesis computation below that. v2 replaces it with v_p(R) = 0 for p > m (Lemmas 3.6–3.7). Its threshold becomes n ≥ 256, with only Cullinan–Hajir needed below.

## Where a flaw would have to sit, and why none does
- Lower bound. Lemma 5.5 is uniform over any k points of [0,1] (Hadamard on the matrix (T_{j−1}(2t_i−1))). |H'(t_i)| is bounded root by root. No assumption about spreading or interlacing.
- 2-adic input: only Wahab's polygon. The digit groups are exactly the Q₂-irreducible factors. So A(1) is exactly 2^{2k−s₂(k)}, a > A(1), and cross-group pairs have unequal v₂(1−t). Hence v₂(t_i−u_j) = min(...) < 2.
- Odd primes. The per-root cap is an identity, not an estimate: v_p(P_n'(α)) = −v_p(U_n(α)) when v_pα ≥ 0. It is attained (776 of 6023 factor checks).
- Large primes (essential). Integrality of F(0) and lc F, together with v_p(M_n(0)), v_p(lc M_n) ≤ 1, forces all p-adically small roots into one factor and all large ones into one factor. Unit roots are simple mod p.
- Mechanism. The archimedean size plus the 2-adic cap leave 𝓡 ≳ (2^{v₂(n!)})^k ≈ 4^{km}. Odd primes ≤ m with capped exponents give at most about e^{km}.

## Reproduction (PARI/GP; each run < 2 min, < 100 MB)
- Per-root test: for each f in factorpadic(c·P_n, p, 60)[,1], with c the denominator of P_n:
  va = (v(f(0)) − v(lc f))/deg f, vd = (v(polresultant(f, c·P_n')) − (n−1)·v(lc f))/deg f, loc = vd − (n−2)·min(va,0).
  Require: loc ≤ logint(2n−1,p); vd ≤ logint(n,p) if va ≥ 0; loc ≤ L_p + va if va < 0; vd = 0 if p > m and va = 0. (f = x gives va = ∞ and loc = vd.)
- 2-adic test: factorpadic(M_n, 2, 80). Check each factor is even, v₂(f(1)) − v₂(lc f) = deg f − 1, and the degree multiset equals 2×(digits of m). For every split, v₂(polresultant(A,B)) ≤ e(2k−s₂(k)) after removing the lc terms.
- Lemma 5.4: take roots from polrootsreal and evaluate P_n, P_n' by the three-term recurrence at realprecision 400. Monomial evaluation at 80 digits cancels catastrophically and gives a false 0.04 at n = 270.
- Margins: LB(n,k) = (2e−s₂(e))log2 + log γ_n + 2(k−1)²/k·log2 − log k, against log D_n = Σ_{odd p ≤ m} logint(2n−1,p)·log p.
- Small n: polisirreducible on P_n (n even) and P_n/x (n odd), n ≤ 255, takes 7 s.

## Operational notes
- `perl -e 'alarm T; exec @ARGV' gp …` does NOT stop gp: gp traps SIGALRM. Use fork plus kill -9 (as in the run.sh wrapper used here).
- Footprint: gp peak resident 81 MB; total CPU about 6 min.

## Verdict
No refutation. The proof appears complete and correct for n ≥ 256 and, with exact D_n, for every n ≥ 64. Small n follow from Cullinan–Hajir and our PARI check. It is not ours and is unrefereed, so it needs expert review.
