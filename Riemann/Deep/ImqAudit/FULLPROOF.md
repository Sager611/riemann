# ImqAudit FULLPROOF: the chain of class-number-10000 with source line numbers

Notation: d = |disc K|, D = -d fundamental, X = floor(sqrt(d/4)), chi = (D|.), h = h(D). Files: I2 = Riemann/Deep/ImqClass2,
LW = Riemann/Deep/Littlewood78. Line numbers refer to the files as of 19 Nov 2026 (isieve2.c sha256 checked against
I2/sieve_logs/build.sha, see REPORT.md).

## Step 1. h(D) >= M(d) := sum_{a <= X} r(a) for fundamental D < -4 (isieve2.c l.4-5)

1a. r(a) := #{b mod 2a : b^2 = D mod 4a} is multiplicative with r(p^k) = 1 + chi(p) for odd p not dividing D
    (2 split, 0 inert; any k >= 1), r(p) = 1, r(p^k) = 0 (k >= 2) for p | D (p || D as D fundamental). p = 2:
    D = 1 mod 8: r(2^k) = 2; D = 5 mod 8: 0; D = 0 mod 4 (D/4 = 2, 3 mod 4): r(2) = 1, r(4) = 0. Hensel/CRT, standard.
    This is exactly isieve2.c's r (l.5) and the code's pattern weights (l.96-103: split weight 2 for every power,
    ramified weight 1 first power only) and dfs_M (l.120-135).
1b. For a <= sqrt(d/4) and b in (-a, a] with b^2 = D mod 4a, put c = (b^2 - D)/(4a). Then 4ac = b^2 + d >= d >= 4a^2,
    so c >= a; c = a forces b = 0 and d = 4a^2, impossible for fundamental D < -4 (D/4 = -a^2 not allowed).
    b = -a is excluded by b in (-a, a]. So (a, b, c) is reduced. gcd(a, b, c) = f gives f^2 | D, so f = 1 (D
    fundamental: odd f^2 never divides D; f = 2 would give D = 0 mod 16): the form is primitive.
1c. Distinct (a, b) give distinct reduced forms, and every proper class of primitive positive definite forms of
    discriminant D < -4 contains exactly one reduced form. Hence h(D) >= #{(a, b)} = sum_{a <= X} r(a) = M(d).
    (Equivalently M(d) = #{reduced forms of D with a <= sqrt(d/4)}, so h - M = #{reduced forms with
    sqrt(d/4) < a <= sqrt(d/3)} >= 0; ImqAudit/hcount.c computes both sides independently.)
    X = isqrt(floor(d/4)) = floor(sqrt(d/4)) (isieve2.c l.68): n^2 <= x iff n^2 <= floor(x).
1d. Candidates: D = -d fundamental iff d = 3 mod 4 squarefree (d = 3, 7, 11, 15 mod 16) or d = 4m, m squarefree,
    -m = 2, 3 mod 4, i.e. m = 1, 2 mod 4 (d = 4, 8 mod 16). So the six classes cover every fundamental d; the
    other classes (d = 0, 1, 2, 5, 6, 9, 10, 12, 13, 14 mod 16) contain no fundamental d. Non-fundamental d in the six
    classes are swept too; for them the sieve's verdict is irrelevant, and tier 3 drops them (fundamental(), l.148-151:
    u = d (d odd) or d/4 (d even, then u = 1, 2 mod 4), squarefree iff no odd p^2 | u, trial division by odd p^2 <= u).

## Step 2. The sieve bound (isieve2.c l.4-18, l.363-405) and its implementation

2a. a <= Xs is split uniquely as a = m n, m 13-smooth, n 13-rough; r(a) = r(m) r(n); r(m) = 0 unless m in G_c
    (split p <= 13 to any power, ramified p <= 13 to the first power; built in build_patterns l.89-112 with weight
    2 per split prime and 1 per ramified prime; prefix sums Fp l.113-116 = F_c(t) = sum_{m in G_c, m <= t} r(m)).
    The class c = d mod QS, QS = 240240 = 16*3*5*7*11*13, fixes chi(2) (l.325: d = 3, 11 mod 16 inert, 7, 15 split,
    4, 8 ramified; others skipped) and chi(p), p = 3..13 (l.327: kron(c, p) = (-c|p) = (-d|p)).
2b. Family of distinct 13-rough n with r(n) > 0 (p, q split primes >= 17): n = 1; n = p^k (17 <= p <= Xs/17,
    k >= 1, p^k <= Xs; l.381); n = pq (p < q, pq <= Xs, so q <= Xs/17; l.334-336, prefix by product l.375,
    sanity check l.373); n = p (Xs/17 < p <= min(Xs, 250000), p >= 17 as pmed starts at 17, l.284-287; l.401).
    (Without p >= 17 the display would be false for tiny Xs, e.g. d = 23; production Xs >= 83,666.) Distinct n, distinct a = mn, so
    M(d) >= sum_n r(n) F_c(floor(Xs/n)) >= sum_n r(n)' F_c(...)  (dropping ramified p > 13 and every other n only
    lowers the sum). Xs = X(segment start) <= X(d) for every d of the segment (l.363-365), F_c nondecreasing.
2c. Rounding round2 (l.173-174): cap 32767, then the largest 2^k or 3*2^(k-1) <= w: always <= w. Phase A runs: each
    prime gets round2(2 sum_k F_c(Xs/p^k)) of its own (runs merge only equal rounded weights, l.382-384). Pair runs
    (l.386-399): t = Xs/(pq) runs from tmax = Xs/323 down to 1; pairs with Xs/(pq) = t are indices
    [posT[t+1], posT[t]); a run [posT[tb+1], posT[t+1]) is emitted when round2(4F(t)) changes, with the weight of
    t' in [t+1, tb] -- each pair gets round2(4 F_c(Xs/(pq))). Phase B: round2(2 F_c(Xs/p)) per prime (l.401-405).
2d. Counting: T_p bit x = [-QS x nonzero square mod p] (l.293-307), built directly for x < p + 64VW + 64, so every
    window [o, o + 64VW + 64), o < p, is periodic-correct. Lane j of class c has x = c QS^{-1} + j mod p (base_init
    l.223, goff l.222; QS^{-1} checked l.297). Hence bit = 1 iff (-d|p) = 1 (split); ramified p (x = 0) gives 0.
    Self-check l.316-320 (400,000 random (d, p) against Jacobi symbols). Harley-Seal CSA (l.189): a+b+c = 2h + l,
    invariant o + 2t + 4f + 8e + 16 c16 = number of set items among <= 256 per chunk (c16 is a 5-bit counter, max 16).
    addsh (l.154-161): exact ripple add of count << b over 16 planes; any bit at plane >= 16 or final carry sets ovf
    (true sum >= 2^16 > N + 1). geq (l.162-165): MSB-first exact comparison acc > needv - 1. Items are masked by the
    live lanes only; a certified lane stays certified (its acc no longer changes). need = N + 1 - F_c(Xs) (l.377):
    need <= 0 certifies by F_c alone (l.414).
2e. Lanes: j in [ceil((seg_lo - c)/QS), floor((seg_hi - c)/QS)] (l.408-409) = exactly the d = c mod QS in
    [seg_lo, seg_hi]; valid-lane mask GVM l.420-425; segments [LO + i LSEG, min(.. + LSEG - 1, HI)] (l.363-364), no gap.
    Uncertified lanes: tier 3 = exact M(d) by depth-first enumeration (l.120-147; X = Xof(d); every a <= X with
    r(a) > 0 visited once, primes ascending, powers of split primes), capped at N + 1; M >= N + 1 passes; else
    non-fundamental lanes are dropped (fundamental() l.148-151) and the rest go to tier 4 (l.476-480).
2f. Production logs (sieve_logs/p_*.log, 38 chunks): every header N=10000, pairs on, PT=250000; chunks consecutive
    from LO = 28,000,000,000 to HI = 762,700,000,000 (no gap, no overlap); for each chunk lanes = the number of
    d = 3, 7, 11, 15, 4, 8 mod 16 in [LO, HI] (own count, ImqAudit/logcheck.py) = certF + cert1 + t3, t3 = t3pass +
    nonfund + t4, bad 0, no FAIL/MISMATCH line; total 275,512,500,000; t4 = 33 (chunk x5 only). isieve2.c sha256
    be745424... = sieve_logs/build.sha (source of the production binary).


## Step 3. Tier 3 and tier 4

3a. Tier 3 = dfs_M (Step 2e) is the exact M(d), capped at N + 1; symbols for p <= 13 and p > PT by kron (Jacobi,
    l.75-77: p = 2 gives 0 for d even, +1 for d = 7 mod 8, -1 for d = 3 mod 8, matching Step 1a), table primes by
    kron_tab (= kron, self-checked l.316-320). M >= N + 1 → h >= M > N. Fundamental d with M <= N → tier 4.
3b. Tier 4 = t4.gp = classno_upto(d, N) of cln-imq-2k v1.4 qsieve/postproc.gp (71 lines, sha256 17834364db42...,
    fetched 19 Nov 2026; ImqAudit read l.25-48): starts from the identity, and for each prime p <= floor(sqrt(d/3))
    with (D|p) != -1 adjoins g = qfbred(qfbprimeform(D, p)) by appending the disjoint cosets G g^k until g^m in G
    (disjointness asserted); stops with N + 1 as soon as #G > N, else returns #G. Unconditional: every class has a
    reduced form with a <= sqrt(d/3), whose primitive ideal has prime factors of norm p | a, so these prime classes
    generate Cl(D); qfbcomp/qfbred return the unique normalized reduced form, so (a, b) keys distinct classes.
    Returns min(h, N + 1) exactly. Used only on fundamental d > 4 (tier 4 is fundamental by construction).
3c. ImqAudit recomputed the 33 tier-4 d of chunk x5 independently (hcount.c, own reduced-form count): M agrees
    with x5_t4.txt for all 33; h <= 10^4 only for d = 28,154,696,707 (h = 9819), the other 32 have h >= 10,050.

## Step 4. Range (i): d >= e^27.36 (threshold.py, LW/cert.py)

4a. h(D) = sqrt(d) L(1, chi)/pi for d > 4 (w = 2, Dirichlet), chi primitive mod d.
4b. Cells [l, l + 0.01] of l = log d from 27.36 to 27.53 (cert.py class_threshold l.246-268): certified when
    exp(l/2) exp(F(l + 0.01, P))/pi > 10^4 in arb balls (l.259-260). F (l.42-61) depends on l only through
    W = ... + (l - log pi)/2 (l.54-56) with val = ... - kp W, kp = k0/(1 - k1) > 0 (k1 < 1 checked, l.47): F is
    decreasing in l, so log L(1, chi) >= F(l, P) >= F(l + 0.01, P) on the cell, and sqrt(d) >= exp(l/2).
4c. Beyond l_end = 27.53: claim (A) |L(1, chi)| >= 0.11/log log d (log log d >= 3), and exp(l/2)/log l is increasing
    for l >= 27.53 (derivative 1/2 - 1/(l log l) > 0): h > 10^4 (l.250-253). Float check: 0.11 e^(l/2)/(pi log l) =
    10,041 at l = 27.53 and 9,992 at 27.52, so l_end = 27.53 is the first hundredth that works; 2753 - 2736 = 17 cells.
4d. ImqAudit reran threshold.py 10000: "10000 27.36 7.627e11 27.53" (identical); e^27.36 = 762,600,393,219.6 <
    762,700,000,000 = sieve top (decimal check). LittlewoodAudit (PASSED) refereed F and re-certified the 10^4
    threshold with its own cells.

## Step 5. Range (iii): d < 111 * 2^28 = 29,796,335,616 (mjx.c)

5a. Input: Mosunov-Jacobson, Math. Comp. 85 (2016) 1983-2009, all fundamental |D| < 2^40, unconditional
    (class numbers verified with the Eichler-Selberg trace formula); trusted as published (page Caveat).
5b. mjx.c process_file (l.78-153): reads delta-coded d (d += m * vals[0], l.116) and requires each d to be the next
    fundamental d of its class in the slab (expected e from a sieve of odd p^2, l.41-66, l.127-129), and no
    fundamental d missing at file end (l.149): each file lists exactly the fundamental d of its class and slab;
    invariants multiply to h, divisibility chain (l.119-125); gzip CRC (l.145-146). Classes 3, 7 mod 8 and 4, 8 mod 16
    = all fundamental d. merge_update (l.158-180) merges the four classes in increasing d (strict order checked).
5c. ImqAudit (slabcheck.py on mjx_slabs.log): 444 file lines, slabs k = 0..110 contiguous, all four classes in each,
    first/last d inside each slab; 9,056,999,962 entries; slab 104 (ends 28,185,722,880) has min h 9736; slabs
    105-110 min h >= 10,052.


## Step 6. Assembly (page Claim)

6a. Fundamental d >= e^27.36: h > 10^4 (Step 4). 2.8e10 <= d <= 762,700,000,000 (>= e^27.36): every fundamental d has
    M(d) >= 10^4 + 1 by tier 1 or tier 3 (Steps 2-3), except the 33 tier-4 d, decided exactly (Step 3b, 3c): only
    d = 28,154,696,707 (h = 9819). d < 29,796,335,616 (> 2.8e10): the tables (Step 5). The three ranges cover all d > 0
    with overlaps at both joins; d = 3, 4 (h = 1) are table entries (the M(d) lemma is only used for d >= 2.8e10).
6b. In the overlap [2.8e10, 29,796,335,616) the sieve's tier 4 found exactly the table fields with h <= 10^4 (one:
    28,154,696,707), t4check.py; ImqAudit: slab 104 min h 9736 (d = 27,993,568,387 < 2.8e10), slabs 105-110 min h >= 10,052.
6c. Hence the fields with h <= 10^4 are the table fields with h <= 10^4: 410,810,155 (perh total; = row N = 10^4), the
    largest 28,154,696,707 (h = 9819; hcount exact), every h <= 10^4 with count >= 1.
6d. Inputs not re-derived here: Mosunov-Jacobson tables (published, unconditional); Littlewood78 (paper proof, refereed by
    LittlewoodAudit, PASSED); the Lean-checked half-plane Re s > 0.87338196661320752 (HPAudit).

