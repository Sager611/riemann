RECORD: Σ_{χ mod q} N(σ, T, χ) ≪ (qT)^{4(1−σ)/(6σ−3)+ε} and Σ_{q≤Q} Σ*_χ N(σ, T, χ) ≪ (Q²T)^{4(1−σ)/(6σ−3)+ε} for 5/6 < σ ≤ 7/8, below the old record 3/(3σ − 1) (Huxley 1975, as stated by Chen–Gupta–Li 2025; Jutila 1977 lists only 2 for families); unconditional (Petrow–Young's hybrid Weyl bound + Jutila's family large-value theorem + Ivić's 3/(2σ) for the zeta part). The Lean-checked half-plane does not improve these uniform exponents; it gives a longer family Montgomery range and a better large-sieve exponent only when Q^{0.7..0.8} ≤ T ≤ Q^{13.8..14.7}.

# DensityQ report (2026-11-05)

Proofs: [PROOF.md](PROOF.md). Certificate: [verify.py](verify.py) (exact rationals, standard library, ALL PASS, 19 s,
29 MB). Paper proofs, not refereed, no Lean.

## In words
- Zeros of Dirichlet L-functions in families are counted by large values of Dirichlet polynomials. Two kinds of
  polynomials matter: long ones (handled by large-value theorems, where the half-plane helps) and L-function partial
  sums of length about X^{1/(6σ−3)} (X = qT or Q²T), whose large values are controlled only by pointwise bounds for
  L(1/2 + it, χ). The second kind caps every uniform family exponent at 4/(6σ − 3) on [0.84, B), with or without the
  half-plane (Proposition B).
- The best stated family bound here is 3/(3σ − 1) (Huxley 1975's large-value estimate, as stated by Chen–Gupta–Li
  2025; Jutila 1977 lists only the density hypothesis 2 for families). Petrow–Young's L(1/2 + it, χ) ≪
  (q(1 + |t|))^{1/6+ε} (all q, Duke 2023) excludes the partial sums, Jutila's 1977 family large-value theorem handles
  the rest, and the principal characters (that is, ζ) are covered by Ivić's 3/(2σ): A ≤ 4/(6σ − 3) for both families on
  (5/6, 7/8] (Theorem C), below 3/(3σ − 1) exactly for σ > 5/6. This is the record; the half-plane is not needed.
- The half-plane does give (Theorem A) the Montgomery conjecture for family polynomials with X ≤ N^τ, τ < τ_B(σ) =
  (2σ − 2 + B)/(B − 1/2), longer than Jutila's 3σ − 1 for σ > 0.856092, and (Theorem D) the zeta-type bound
  3(B − 1/2)/((5/2 − B)σ − (5/2 − 2B)) for the large-sieve family when T is a moderate power of Q.

## Table
| σ range | family | old record (source) | new bound | derivation |
|---|---|---|---|---|
| (5/6, 7/8] ⊃ [0.84, B) | one modulus, A_q | 3/(3σ − 1) (Huxley 1975 as stated by Chen–Gupta–Li 2507.08296, eq. zeeshmh; Jutila 1977 Cor. (1.5): 2) | 4/(6σ − 3), unconditional | Thm C: Petrow–Young removes class-II zeros for Y = X^{1/(6σ−3)+} (non-principal primitive χ); class I by Jutila (1.2), k = 2, and squares/cubes; ζ (principal characters) by Ivić's 3/(2σ) |
| same | large sieve, A* | 3/(3σ − 1) (same Huxley 1975 estimate; Jutila 1977 Cor. (1.6): 2) | 4/(6σ − 3), unconditional | Thm C with Jutila (1.3); ζ by Ivić's 3/(2σ) |
| [0.84, B) | both, with the half-plane | — | 4/(6σ − 3), no further gain | Prop. B: the cut 6σ − 3 is not moved by the half-plane; with any Montgomery range τ_M ≤ 8σ − 5 the scheme gives exactly 4/(6σ − 3) |
| (0.856092, B) | large sieve, λ = log(QT)/log(Q²T) ∈ [λ_a, λ_b] (T from Q^{0.69..0.79} to Q^{13.8..14.7}; the bound itself holds from T ≥ Q^{0.15..0.18}) | 3/(3σ − 1) (stated); 3/(2σ) (same method without the half-plane, unpublished) | 3(B − 1/2)/((5/2 − B)σ − (5/2 − 2B)) | Thm D: Thm A + t-subdivision + Tao–Trudgian–Yang Cor. zero-large-cor3, class II absent for λ ≤ λ_b; ζ by DensityB |
| (0.856092, 0.900522) | large values, both families | Montgomery range τ ≤ 3σ − 1 (Jutila 1977 (1.2)/(1.3), k = 2) | τ < τ_B(σ) | Thm A: Halász–Montgomery with off-diagonal character sums ≪ N^{1−B}X^{B−1/2+ε} (reflected half-plane, conductor-uniform) |
| [B, 1) | all | ... | N = 0 | the half-plane |

Values (A): σ = 0.84: 1.96078; 0.85: 1.90476; 0.86: 1.85185; 0.87: 1.80180; B−: 1.78548 (old 3/(3σ − 1): 1.97368,
1.93548, 1.89873, 1.86335, 1.85169; density hypothesis 2).
Theorem D at σ = 0.86, 0.87, B−: 1.73490, 1.69226, 1.67832 (routine comparison 3/(2σ): 1.74419, 1.72414, 1.71746).

## What the half-plane does and does not do here
- Does: Lemma 1 (L(1 − B − η + iv, ψ) ≪ (k(|v|+2))^{B+η−1/2+ε}, from conductor-uniform Lindelöf beyond B and the
  functional equation) makes every off-diagonal sum in Halász–Montgomery ≪ N^{1−B} X^{B−1/2+ε}: Theorem A.
- Does not: the class-II cut. Under the half-plane and Petrow–Young the conductor-aspect exponent of an L-sum of length
  N at X = N^τ is min(τ/2, 1 − B + (B − 1/2)τ, 1/2 + τ/6, B) = 1/2 + τ/6 for 2 ≤ τ ≤ 6(B − 1/2); interpolated
  μ-lines and interpolated moments give nothing better (PROOF §5, verify B1–B4). So the polynomials of length
  X^{1/(6σ−3)} still need squaring, which costs exactly 4/(6σ − 3).
- Earlier tree rounds R299/R304 said the half-plane cannot help Halász–Montgomery or density below it; DensityB corrected
  that for zeta; here the large-value part is corrected for families too (Theorem A), while the family analogue of
  R304's density claim holds for the uniform exponents (class-II barrier, Proposition B).

## verify.py output (uv run --no-project python3 verify.py; 19 s, max RSS 29 MB, one process)
```
B = 0.8733819666132075  s* = 0.8560920118872971
PASS A1 gB(tauB(s)) = 2s - 1 (HM threshold), exact at 3 points of an affine identity
PASS A2 tauB(3/4) = 1 and tauB((1+B)/2) = 2
PASS A3 tauB(s*) = 3s* - 1 (crossing with Jutila k=2 family range)
PASS A4 slope of tauB = 2/(B-1/2) > 3, so tauB > 3s-1 exactly for s > s*
PASS A5 tauB > 4s-2 (Huxley family range) on [0.84, B]
PASS A6 tauB(sPY) = 12 sPY - 9 (Petrow-Young line crossing)
PASS A7 slope tauB < 12, so tauB > 12s-9 exactly for s < sPY
   sPY = 0.9005218234570964  (Theorem A beats every known family range on (s*, sPY))
PASS A8 hull: gB(t) <= 1/2 + t/6 iff t <= (B-1/2)/(B-2/3) = 1.8063; tauB(sPY) <= that
PASS A9 hull: gB(t) <= t/2 for t >= 1 (gB(1) = 1/2, slope B-1/2 < 1/2)
PASS A10 Petrow-Young point (1/2,1/6) lies below the (H) chord from (1-B,B-1/2) to (B,0)
PASS B1 on [2, 6(B-1/2)] the hull minimum is the Weyl vertex 1/2 + t/6 (endpoints; pieces affine)
PASS B2 cut 6s-3 <= 6(B-1/2) for s <= B, so exclusion holds exactly for tau < 6s-3
PASS B3 (s-u)/mu*(u) is maximal at u = 1/2 (decreasing in u for s < B; 50 grid points)
PASS B4 at u = 1/2 that cut equals 6s-3
PASS C1 unconditional: sup_{[6s-3,12s-6]} F_Jutila <= (4-4s)/(6s-3) on 22 points of (5/6, 7/8]
PASS C1a 9s-9/2 >= 2(3s-1) iff s >= 5/6
PASS C1b 2 + (8-16s)/(9s-9/2) = 2/9 identically
PASS C1c 2/9 <= (4-4s)/(6s-3) iff s <= 7/8
PASS C1d (19.5-23s)/(9s-4.5) <= 2/9 for s >= 0.82
PASS C1e Jutila k=2 Montgomery range is 3s-1 (14s-10 >= 3s-1 for s >= 9/11)
PASS C2a tauB <= 8s-5 on [0.84, B] (affine)
PASS C2b 3s-1 <= 8s-5 for s >= 4/5
PASS C2c s=0.840000: min_a G with tauM = tauB equals (4-4s)/(6s-3) = 0.313725 (got 0.313725)
PASS C2c s=0.850000: min_a G with tauM = tauB equals (4-4s)/(6s-3) = 0.285714 (got 0.285714)
PASS C2c s=0.860000: min_a G with tauM = tauB equals (4-4s)/(6s-3) = 0.259259 (got 0.259259)
PASS C2c s=0.870000: min_a G with tauM = tauB equals (4-4s)/(6s-3) = 0.234234 (got 0.234234)
PASS C2c s=0.873381: min_a G with tauM = tauB equals (4-4s)/(6s-3) = 0.226077 (got 0.226077)
PASS C3a tauB >= (11/2)s - 13/4 on [0.84, B] (affine)
PASS C3b tauB < (9/2)s - 9/4 at 0.84 and B (no-subdivision variant is weaker: 3/tauB)
PASS C3c Huxley+PY gives only 3/(3s-1): (4s-2) < (11/2)s - 13/4 for s > 5/6

 sigma     lam_a    lam_b    T from Q^a to Q^b      old 3/(2s)  new 3/tau0   uniform 4/(6s-3)  DH
  0.860000  0.6279   0.9368   Q^0.688 .. Q^13.83     1.74419     1.73490      1.85185        2
  0.865000  0.6329   0.9380   Q^0.724 .. Q^14.14     1.73410     1.71332      1.82648        2
  0.870000  0.6379   0.9392   Q^0.762 .. Q^14.45     1.72414     1.69226      1.80180        2
  0.872500  0.6404   0.9398   Q^0.781 .. Q^14.61     1.71920     1.68193      1.78971        2
  0.873382  0.6413   0.9400   Q^0.788 .. Q^14.66     1.71746     1.67832      1.78548        2
PASS D1 lam_a < lam_b on {0.86, 0.87, B-} (nonempty regime)
PASS D2 s=0.86000 lam=0.6279: new <= (3-3s)/tau0 = 0.242886 < old(subdiv) = 0.244186 <= old(Jutila form) = 0.244186
PASS D2 s=0.86000 lam=0.7824: new <= (3-3s)/tau0 = 0.242886 < old(subdiv) = 0.244186 <= old(Jutila form) = 0.254545
PASS D2 s=0.86000 lam=0.9368: new <= (3-3s)/tau0 = 0.242886 < old(subdiv) = 0.244186 <= old(Jutila form) = 0.254545
PASS D2 s=0.87000 lam=0.6379: new <= (3-3s)/tau0 = 0.219994 < old(subdiv) = 0.224138 <= old(Jutila form) = 0.224138
PASS D2 s=0.87000 lam=0.7886: new <= (3-3s)/tau0 = 0.219994 < old(subdiv) = 0.224138 <= old(Jutila form) = 0.232836
PASS D2 s=0.87000 lam=0.9392: new <= (3-3s)/tau0 = 0.219994 < old(subdiv) = 0.224138 <= old(Jutila form) = 0.232836
PASS D2 s=0.87338 lam=0.6413: new <= (3-3s)/tau0 = 0.212505 < old(subdiv) = 0.217462 <= old(Jutila form) = 0.217462
PASS D2 s=0.87338 lam=0.7906: new <= (3-3s)/tau0 = 0.212505 < old(subdiv) = 0.217462 <= old(Jutila form) = 0.225640
PASS D2 s=0.87338 lam=0.9400: new <= (3-3s)/tau0 = 0.212505 < old(subdiv) = 0.217462 <= old(Jutila form) = 0.225640
PASS D3 at s*: tau0 = 2s* (no gain at the left end)

 sigma     DH   Huxley+PY 3/(3s-1)   uniform 4/(6s-3)   tauB      3s-1   4s-2
  0.840000  2.0   1.97368              1.96078            1.4821   1.5200  1.3600
  0.850000  2.0   1.93548              1.90476            1.5356   1.5500  1.4000
  0.860000  2.0   1.89873              1.85185            1.5892   1.5800  1.4400
  0.870000  2.0   1.86335              1.80180            1.6428   1.6100  1.4800
  0.873382  2.0   1.85169              1.78548            1.6609   1.6201  1.4935
PASS E1 4/(6s-3) = 2 at s = 5/6 and is decreasing, so < 2 exactly for s > 5/6
PASS E2 4/(6s-3) < 3/(3s-1) iff 12s-4 < 18s-9 iff s > 5/6 (equality at 5/6)
PASS F1 3/(2s) <= 4/(6s-3) iff 9 - 10s >= 0: equality at 9/10, holds on (5/6, 7/8]
PASS F2 Ivic 3831/4791 and Bourgain 3734/4694 lie below 5/6

runtime 18.67 s
ALL PASS
```

## Novelty pass
- arXiv math.NT listings 2026-11 (140 entries) and 2026-10 (854), titles scanned for zero density, Dirichlet L,
  L-function, large sieve, half-plane, zero-free, large values, Dirichlet polynomial, subconvex, Weyl, arithmetic
  progression, Linnik, Bombieri (39 hits). Read: 2610.12234, 2610.15291, 2610.26835, 2611.00924 (half-plane constants),
  2610.18491 (quasi-RH applications to Chebotarev and class groups), 2610.15989 (zeta, short strips), 2610.18923
  (large values, trilinear sums), 2610.37371 (sym² weight aspect), 2610.04045 (large sieve, fixed order). None gives a
  family zero-density bound below 2 on [0.84, B).
- arXiv search (13 queries: family zero density, Guth–Maynard with characters, Jutila, Huxley 12/5, Heath-Brown,
  Petrow–Young and zero density, hybrid subconvexity, density hypothesis for families): newest family paper is
  Chen–Gupta–Li 2507.08296v2 (A = 7/3 uniformly; aimed at σ ∈ [0.7, 0.8]; they quote Huxley's 3/(3σ − 1) for
  σ ≥ 0.8, which needs a hybrid Weyl input to go below 2 and in any case is beaten by Theorem C on (5/6, 7/8]);
  Pintz 1804.05552 (σ near 1, log-free); 2211.00260 (sparse moduli); 2507.15334 (uses A = 7/3). No paper applies
  Petrow–Young to family zero density.
- Jutila 1977 read in full text of the theorem page (open scan, matwbn.icm.edu.pl aa3216.pdf): its family exponents
  are max(2, ...); p. 57 lists only density-hypothesis results for families (σ ≥ 5/6 earlier, 21/26 by him).
- Not read: Heath-Brown, J. London Math. Soc. 19 (1979) 221–232 (paywalled), which also treats Dirichlet
  L-functions. In this scheme his hybrid exponent 3/16 gives the cut (16/3)(σ − 1/2) < 2 for σ < 7/8, so a bound
  below 2 on (5/6, 7/8) there would need an L-sum estimate of twelfth-moment strength; reviewer please check.
- Zenodo API ("zero density" Dirichlet; "zero-density" "L-functions"; Weyl bound zero density characters; 25 newest
  each): nothing relevant. GXL evidence API: daily quota exhausted (HTTP 429), not used.
- ANTEDB (main, 2026-11-05): no Dirichlet L-function chapter. Our files: COROLLARIES (R354: conductor-aspect Lindelöf
  and the Petrow–Young μ-line, used here), DensityB (zeta), NOGO R304; nothing on family density below 2.

## Self-review (two passes; details in PROOF §7)
- Hidden hypotheses: Jutila's theorem is used as stated (k = 2, T ≥ 2, both families, no subdivision; all lengths
  used have N ≤ X); Petrow–Young is used for primitive characters at heights ≤ T + log²T; the zero detection is the
  standard one with a short mollifier; Theorem D needs t-subdivision (lower end λ_a) and class II absent (upper end λ_b).
- Newer records: none found (above). Gain only from an unproved input? No: Theorem C is unconditional; the half-plane
  results (Theorems A, D) are conditional on the Lean-checked half-plane only.
- Weak points: Theorem C is a routine combination (2023 input + 1977 machinery); its novelty rests on the searches
  above and on Heath-Brown 1979 not containing it.
- Review (2026-11-05), both gaps repaired: (1) principal characters are ζ, not O(log²) zeros: Lemma Z now runs the
  detection for non-principal primitive χ only (no residue) and counts ζ's zeros by Ivić's 3/(2σ) ≤ 4/(6σ − 3) on
  (5/6, 9/10] (verify F1, F2), or by DensityB's 3/τ0 in Theorem D; (2) the old record is now 3/(3σ − 1) as stated by
  Chen–Gupta–Li, which Theorem C beats exactly for σ > 5/6 (E2).

## Leads
- Break the class-II barrier: a hybrid twelfth moment Σ_χ ∫ |L(1/2 + it, χ)|^{12} dt ≪ (qT)^{2+ε} would turn the cap
  4/(6σ − 3) into 3/(τ_M + 1 − σ) (3/(2σ) unconditionally, the zeta-type 3/τ0 with the half-plane). It is known in
  the q-aspect for smooth or prime-power moduli (Nunes 1802.05180, Milićević–White 1908.04833): check whether their
  proofs give the t-hybrid form; if so, a single-modulus result for those moduli follows at once.
- Below σ = 5/6 the cut 6σ − 3 is < 2 and this scheme cannot beat 2; there the remaining lever is again the L-sums
  (hybrid moments), not large values.

RH relevance: N7/N11 (zero density in families of Dirichlet L-functions).

Runtime of this agent: about 1 h (20:15–21:15 UTC). Downloads deleted.
