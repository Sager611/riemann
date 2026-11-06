RECORD (routine corollaries, not for DISCOVERIES): first explicit Landau–Siegel constant over base fields k ≠ Q, using only the k = Q exclusion (Lean-checked Dirichlet half-plane) and Stark 1974: (1 − β) log d_F ≥ (3 − 2√2)/2 = 0.08578 for every quadratic F/k and every real zero β of L(s, χ_{F/k}), k = Q(i), Q(√−3) or any quadratic field, and ≥ (6 − 4√2)/(2^{n−1} n!) for [k:Q] = n; explicit effective Brauer–Siegel over quadratic k: h_K R_K > 3.975·10⁻³ d_K^{1/2}/((2π)^{[K:Q]/2} log d_K) for every Galois K/k (Stark's residue constant made explicit: c_4 = 43.1621). Cureton–Thorner's Lean-checked 31/32 gives the larger (1/32) log d_F ≥ 0.155 (Q(i)), 0.149 (Q(√−3)) in one line.

# Explicit c(k) over a base field and explicit effective Brauer–Siegel (node N11)

Agent bsexpl, 2026-11-06, 01:25–02:05 UTC; reviewer pass on DERIVATION.md §§2–5: CORRECT. Files here: DERIVATION.md (all proofs), cert.py + cert_output.txt (arb
certificate), adversarial.gp + adversarial_output.txt (PARI checks).

## What this is
Banks (arXiv 2610.31402v2) proves that for every number field k there is an effective c(k) > 0 with (1 − β) log d_F ≥ c(k)
for every quadratic F/k and every real zero β ∈ (0,1) of L(s, χ_{F/k}) = ζ_F/ζ_k, and deduces h_K R_K ≥ c(k) d_K^{1/2}/
((2π)^{[K:Q]} log d_K) for all Galois K/k. He leaves c(k) unevaluated. We give explicit values, by a shorter route.

## Results
| statement | constant | inputs |
|---|---|---|
| Thm A, k quadratic (Q(i), Q(√−3), every quadratic k) | (1 − β) log d_F ≥ (3 − 2√2)/2 = 0.085786 | Lean half-plane for Dirichlet L (V4, C4 cases); Lemma D (D4 case) |
| Thm A, [k:Q] = n ≥ 2 | ≥ (6 − 4√2)/(2^{n−1} n!): n = 3: 0.01430, n = 4: 1.787·10⁻³ | + Stark Thm 3, Aramata–Brauer, Toyama |
| Thm A with Guo's c_Q = 1/100 instead of the half-plane | ≥ min(n/100, (6 − 4√2)/(2^{n−1} n!)); 0.02 for quadratic k | no 2026 half-plane at all |
| Lemma D: ζ_M (M ≠ Q) with a real zero β of multiplicity ≥ 2 | (1 − β) log d_M > 6 − 4√2 = 0.343146 | Stark 1974 p. 140 |
| Stark's Lemma 4 explicit | κ_K > (1 − β_1)/43.1621 for every K ≠ Q | Stark 1974 Lemmas 2–4 |
| Thm B, k quadratic, every Galois K/k | (1 − β) log d_K ≥ 0.085786 for real zeros of ζ_K; h_K R_K > 3.9750·10⁻³ d_K^{1/2}/((2π)^{[K:Q]/2} log d_K) | Thm A + Stark Thm 3, Lemmas 3, 4 |
| Thm B, [k:Q] = n | c_BS = 2c(k)/c_4: n = 3: 6.63·10⁻⁴, n = 4: 8.28·10⁻⁵ | same |
The form (2π)^{[K:Q]/2} is stronger than Banks's (2π)^{[K:Q]} (2^{r_1}(2π)^{r_2} ≤ (2π)^{[K:Q]/2}, w_K ≥ 2).

## Method (details and proofs: DERIVATION.md)
1. Lemma D. Stark's inequality Σ′ 1/(σ − ρ) < 1/(σ − 1) + ½ log d_M (1 < σ < 2) with a double real zero and σ = 1 + (1 + √2)(1 − β)
   gives (1 − β) log d_M > 6 − 4√2 (Stark's "at most one zero" constant is 1/4).
2. Quadratic k. F/Q is V4, C4 or D4. V4/C4: L(s, χ_{F/k}) is a product of two Dirichlet L-functions, so β ≤ T* = B_2 + 2.2·10⁻²⁵⁴
   and (1 − β) log d_F ≥ (1 − T*) log 9 = 0.278. D4: L(s, χ_{F̄/k}) = L(s, χ_{F/k}) identically, so β is a double zero of ζ_M,
   M = F F̄, and d_M ≤ (d_F/d_k)⁴ (conductor–discriminant over k); Lemma D gives (1 − β) log d_F ≥ (6 − 4√2)/4.
3. General k. β is a zero of ζ_M (M = normal closure of F over Q; Aramata–Brauer). Multiple ⇒ Lemma D, log d_M ≤ 2^{n−1} n! log d_F
   (Toyama; r ≤ n conjugates of F; [M:Q] ≤ 2^n n!). Simple ⇒ Stark's Theorem 3 over Q puts β on a quadratic subfield Q(√D) ⊆ F,
   i.e. L(β, χ_D) = 0, so β ≤ T*.
4. Brauer–Siegel. Banks's argument with Lemma D in place of Stark's Lemma 3 and with Stark's c_4 evaluated:
   c_4 = 2c_3, c_3 = exp(21/8 + c_1/2 − (c_2/8)ψ(½)), c_1 = 0, c_2 = 2/log 3 (needs d_K ≥ 3^{[K:Q]/2}, from Minkowski).
   η_k (the real zeros of ζ_k) is bounded the same way through k̃: η_k ≥ min((1 − T*) log 3, (6 − 4√2)/n!); η = ∞ for Q(i), Q(√−3).

## Banks's own chain (not used)
Not explicit in 2610.31402: the lattice-count errors of lem:count, the thresholds of prop:main and thm:precise (X_0(n_K, H)),
"d_F large in terms of k, ℓ, H, C" in thm:quadratic (for k = Q(i): ℓ = 131, H = 33 280, n_K ≤ 520), the absolute c_1 of
lem:bias, and the small-field step (Hermite + mean value theorem + Rademacher), which gives (1 − β) ≳ d_F^{−O(1)} at an
astronomical threshold; Stark's c_2 = 1/c_4 (now 1/43.16) and Mertens for k (Garcia–Lee) are explicit. Table: DERIVATION.md §1.
Evaluating this chain would give a far smaller constant than §§2–4, so we did not.

## Before
- Over k ≠ Q no explicit c(k) was stated. Banks (2610.31402v2, l. 266, 289): effective, not computed; Frei–Kai (Hexagon
  2610.00466, l. 223): ε_K "in principle effective, but we have not attempted to make it explicit".
- One line away, but not stated: Cureton–Thorner (Hexagon 2610.00200v2, Cor. Hecke, Lean: AxiomMath/DedekindZetaZeroFree)
  L(s, χ, K) ≠ 0 on Re s > 31/32 for every finite-order Hecke χ ⇒ (1 − β) log d_F ≥ (1/32) log d_F (0.1553 for Q(i), min d_F = 144;
  0.1488 for Q(√−3), min d_F = 117) and h_K R_K > 1.59·10⁻³ d_K^{1/2}/((2π)^{[K:Q]/2} log d_K) for every K ≠ Q. Liu (arXiv
  2610.26835, Re s > 0.874887, audited sound, unrefereed): 0.1251 log d_F (0.622, 0.596). Hu (Hexagon 2610.00216, 7/8, not Lean).
  None of the three TeX sources mentions Siegel zeros or Brauer–Siegel. Banks (l. 128–139) and Frei–Kai (l. 223) say these
  half-planes would settle the question.
- k = Q: Guo 1/100 (2610.30476), 2610.21511 5.811·10⁻⁴, Banks 3·10⁻³ (q ≥ 10¹⁰); the Lean half-plane itself gives (1 − T*) log 3 = 0.139.
- Our constants are smaller than the Cureton–Thorner corollary; they are the first that use no general-field half-plane, and the
  D4 case uses no 2026 input at all.

## Finding about 2610.31402
Banks's cor:siegel and thm:bs follow from the k = Q theorem plus Stark 1974, with constants depending only on [k:Q]: pass to the
normal closure of F (or of k) over Q, whose degree is ≤ 2^n n!, and apply Stark's Theorem 3 there. His "Stark's reduction stops
at k" (l. 197) holds for the tower K/k, not for the quadratic step F/k; his uniformity obstruction "we do not address it"
(§uniform, l. 965–982: real zeros of ζ_k for non-almost-normal k) is handled the same way through k̃, for k of bounded degree.
Not a flaw in his proofs; it makes the base-field theorem a corollary of the k = Q theorem.

## Certificate
- cert.py (python-flint arb, 200 bits): 6 − 4√2, (3 − 2√2)/2, 1 − T*, c_3, c_4, the Minkowski step for c_2, the Lemma D
  identity, every table entry, the comparison constants. Output cert_output.txt: ALL PASS.
  Run: OMP_NUM_THREADS=1 uv run --no-project --with python-flint python3 cert.py > cert_output.txt (< 1 s, < 50 MB).
- adversarial.gp (PARI/GP 2.19, one thread): perl -e 'alarm 1500; exec @ARGV' gp -q adversarial.gp > adversarial_output.txt
  (32 s, 16 MB).

## Adversarial checks (adversarial_output.txt)
1. All quartic F ⊃ Q(i) (173 fields) and F ⊃ Q(√−3) (272 fields) with d_F ≤ 20 000 (nflist V4, C4, D4, r_2 = 2):
   L(σ, χ_{F/k}) = ζ_F/ζ_k sampled on [½, 1] (29 points): no sign change, minimum 0.240 (Q(i)) and 0.132 (Q(√−3)); min d_F = 144
   and 117 as used above. No real zero exists near the excluded zone [1 − 0.0858/log d_F, 1); nothing contradicts Theorem A.
2. D4 closures (8 fields): d_M ≤ (d_F/d_k)⁴ in all cases, with equality in 4 (the factor 4 is sharp); ζ_M(2) = ζ_k(2)
   L(2, χ_F)² L(2, χ_{F_3}) to 10⁻³⁸, i.e. the double zero is structural.
3. Stark's p. 140 inequality with the computed zeros (|γ| ≤ 20) of one octic M: holds at σ = 1.05, 1.2, 1.6.
4. Theorem B and Lemma 4 against 10 Galois K/k (degrees 4 and 8; bnfinit, GRH-conditional, sanity only): κ_K and h_K R_K
   exceed the bounds by factors ≥ 300.
5. k = Q specialization: Theorem A with M = F returns min(c_Q, 6 − 4√2) = c_Q (Guo, 2610.21511, or the half-plane's 0.139);
   it is never better than its input, so nothing beats Kerr/Guo except through the Lean half-plane itself (cert §7).
6. Limit of the positivity step (D4): a nonnegative combination of the characters of D4 with a simple pole contains ρ at most
   twice (classes r² and r, s, rs), so ζ_M has the largest possible multiplicity; in the family 1 + χ_k + u(χ_b + χ_c) + (1 + u)ρ
   the bound is 2(√m − 1)²/(3m − 2), m = 1 + u, largest at u = 1 (ζ_M): 0.0858. A smoothed explicit formula does not close the
   gap to the Cureton–Thorner corollary: in the asymptotic LP 2Φ(a + b) − Φ(a) > f(0)/2 (Φ the Laplace transform of a positive-
   definite weight f), the triangle weight gives b ≈ 0.474 and the Bohman weight b ≈ 0.513 instead of 6 − 4√2 = 0.343, i.e.
   c(k) ≈ 0.128 < 0.155 (uncertified estimate, lp_estimate.py, about 10 min).

## Novelty check (2026-11-06)
- python3 .scratch/lst.py 'Brauer.?Siegel|Siegel zero|Landau.?Siegel|exceptional zero' math.NT --zenodo '"Brauer-Siegel"':
  2610.15483, 2610.21099, 2610.30476, 2610.31402 only (all known); Zenodo: 3 unrelated records.
- Wider title scan (Dedekind, class number, Hecke, real zero, zero-free, half-plane, residue, Stark, exceptional, Chebotarev;
  math.NT 2026-10/11, 45 hits): no explicit c(k); 2610.14257 is GRH-conditional.
- Related work of v2: Frei–Kai TeX (not explicit), Cureton–Thorner TeX (half-plane only), Liu 2610.26835 TeX and Hu TeX
  (no Siegel or Brauer–Siegel statement). Our repo: SiegelBanks audit (noted the Q(√−3) half-plane case), R370 harvest.
- Verdict: the explicit statements are not in print, but each is a short deduction from published or Lean-checked results
  (Stark 1974; the half-planes). Routine: COROLLARIES.md, not DISCOVERIES.md.

## Caveats
- Thm A's V4/C4 cases and the general simple-zero case use the Lean-checked Dirichlet half-plane (Bebel, B_2 + 2.2·10⁻²⁵⁴); with
  OpenAI's 7/8 (Lean port, Riemann/Deep/QRH) nothing changes ((1/8) log 9 > 0.0858); with Guo's audited 1/100 the quadratic-k
  constant drops to 0.02.
- Stark 1974 (Lemmas 2–4, Theorem 3) is used as published (refereed; read from the GDZ scan, not kept); we rechecked the
  arithmetic behind c_3 and the case split for zeros outside S. Aramata–Brauer and Toyama are classical (general-n case only).
- Not formalized in Lean. Constants are not optimized beyond the plain-positivity limit (item 6).

RH relevance: N11.
