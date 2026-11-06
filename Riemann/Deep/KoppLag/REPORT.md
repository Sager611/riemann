RESULT: Kopp–Lagarias Problem 1.7 solved: the unit-generated real quadratic orders with Cl(O) = Cl(O)[2] are exactly the 89 orders of their Tables 2A/2B (largest Δ = 68640), and exactly 84 of them have one class per genus; elementary and unconditional for Δ < 10¹⁸, complete for all Δ given the Lean-checked zero-free half-plane (via our explicit Littlewood bound).

# Kopp–Lagarias Problem 1.7 without GRH or Siegel (node N11)

Agent kopplag, 2026-11-06, 12:21–13:08 UTC. Brief: `.scratch/orca/kopplag-brief.md`. Proofs: [DERIVATION.md](DERIVATION.md). Novelty: [NOVELTY.md](NOVELTY.md).

## Claim
Let O_Δ be the order of discriminant Δ = n² − 4 (n ≥ 3) or Δ = n² + 4 (n ≥ 1): these are all unit-generated real quadratic orders (Kopp–Lagarias Prop. 2.3).
1. Cl(O_Δ) = Cl(O_Δ)[2] exactly for the 53 discriminants n² − 4 and the 37 discriminants n² + 4 of Kopp–Lagarias Tables 2A/2B
   (89 distinct orders, Δ = 5 in both families; 60 maximal, 29 non-maximal; largest Δ = 68640 = 262² − 4, f = 2, Cl ≅ (ℤ/2)⁴, Cl⁺ ≅ (ℤ/2)⁵).
   Every Cl and Cl⁺ structure and every conductor in their tables is confirmed.
2. One class per genus (Cl⁺ = Cl⁺[2]) holds for exactly 84 of them: all except Δ = 221, 725, 1517 (Cl⁺ ≅ ℤ/4) and 3965, 7565 (Cl⁺ ≅ ℤ/2 × ℤ/4).
3. Inputs: for Δ < 10¹⁸ only elementary algebra and an exact computation (no GRH, no zero-free region); for Δ ≥ 10¹⁸ the bound |L(1, χ)| ≥ 0.11/max(3, log log q)
   (Littlewood78, from Bebel's Lean-checked half-plane); any constant c ≥ 0.037 in place of 0.11 suffices (the induction base is then k = 16).
Side result: Kopp–Lagarias's text (l.1609) says 86 orders; their tables list 89 distinct real orders (count erratum).

## Before
Kopp–Lagarias (arXiv 2512.11311 v3, Apr 2026) prove finiteness ineffectively (Hua's Siegel-type theorem), list the orders with Δ < 10¹⁰ as
"conjecturally complete", pose Problem 1.7 and write (l.1743) that the complete classification "seems out of reach at this time, in parallel with the
idoneal number problem being unsolved". For maximal orders of extended Richaud–Degert type, Louboutin–Mollin–Williams (1993) had the exponent-2 list
with one possible exception (none under GRH). For non-maximal orders only class number one was classified (KL Thm 1.4); the 2-torsion case was known only
by computation (Δ < 10¹⁰).

## Method
1. **Lemma A (elementary).** For Δ = n² ∓ 4 ≥ 100, a primitive α ∈ O_Δ with |N(α)| < √Δ/2 is a unit. Proof: normalize α by powers of ε_n = (n + √Δ)/2,
   write α = (x + y√Δ)/2; the size bounds force x − ny ∈ {−1, 0, 1}, parity forces x = ny, so α = yε_n.
2. **Corollary B (split-prime sieve).** If Cl(O_Δ) is 2-torsion, no prime p with p² < √Δ/2 splits in O_Δ: else 𝔭² would be a primitive principal ideal of norm p².
3. **Effective bound.** Class number formula for orders h(O)R(O) = (√Δ/2)L(1, χ_{Δ₀})∏_{p|f}(1 − χ(p)/p), R(O) ≤ log ε_n < ½ log(Δ + 4),
   φ(f)/f ≥ 1/P_{ω(Δ)}, and Littlewood78 give h(O) > 0.11√Δ/(log(Δ + 4) max(3, log log Δ) P_{ω(Δ)}); genus theory gives h ≤ 2^{μ(Δ)−1} ≤ 2^{ω(Δ)} for 2-torsion Cl.
   Splitting by k = ω(Δ) (Δ ≥ p₁⋯p_k) gives Δ < D₀ = 20485358601780401 ≈ 2.05·10¹⁶ (k = 14 is the worst; k ≥ 15 impossible by an induction using only Euclid's p_{k+1} ≤ p₁⋯p_k + 1).
4. **Sieve** every n with Δ_n < 10¹⁸ (n ≤ 10⁹) by Corollary B: 82 survivors with n > 1000, all n ≤ 5968.
5. **Exact class groups** (reduced indefinite forms, cycles, composition; integer arithmetic only) for all n ≤ 1000 and the 82 survivors.

## Certificate
`cert.py` (`uv run --no-project --with numpy --with python-flint python3 cert.py`, 112–141 s, max RSS 188 MB (`/usr/bin/time -l`), single thread; two runs, identical output): Part 1 re-derives D₀ with arb balls
(bisection per k, k ≤ 60, induction constant 3.022 > 2.039), Part 2 sieves n ∈ [1001, 10⁹] in both families (58 + 24 survivors), Part 3 computes the
class groups of 2080 orders exactly, Part 4 prints the list and compares it with Kopp–Lagarias Tables 2A/2B: IDENTICAL.
Output `logs/cert_output.txt`; survivor log `logs/survivors.txt` (h between 21 and 192, none 2-torsion).

## Checks
`checks.py` (48 s) and `pari_checks.gp` (gp, 70 s), outputs `logs/checks_output.txt`, `logs/pari_checks_output.txt`:
- C1 Lemma A by brute force over its normalized window: 6061 orders (all n ≤ 3000, 80 random n ≤ 10⁹), no counterexample.
- C2 principal cycles are the predicted 2-cycles for all n ≤ 20000 (Lagrange's form of Lemma A).
- C3 Corollary B directly: 58836 (n, p) pairs, 𝔭∘𝔭 reduced, cycle walked: never principal; sharpness: for n = p², N = −1, p² < √Δ but 𝔭² principal, so √Δ/2 cannot be doubled.
- C4 genus count |Cl⁺[2]| = 2^{μ−1} (Halter-Koch Thm 6.5.2) on all 1998 orders with n ≤ 1000.
- C5 the numpy sieve equals an independent brute-force Jacobi-symbol sieve on n ∈ [1001, 3·10⁶] (both families);
  C9 an independent C sieve (`sievecheck.c`, 3 min) over the full range n ∈ [1001, 10⁹] finds the same 58 + 24 survivors (`logs/sievecheck_output.txt`).
- C6 every listed order with Δ ≥ 100 passes the split-prime test (Corollary B is consistent with the list).
- C8 Lemma A's two endpoint inequalities, arb on Δ ∈ [100, 10⁴⁰]: minima 1.5639 and 3.1440 (> 0).
- PARI (GRH-conditional, cross-check only): quadclassunit agrees with every exact h and every 2-torsion verdict (2080 orders); quadunit = (n + √Δ)/2 for Δ ≠ 5;
  the class number formula for orders holds with |Q − 1| ≤ 1.2·10⁻³⁸ (lfun at 38 digits), which checks the unit index and conductor factors used in step 3;
  L(1, χ_{Δ₀}) ≥ 6.57 × the Littlewood78 bound on all these Δ₀.
- Constants re-checked by hand in DERIVATION §2.2, §3.2–3.4 (ε_n bounds for both signs, the −/+4 cases, μ − 1 ≤ ω, P_k ≥ f/φ(f), monotonicity of g).
- Orchestrator reviewer (6 Nov, DERIVATION §§0–3 inline): CORRECT, SIGNIFICANT; its c = 0.037 base note is in DERIVATION §3.5.
- Two reviewer passes (verdict-only, excerpts inline). Pass 1 (Lemma A, Cor. B, §3): GAP = missing hypothesis n ≥ 1 in Lemma A's statement (it is in §0) and a
  rounding (3.0221 → 3.0220, the step needs only 2.039); both fixed. Pass 2 (§4–5, forms): GAP = the printed definition of a reduced form omitted 2|a| < √Δ + b
  (the code and §4.1 use the full definition); fixed. No mathematical gap found.

## Novelty
arXiv listings math.NT Oct–Nov 2026, arXiv API (Kopp, Lagarias, unit-generated, one class per genus, Richaud–Degert, exponent two, idoneal), Semantic Scholar
citations of 2512.11311 (2, unrelated), Zenodo quoted phrases (3 AI summaries without theorems), HAL (0), GXL Evidence API (only KL and LMW 1993):
no solution of Problem 1.7 and no unconditional exponent-2 list for these orders. The idoneal-number analogue (imaginary) was completed in Oct 2026 from the
same half-plane (OpenAI QRH note, Jaswal on Zenodo); the real unit-generated case was not treated. Details: NOVELTY.md.

## Caveat
- The range Δ ≥ 10¹⁸ rests on Littlewood78 (reviewed CORRECT in this project, not refereed) and through it on the Lean-checked half-plane; Littlewood78 itself
  is not formalized. Robustness: c ≥ 0.037 suffices, a third of its constant.
- F2 (class number formula for orders) and F3 (genus theory for non-maximal orders) are cited (Neukirch I.12.12; Halter-Koch Thm 6.5.2), not reproved;
  both are checked numerically on 1998–2080 orders.
- The exact class-group code is ours (forms.py); PARI agrees, but PARI's class groups are GRH-conditional, so the agreement is a check, not a proof.
- Lemma A is classical in spirit (principal-cycle argument of Richaud–Degert class-number-one criteria); the new part is its use for the whole problem.

**RH relevance: N11** (GRH families, Siegel zeros, class numbers): an open problem whose only obstruction was the Siegel-type ineffectivity of L(1, χ),
closed by the explicit half-plane bound. It is the real-quadratic analogue of Euler's idoneal numbers, which the same input completed in October 2026.

## Files
- `cert.py`: certificate (D₀, sieve, exact class groups, comparison with KL tables); `forms.py`: exact form arithmetic (reduction, cycles, composition, structures).
- `sievecheck.c`: independent C sieve over the full range (cc -O2).
- `checks.py`, `pari_checks_body.gp`, generated `pari_checks.gp`: independent checks C1–C8 and PARI cross-checks.
- `logs/cert_output.txt`, `logs/survivors.txt`, `logs/checks_output.txt`, `logs/pari_checks_output.txt`, `logs/sievecheck_output.txt`.
- `DERIVATION.md` (Lemma A, Corollary B, D₀, algorithms, proof), `NOVELTY.md`.
