# PROBES — falsification probes on Masullo, *The Music of the Primes*

Code: `code/`; logs: `logs/`. All runs: `bash .scratch/memcap.sh 290 <s> uv run --no-project --with cypari2 python3 …` (PARI 2.17.2), single process, peak memory in each log (≤ 30 MB).

## P1–P5. Theorem 12.5.1.1 with a wide strip = Weil's explicit formula in Connes' form (constants check)
`code/p1_explicit.py`, `logs/p1_explicit.log`. Trivial character, φ(g) = f(log|g|), f a C^∞ bump on (0.3, 2.3) (so φ(C_F^1) = 0, support |g| ∈ (1.35, 9.97)). Left side of (12.5.1.1)–(12.5.1.2) with the manuscript's constants ½ (real), 1/2π (complex), log q_v/μ_v(O_v^×) (finite); right side (12.5.1.3) with vol(C_F^1) = 1, dt/t on R_{>0}, all zeros with 0 < γ < 260 (zeros taken from the factors ζ·L(χ_D), resp. ζ·(ζ_F/ζ); PARI's `lfunzeros` on the Dedekind zeta itself misses close pairs: 278 vs 284 for Q(i), 264 vs 272 for Q(√−3)).

| field | zeros < 260 | LHS | RHS | LHS − RHS |
|---|---|---|---|---|
| Q | 114 | 1.112861987530 | 1.112861978532 | 9.0e-09 |
| Q(i) | 284 | 1.086628444864 | 1.086628429411 | 1.6e-08 |
| Q(√−3) | 272 | 1.102987814575 | 1.102987799629 | 1.5e-08 |
| Q(x³−x−1) | 469 | 1.057709678770 | 1.057709654536 | 2.4e-08 |
| Q(√5) | 293 | 1.036022767994 | 1.036022751118 | 1.7e-08 |

Residuals = truncation of the zero sum. Control: factor 1 instead of ½ at the real place breaks Q by 0.27. **Verdict: the statement's constants and Mellin normalisation are right; with I ⊇ critical strip, Theorem 12.5.1.1 is exactly Connes' form of Weil's explicit formula (KNOWN).** Hence for narrow I the theorem is equivalent to GRH; all content lies in its proof.

## P6. Artin-forced multiplicities (Galois fields) — `code/p2_mult.py`, `code/p6b_a4.py`, logs `p2_mult.log`, `p6b_a4.log`
- S₃: K₆ = splitting field of x³−x−1 (x⁶−3x⁵+5x⁴−5x³+5x²−3x+1). ζ_{K₆} = ζ·L(χ₋₂₃)·L(ρ₂)². At ½ + 5.115683…i (first zero of L(ρ₂) = ζ_{F₃}/ζ): |ζ_{K₆}| = 1.1e−67, |ζ′_{K₆}| = 2.5e−33, |ζ″_{K₆}| = 28.1 → double zero, 2 ≤ [K₆:Q] = 6.
- A₄: F₄ = Q[x]/(x⁴+8x+12), disc 5184, L(ρ₃) = ζ_{F₄}/ζ (conductor 5184; lfuncheckfeq −98 bits). Zeros at ½ + i·{1.81488, 2.70244, 3.72384, …} (11 below height 10); |L(ρ₃, ½+1.81488i)| = 4.6e−30, |L′| = 3.78 (simple in L(ρ₃)). In K₁₂ = Galois closure (degree 12, disc 2¹⁸3¹⁶), ζ_{K₁₂} = ζ_{k₃}·L(ρ₃)³: each is a zero of order ≥ 3 ≤ 12.
- **Verdict:** Artin-forced multiplicities d ≤ √[K:Q] < [K:Q] never violate the claimed bound N ≤ [F:Q]: no contradiction. (Used in P12.)

## P7. F = Q consequence of Thm 13.3.1.1 + Weil: all non-central zeros of ζ and of Dirichlet L simple — `p2_mult.log`
Zeros 0 < t < 100 of L(s,(D/·)), D = 1, −4, −3, 5, −7, 8: 29, 50, 46, 54, 59, 61 zeros; min |L′(ρ)| = 0.79, 1.31, 1.17, 0.99, 0.84, 1.67. Consistent (no known multiple zero). **The manuscript claims more than GRH** (p. 528 l. 31915–31925: N(ρ) ≤ [F:Q] off the central point; for Q, the Simple Zeros Conjecture for ζ and all Dirichlet L-functions); not refutable by data; a red flag of strength, not a contradiction.

## P8. CM central zeros beyond [F:Q] — `p2_mult.log`
E₁₂₅₄: y² = x³ − 1254²x (CM by Z[i]); root number −1; ellrank [3, 3] with points (−912, 25992), (−528, 26136), (1650, 43560), regulator 22.69; analytic rank 3. Gross–Zagier–Kolyvagin: analytic rank odd and ≠ 1 ⇒ ord_{s=1} L(E,s) ≥ 3 unconditionally. L(E,s) = L(s−½, ψ), ψ the unitary Hecke character of Q(i), trivial on τ(R_{>0}) (x₀ = 0): a zero of multiplicity ≥ 3 > [Q(i):Q] = 2 exactly at the central point s = ½ − ix₀. (13.3.1.2) forces M(φ) to vanish there to order 2, so both the claimed formula and Weil's lose this term: **no contradiction**; the manuscript itself concedes this (Prop. 14.1.2.5, p. 530–532, conditional on BSD-type (14.1.2.4)); P8 is an unconditional instance.

## P9. Definition of the central point (p. 528–529) — `p2_mult.log`
Text: c₀ is "the unique point of each connected component allowed to occur as a pole of multiplicity higher than 1 of the partial Euler products Π_{v∈S} L_v(χ_v|·|^s)". F = Q(i), χ = 1, S ⊇ places above 5: Π_{v|5} L_v = (1 − 5^{−s})^{−2} has double poles at s = 2πik/log 5 for every k ∈ Z (at s₁ = 3.9040i: |g| = 0, |g′| = 1.6e−38, |g″| = 5.18 for g = (1−5^{−s})²). **Literal statement false** (infinitely many non-central double poles); intended reading presumably "a pole of L_v for infinitely many v" (Thm 1.7.1.1(2)–(3)). Minor.

## P10. Local multiplicity m_p (p. 533) — `p2_mult.log`
ψ = Hecke character of Q(i) of y² = x³ − x: for the 4783 split primes 5 ≤ p < 10⁵ the two poles of L_v L_{v̄} never coincide (a_p² ≠ 4p) → m_p(ψ) = 1, m_F(ψ) = 1 < [F:Q] = 2; χ = 1: m₅ = 2 = [F:Q] (P9). Matches the text's heuristic (μ(U) = 0 unless Galois type). Consistent.

## P11. The §1.4.5 approximation (brief item (c)) — `code/p11_approx.py`, `logs/p11_approx.log`
φ_z(t) = t^{−z}1_{t>1} gives M(φ_z)(ρ−½) = 1/(z − (ρ−½)). Over zeros of ζ (z = 3/2): Σ_{0<γ<T} 2/|z − iγ| = 1.18, 2.34, 4.06 for T = 100, 300, 1000 (diverges like (log T)²/2π); symmetric sum Σ 2Re 1/(z − iγ) = 0.0511, 0.0613, 0.0662 (converges). So φ_z is not reached in any topology where the spectral side is absolutely continuous: the limit ε → 0, h → 1 is a conditionally convergent limit over zeros, and the strip-restricted sub-sum converges only via zero-density estimates for the off-strip part (Ingham-type N(σ,T) ≪ T^{1−c(σ−½)}). §1.4.3 claims it by "a similar reasoning" with Montel's theorem (p. 29), which needs local uniform bounds that are not given. **GAP, standard-fixable** (classical explicit-formula estimates); Ch. 13, where it was to be done, is empty (fn 3 p. 522). Not where the difficulty lies: subtracting Weil's formula from (12.5.1.1) for compactly supported φ already gives GRH without any approximation (Laplace-transform uniqueness).

## P12. Sobolev-degree test of the gluing step (★) — logic on P1, P6 (see CRUX §5)
I = {0}, s = 4 (all cited hypotheses allow it), F = K₁₂ (A₄ closure), χ = 1: spectral side caps multiplicity at min{N, ⌊(s+1)/2⌋} = 2 (p. 440, (10.2.1.6); same as Connes' n < (1+δ)/2 for even s), the zeros of L(ρ₃) have multiplicity ≥ 3 in ζ_{K₁₂}, the p∉S strata traces vanish for every cap (one pole progression per p), Thm 12.4.1.1 has no s. Chain as written ⇒ Σ_{L(ρ,ρ₃)=0} (N−2) M̂(ρ) = 0 for all φ: false. **(★) is not a formal identity: it fails at (α, s) = (0, 4) on K₁₂; any proof must be s- and α-sensitive. Decisive for the structure, not for a specific written line (the line is missing).**
- Side remark: for odd s, (10.2.1.6) gives ⌊(s+1)/2⌋, one more than Connes' count (largest n < (1+s)/2) and than the H^{s/2} Sobolev count m < (s−1)/2 if W^{s,2} is the (1+log²|g|)^{s/2}-weighted space; to check against Ch. 9 (not load-bearing at s = 2[F:Q]).

## P13. Thm 12.5.1.1 literally at α = 0 — `p1_explicit.log`
For α = 0 the interior of the strip is empty, so (12.5.1.1) would read local terms = M(φ)(|·|^{1/2}) + M(φ)(|·|^{−1/2}). For Q: 1.112862 vs 0.867415 + 0.236398 = 1.103813; difference 0.009048 (zero sum), numerical error 1e−8. **Statement false at α = 0 as written** (the capped boundary terms (10.2.2.5) are omitted; the intended α = 0 version is 13.3.1.1). Minor statement error.

## P14. Thm 1.7.1.1 (1), (3) (polar multiplicities) — by hand
(1) Poles of L_{v_i} = (1 − χ(π_{v_i}) q_i^{−s})^{−1} (unitary χ) lie on i(θ_i + 2πZ)/log q_i. Two common points ⇒ log q₁/log q₂ ∈ Q ⇒ p₁^a = p₂^b: impossible for distinct primes. True. (3) χ of finite order n: s = it is a pole of L_v iff χ(π_v) = q_v^{it}; for t ≠ 0 this holds for infinitely many v only if t log q_v ∈ (2π/n)Z for infinitely many primes, again impossible; s = 0 occurs for the infinitely many v with χ(π_v) = 1 (Chebotarev). True. KNOWN-level, not load-bearing.


## P12 (refined after review). Sobolev test, internal form — `code/p6c_certify.py`, `logs/p6c_certify.log`
Hardy Z of L(ρ₃) (root number +1, feq to 2⁻⁹⁸): +0.4030 at t = 1.80, −0.4060 at t = 1.83 → a zero on Re s = ½ (certified by sign change); in ζ_{K₁₂} = ζ_{k₃}L(ρ₃)³ it has multiplicity ≥ 3. At α = 0 the spectral side (10.3.2.3) (p. 454) counts it min{N,2} = 2 times at s = 4 and 3 times at s = 6 (both admissible); the asserted geometric side (p. 520: Thm 11.3.2.1 + Thm 12.4.1.1 + pole terms) is the same for both. **Reviewer (independent, 2nd round): D1 CORRECT** — the parameter-uniform gluing (★)+(★★) is false for some admissible s; D3: verdict GAP correct; D2 corrected to "implies GRH" (not "equivalent"). Same contradiction via Weil at s = 4 alone (no on-line zero needed). S₃ analogue (K₆, double zeros; Hardy Z of L(ρ₂): +0.133 at 5.10, −0.121 at 5.13) needs s = 2 vs 4, excluded by Def. 10.4.3.7 (s ≥ 3): hence A₄.


## P15. Hidden GRH input through 1/L? (scout + own spot check, music.txt l. 16464–17119, 22417–22746, 38792–39547)
Strong Co-primality (Thm 7.5.0.1, p. 374–378): γ_c^{(k)}(ϕ) = ∂_s^k(Z_S(F_S(ϕ₀), χ|·|^s_{F,S}) L_S(χ, s))|_{s=1−s_c} (l. 22611–22613), s_c a pole of a local factor (Re s_c = 0 for unitary χ), so L_S is evaluated on Re s = 1; non-vanishing "a consequence of the non-vanishing of the Hecke L-function L(χ, s) on the 1-line" (l. 22665). B.7: reciprocal partial Euler products only on Re s = 1; the one contour shift (for −ζ′_F/ζ_F·x^w/w) stays "contained in the classical zero-free region" (l. 39074–39080). §6.1: only local L_v. **No contour crossing unknown zeros; no hidden GRH input found in these ranges** (negative result). Ch. 11.2 not checked.

