AUDIT PASSED: no false or unjustified step in the zero-exclusion core of 029 Theorem 1.2; every internal step (§§3–8, App. A) re-derived, all exponents exact, local identities brute-forced (planted errors fail); KP84 + erratum and Goldmakher–Louvel supply exactly what is used (checked at statement level; their refereed proofs not re-proved).

Source: openai/math preprints/Primitive-roots-for-every-admissible-integer-base-October-4-2026/build/sections/*.tex (fetched 14 Nov 2026). Chain with line numbers: [FULLPROOF.md](FULLPROOF.md).

## Short answer
- Theorem 1.2 (01-introduction.tex l.38–45) holds as written: F **cyclotomic** with μ12 ⊂ F, every finite-order η: no zeros of L_F(s,η) in Re s > 1 − 10⁻⁶. The proof gives 1 − 1/20000 (06-zero-free l.716–717, 755–757).
- Mechanism: one probe I_Z = I(Z^{9/10}, Z^{133/1000}, Z) evaluated twice. Reflection: Cauchy–Schwarz + reflected second moment (P4.1) ⇒ I_Z ≪ Z^{933/2000+ε} = Z^{7/15 − 1/6000 + ε}. Poisson: I_Z = c0 f_η(Z) + O(Z^{7/15 − .0004}), f_η the inverse Mellin transform of e^{(s−5/6)²}H_η(s)/L^S(s,η), c0 = R_S²M(1/6)Ŵ1(1)/6 ≠ 0. Hence f_η ≪ Z^{7/15 − 1/20000}; Mellin + identity theorem ⇒ 1/L^S holomorphic on Re s > 1 − 1/20000.
- Consequence (as kusasyu36 notes): F = Q(μ12), η = 1 gives ζ(s) ≠ 0 for Re s > 1 − 1/20000; F = Q(μ_{12q}) gives the same strip for every Dirichlet L-function mod q.

## The open items of R470 — status
| Item | Where | Status |
|---|---|---|
| Unfolding/Fubini normalization | L6.1 (04-reflection l.276–303), L6.3 (l.811–848), Poisson orbits (05-poisson l.303–321), triple Mellin (l.353–389) | Resolved: each is a measure-preserving substitution (ℓ_Tℓ_U = ℓ_{TU}, y ↦ y⁻¹, quotient measure with counting measure on U_S⁶, free action on h ≠ 0), dominated by (6.4), (6.18), (6.24) and the sheet bound O(T^{−σ}). A normalization constant multiplies every term of I_Z alike; only c0 ≠ 0 enters. |
| Lemma 6.3, local data above 2 and 3 | 04-reflection l.505–931; 02-arithmetic l.18 | Places above 6 lie in S, so 2 and 3 enter only via H_S, E ⊕ B, C_S = H_S² (03-theta l.233–251) and the S-Whittaker bounds. Checked: alternation (−1 = ζ12⁶), dim V_ℓ = 2s_f + 2r, Lagrangian E, exact reciprocity; local units at p∣d, p∣f, q∣c'; pairing χ_p(q)^{1−4} = χ_p(q)³; no power of q_d in the reflected coefficients (matches a GL2-twist functional equation, dual length cond²/Z). |
| Lemmas 3.2–3.4, 5.4–5.5 | 02-arithmetic l.134–625; 03-theta l.503–775 | Checked line by line; L3.4 identities also numerically. |
| Exponents and both margins | 02a l.69–76; 06-zero-free l.371–683 | Exact rationals ([exponents.log](exponents.log)): reflection margin 1/6000 = (4 − 4a − 3b)/6; good central rows −377/300000; slack after the 10⁻⁴ losses 227/300000 above .0004; without the ξ-move to .998 the good rows sit at +22303/(3·10⁷), so P8.2 is load-bearing and correct (−3.46, U^{.088+ρ}, θ = 0.7505). |
| Global lift of w, Weyl-section compatibility | 04-reflection l.133; App. A l.200–248 | Resolved: in KP's cocycle model the rational splitting is γ ↦ (γ_v, 1)_v (Π_v σ_v = 1 on rational pairs; KP Prop 0.1.3 at tame places), the identity sections used locally; App. A's lift wn(H) = n(H⁻¹)a(−H⁻²)n̄(H⁻¹) re-derived. Any residual discrepancy is a μ3 constant, invisible in Σ|B_m|². |
| Wild-prime Whittaker bounds | 03-theta l.344–396 | Resolved: KP Thm I.2.9(e) (Jacquet module of the exceptional representation = ind(w0ω'), no tameness hypothesis) gives |W(y)| ≪ |y|^δ, δ < 1/3, at every finite place; P4.1 needs only δ > 0. |
| Global reflection | 04-reflection l.333–337, 659–848 | D⁺ = D⁻ by rational automorphy; constant term killed by the angular L⁻¹ (complex inducing character |y|^{1/3}, radial); checked. |
| KP84 + erratum | App. A l.14–167 | Statement level (Numdam scans): Prop 0.1.3 (p.44), Thm I.2.9 (p.72), Cor I.3.6 (p.79: dim Wh = 1 iff r = n−1, 2(c+1) ≡ 0 mod n, or r = n; tameness "not needed"), Thm II.2.1 (p.118), Thm II.2.5 (p.122). Erratum (Publ. IHES 62, p.203): pp.119–120 assertion false; c(η) = lim C_S(η)T(S), T(S) = Π(1+q_v⁻¹) in rank 2; Thms II.2.2–II.2.3 and p.130 modified. 029 normalises each exterior functional at φ°_v and never uses the withdrawn assertion. |
| Goldmakher–Louvel | 04-reflection l.941–1025 | arXiv 1112.1642v2 Def. 1: primitive Hecke characters χ_a of trivial infinite type indexed by squarefree a ∈ I(c), (1) order ∣ n, (2) χ_a(b) = χ_b(a)C([a],[b]) via a homomorphism to a finite group, (3) [a] = [b] ⇒ χ_aχ̄_b primitive mod ab. Thm 1.1: Σ*_{Na≤M}|Σ*_{Nb≤N} λ_b χ_b(a)|² ≪_{k,c,ε} (MN)^ε(M+N)Σ*|λ_b|². 029's family Ξ_c (F(√c)/F, c = q_S, C ≡ 1, cl2) meets (1)–(3); rows d, columns c as used. |

## Falsification (smallest cases; [checks.log](checks.log), [exponents.log](exponents.log))
- P7.3 local sums: C_{t,k}(J) from its definition (05-poisson l.498–507) by brute force over Z/13ⁿ and Z/37ⁿ against the closed values (l.571–578): 136 cases, max relative difference 1.6·10⁻¹². Planted error (τ1 → τ5 in C_{t,1}) fails: 1.95.
- P7.3 Euler factor: series (l.538–548) vs closed form (l.644–654) at 12 random points, p = 13, 37, η built as −γ1γ2ϑ from actual Gauss sums: 4.4·10⁻¹⁶; H(E1 = 0) − 1 formula (l.690–694): 4·10⁻¹⁰. Planted sign ϑ → −ϑ (kills the −D1 = first term of 1/L) fails: 7·10⁻².
- L3.4 finite-field identities at degree-one primes 13, 37, 61, 73, 97: χ(−1) = 1, τ1τ5 = Q, τ0 = −1, (Q^{−1/2}G2)³ = Q^{−1/2}J(ρ,ρ), G1G2 = χ̄(4)J(ρ,ρ)G3: all ≤ 2.4·10⁻¹⁵.
- Footprint: max RSS 16 MB (/usr/bin/time -l); memcap peaks below its 2 s sampling (runs < 10 s).
- Not testable at feasible sizes: the asymptotic second moment P4.1 itself. Earlier passes (R470) stand: kusasyu36's sextic-twist reflection identity for Patterson's theta over Q(ω) (8 primes, 2.5·10⁻¹³) and L6.2's algebra at (ζ12 − 2).

## Reviewer (delicate step pasted inline; two runs)
- Step: L6.2 at a valuation-one row prime (R(t) = ω_pχ_p(t), |ω_p| = 1) and the row–column pairing χ_p(q)χ_q(p)⁻⁴ = χ_p(q)³ with no q_d power.
- Verdict: "GAP: Claim A is correct, but Claim B's unimodular coefficient formula requires (n', d) = 1; a common prime gives a zero masked row." My paste had omitted that support condition; the paper imposes it ((c'n', df) = 1 with the frequency set (6.12), 04-reflection l.467–475; row support (d, hn') = 1, l.591–592) and the reviewer's own computation (R_{p,1,e} = 0 for e > 0) is the paper's reason. Otherwise: "no sign or orientation error occurs".
- Rerun with the support condition stated: "VERDICT: CORRECT" (ω_P = G3/G_{−1}; χ(t)^{−5} = χ(t); unit exponent 2(3k+1) ≡ 2; cross term z·z⁻⁴ = z³; modulus q_{c'}^{−1/2}q_{n'}^{−1}, no q_d power); its extra phase Π_{P∣d}χ_P(n')³ sits in the paper's row factor ρ (n' is frozen there).

## What is not verified (residual risk)
- Proofs of the cited refereed theorems: KP84 (incl. §0.2 rational splitting and §II.1 constant term, used through 029's App. A paraphrase; Prop II.1.2 not located in the OCR scan), the erratum's corrected argument, Goldmakher–Louvel Thm 1.1; classical inputs (Tate's local functional equation, Hecke convexity, Weil-index product formula, Davenport–Hasse, Hilbert reciprocity).
- Reading risk: an AI-written proof read by AI (here and in kusasyu36's audit). The construction is thin by design: an unnoticed exponent loss of 1/6000 on the reflection side, or of 377/300000 ≈ 1/796 on the good rows, would break it (any smaller loss only narrows the width). None found.
- Out of scope: Theorem 1.1's other input, the marked Type II estimate of companion [21] (our PrimRoot REPORT2).

## Remarks for the tree
- TREE N11 says "fields containing μ12"; the TeX says cyclotomic fields containing μ12. The cyclotomic hypothesis seems unused (only μ12 ⊂ F: complex places, χ_p(−1) = 1, −1 ∈ F^{×6}, Q(ω) ⊂ F), agreeing with kusasyu36 §5.18; not claimed by the paper.
- R470's "0.00126 − 0.0004 ≈ 0.00086" ignores the paper's 10⁻⁴ loss budget; exact slack 227/300000 ≈ 0.00076.
- kusasyu36's gaps 3 (global w lift) and 5 (wild Whittaker bounds) close as above; gap 4 (unfolding) closes; gap 1 (KP84) reduces to the quoted statements; gap 2 (GL) to its refereed proof.

## Files
- [FULLPROOF.md](FULLPROOF.md): the chain with TeX line numbers and each check.
- [euler_check.py](euler_check.py), [euler_factor.py](euler_factor.py), [checks.log](checks.log): local sums, Euler factor, L3.4, planted errors.
- [exponents.py](exponents.py), [exponents.log](exponents.log): exact exponent bookkeeping.
- Downloads (TeX, Zenodo texts, KP84/erratum scans, GL source) deleted after use.

RH relevance: N11
