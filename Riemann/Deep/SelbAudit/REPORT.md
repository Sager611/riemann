AUDIT PASSED — every new step (1)–(8) of the SelbergN4 record (κ_F ≥ 6.8158·10⁻⁴ at N = 4, order c/N) re-proved line by line with no error found, given Dousselin's analytic inputs; §10 needed explicit uniformity statements, now written and reviewer-confirmed.

# selbaudit — referee report on discoveries/selberg-linear-combinations.md (Riemann/Deep/SelbergN4)

## Verdict
PASSED. Complete proofs: FULLPROOF.md (i)–(viii), §9, code audit, numerics. No false step, no constant changes. Three gaps in the PRESENTATION of DERIVATION.md were found and closed:
1. §10 states that every δ-sum has the one-factor form Σ τ(δ)τ(δ')f(u_δ)g(u_δ')…; the large-v sum (§8) does not (J couples δ1, δ2). Closed: J(δ1,δ2) = J(δ_{11},δ_{21}) + o(1) uniformly for small-prime parts ≤ K,
   since J depends on them only through u (T-independent continuous profile d(x;u) on [1,Aϰ]) and a phase (δ_{10}δ_{20})^{−iv}, v ≤ H; the local sums then separate (FULLPROOF (viii)).
2. §10's "by uniform continuity … o(1/φ(d_1))" needs (R3) a T-independent modulus for the NORMALISED profiles, (R4) domination of the freezing errors by the old absolute (summable) majorant,
   (R6) the tail normalised as L·Σ1,tail ≪ tail(D_0) + o(1). All hold (the profiles √L·b, √L·δ^{−iv}s, F4_X, d(x;u) are explicit T-free continuous functions on compacts). Reviewer 1 (on the claim as
   originally stated): GAP on exactly these points; reviewer 2 (with the closure): CORRECT, "in the stated iterated limits — not an unrestricted simultaneous limit" (logs/review_viii.log).
3. Minor: |b(δ)| = τ_{1/2}(δ)… should read ≤ (|τ_{−1/2}(p^a)| = τ_{1/2}(p^a)/(2a−1)); §4.4 quotes the derivative error once normalised and once unnormalised (consistent); Dousselin's lastone
   (l.1671–1676) moves ∏_{p|q}(1−p^{−1−iv})^{−1} out of the modulus of a v-integral (not an identity; harmless since it is (q/φ(q))(1+O_q(v)), cost O_q(A³/log T)).

## Step by step (details: FULLPROOF.md)
| Step | Claim | Result |
|---|---|---|
| (i) §1 | m ≥ meas(G)/H, main term 1/A | CORRECT: G ⊂ ∪(γ_i − H, γ_i); Dousselin's pairing (l.736–741) loses exactly the factor 2 in the main term |
| (ii) §2 | N(√C1+√C2)² | CORRECT: (I*), (J*) are pointwise on S_j^* (l.886–898), so |I_j| + |M_j| ≥ H(1−o(1)) on S'_j; Cauchy–Schwarz; c_j fixed |
| (iii) §3 | C2 = 4(C3−1)/(θϰ)² | CORRECT: Minkowski; first moment ∫_T^{2T}Lη² = T + o(T) by the rectangle [1/2,3/2]×[T,2T] (convexity bound suffices); Dousselin's Gabriel step is monotone in its input |
| (iv) §4 | exact main terms, C5' = 2/√(πϰ(1−θ)) | CORRECT: closed forms re-derived; exact δ-structure B_d(δ) = τ_{−1/2}(δ)G_{qd}(ξ/δ); Im part via ∂_u and iH_2' = LH_2 − 2H_3 (error o(√log T) uniformly) |
| (v) §5 | Lemma Φ, Lemma W | CORRECT; Φ_BB(t) = min(1,(1−t/2)/(1−θ)) exactly (cert's step bounds are above it) |
| (vi) §7 | local factors cancel (P2 → P1) | CORRECT: identity verified term by term; ε_q(v) = 1 + O_q(v) |
| (vii) §8 | IBP twice, C1 linear in A | CORRECT: exact identity; d bounded, |d'| ≤ 2Psimax/x², log ξ/Λ ≤ ϰr(t_d); this removes Dousselin's ∫(H−v)/v = A log A |
| (viii) §10 | signed local sums −1/p, 0; P1 → P1eff | CORRECT after the closures above; order: T → ∞, then K, then D_0 → ∞, P_0 fixed |
| §9 | negative terms Sb, −aS1 | CORRECT: they multiply (A − X/ϰ) > 0 and (Aϰ − 1) > 0 with the same upper bound w̄ in both parts (code l.316–317) |

Uniformity: q, χ_j, c_j, N, θ, ϰ, A fixed; every o(1) is T → ∞; O_q- and O_{qP}-constants (P = ∏_{p≤10⁶}p) are fixed before the limit. κ_F is a liminf; no uniformity in q (as the page says).

## Falsification attempts (none succeeded; logs/)
- σ_p(k,0) = −1/p, 0 to 6 digits (p ≤ 13, k ≤ 4); drift ∝ v (deltasums.log).
- Exact A_d/A_1 at ξ = 10⁴, 10⁵, 10⁶, q = 1 and 13: d = 2: −1.28, −1.24, −1.21 (limit −1; the absolute bound allows 4.53); d = 4, 8, 9, 12 → 0; d = 3, 5, 6, 10, 30 converge slowly
  (still 20–60% off at ξ = 10⁶). Limit forced by dominated convergence. Rate caveat: the asymptotic regime is far away; irrelevant for a liminf.
- Exact Σ1·log ξ (validated against brute force): 0.818 (θ = 0.15, ξ = 2000; 0.880 for q = 13) against the new bound 1.945 and the old 22.4 (sigma1.log).
- J(B) with exact d(x;u), B ≤ 1290, ρ ∈ {2,3,6}: |J| ≤ 0.454 × IBP bound; |J|/B → constant (ibpcheck.log). Exact finite Dirichlet sums vs the profile: ≤ 8.6% at ξ = 10⁶ (exactsums.log).
- Counting toy: m ≥ meas(G)/H always; pairing up to 1.8× worse (counting.log). Special-function constants (mpmath) lie below the certified ones.

## Certificate
Rerun on a copy in /tmp (memcap 290 MB, one thread): identical cert.log up to timings, identical result.json; 262 s, peak 95 MB (logs/cert_rerun.log). Part 2 read against DERIVATION §4–§11:
formulas match, sups over boxes, upward rounding, κ_F = lower end of the Arb ball. Only non-rigour: float midpoints decide the peak test and Gsin's branch (wrong only within ~1e−30 of
(n+½)π or nπ; effect O(radius²)). Peaks of the audit scripts: deltasums 189 MB, sigma1 156 MB, exactsums 92 MB, counting 33 MB, ibpcheck < 100 MB (all < 290 MB).

## Not covered (as on the page)
Dousselin's analytic inputs (dominance sets and the Tsang-type CLT, Bauer's twisted moment, eqIntfacile, Lemma 12/eqintero, Gabriel step); odd characters ("similarly", l.640). R416 authenticated
the quoted inputs; this audit did not repeat that. If §10 were rejected the certified chain would give 9.73·10⁻⁵ (N = 4), still 2.6·10⁴ × Dousselin — it is not rejected.

## Files
Riemann/Deep/SelbAudit/: REPORT.md, FULLPROOF.md, deltasums.py, sigma1.py, ibpcheck.py, counting.py, exactsums.py, logs/ (cert_rerun.log, *.log, review_viii.log).

RH relevance: N12, N8
