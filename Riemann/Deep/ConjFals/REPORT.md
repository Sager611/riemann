NO COUNTEREXAMPLE: Jerby's Gram-point repulsion μ > 4 (2406.06548, 2310.14415; all 26,752,872 isolated bad Gram points with n < 2.54·10^8; smallest μ = 4.2155 at n = 69,784,845), Johnston–Leong–Tudzi's W_a conjecture (2408.04143; 17 values 2 ≤ a ≤ 32, x ≤ 2^30), Farr–Pauli–Saidak's Stieltjes bound (2301.11209; n ≤ 1000, a = j/16), Michałowski's global PF₄ (2602.20313; 1.1·10^7 random 4×4 configurations, local PF₄ on [0, 0.99]), Henry's inequality (2601.00631; 0 < x < 1), the 2511.02106 Lemma for 1 < κ < 3/2. Barrett–Jarvis's Redheffer conjecture, called open in 2511.13627 and misstated there, was disproved by Cardon (arXiv 0809.0076).

# ConjFals — falsification of stated conjectures on the RH tree

Orca agent conjfals, 11 Nov 2026, 3 h. Every computation single-threaded, < 300 MB, one at a time.

**Harvest** ([HARVEST.md](HARVEST.md)). 98 arXiv API queries → 441 papers 2023–2026; 105 e-prints scanned in TeX for conjecture/question/problem environments and conjecture phrases; Zenodo/HAL queries; one GXL evidence query for the Redheffer item. 12 rows of candidates ranked by chance × RH relevance: 8 finite, numerically testable statements tested below; the rest asymptotic, RH-conditional or with unspecified constants. Two conjectures found already settled by others (Barrett–Jarvis by Cardon 2008; Berndt–Straub's S_{2n} conjecture by 2610.18065).

## Summary

| Item | Range tested | Result |
|---|---|---|
| Jerby repulsion μ(g_n) > 4 at isolated bad Gram points | every Gram point n < 2.54·10^8 (41,773,404 bad, 26,752,872 isolated bad) | holds; smallest μ = 4.21553 at n = 69,784,845 |
| Sun \|W(x)\| < x, x ≥ 3078 | known to 2.5·10^14 (Mossinghoff–Trudgian) | holds; our octave sup 0.9777 at 2^30 |
| W_a: s_a < 1 for real a ≥ 2; lim s_a = 1 | a ∈ {2, 2.02, …, 32}, x ≤ 2^30 | holds (0.79 ≤ s_a ≤ 0.98); lim s_a = 1 proved (§2) |
| Barrett–Jarvis, as restated in 2511.13627 | n ≤ 262,144 | restatement false at n = 300; original disproved 2008 (Cardon) |
| Farr–Pauli–Saidak bound | n ≤ 1000, a = j/16 | holds, max ratio 0.627 (n = 1, a = 1) |
| Michałowski global PF₄ | 1.1·10^7 random configurations; H₄(u), u ∈ [0, 0.99] | no negative minor |
| Henry's two-sided inequality | x ∈ (0, 1), step 1/4000, 30 digits | holds; lower margin 0.00306 at x ≈ 0.485 |
| 2511.02106 Lemma x₂ < √(2x) for κ > 1 | κ ∈ [1.001, 1.49], y ∈ (2, 10^8] | holds |
| Spira's RH for sections | Lehmer pair γ₆₇₀₉, γ₆₇₁₀ | section has both zeros |

## 1. Gram-point repulsion (Jerby)

**Statement.** 2406.06548 (Nonlinearity 2024) conj:rep, TeX l. 313 and l. 777: if g_n is bad with good neighbours g_{n±1}, then \|Z'(g_n)\| > 4\|Z(g_n)\|. 2310.14415 Conjecture "Repulsion G-B-G", l. 266 (viscosity μ(g_n) = Z'(g_n)/Z(g_n), "non-sharp bound \|μ(g_n)\| > 4"), and Conjecture G-B-G, l. 753 ("corrupt" bad points, μ below the bound, are never isolated). Restated in 2610.21952 l. 323 as holding "to height of up to g_n of index n ~ 2.5·10^7", smallest value 4.26178 in the first 250,000 Gram points. Definitions re-read in the TeX: θ(g_n) = πn; g_n good iff (−1)^n Z(g_n) > 0; isolated = middle point of a Gram block of length 2.

**Test.** At a Gram point ζ(½ + ig_n) is real, Z(g_n) = (−1)^n ζ(½ + ig_n) and \|Z'(g_n)\| = \|Im ζ'(½ + ig_n)\|. (i) n < 10^6: Arb ball arithmetic at every Gram point (python-flint acb.zeta, acb_series; `gram_rep.py`). (ii) 10^6 ≤ n < 2.54·10^8: double-precision Riemann–Siegel main sum + C₀ (`gram_fast.c`, `gram_fast2.c`); on 2·10^4 ≤ n < 10^6 it reproduces the Arb run exactly (same 128,264 bad and 98,523 isolated bad points, same 123 points with μ < 5, \|Δμ\| ≤ 4.4·10⁻⁴). (iii) Range minima re-certified in Arb at 200 bits (interval Newton for g_n, enclosures of Z and Z') and recomputed independently with mpmath siegelz (`certify_gram.py`, `certify_minima.log`).

**Result.** No isolated bad Gram point with μ ≤ 4:

| range of n | isolated bad | min μ | at n |
|---|---|---|---|
| [−1, 2·10^4) | 1,680 | 4.26161122355266 | 8,569 |
| [2·10^4, 10^6) | 98,523 | 4.40850 | 45,924 |
| [10^6, 2.5·10^7) | 2,502,603 | 4.39258010410489 | 4,579,680 |
| [2.5·10^7, 10^8) | 7,890,681 | 4.21553180367418 | 69,784,845 |
| [10^8, 1.1·10^8) | 1,054,449 | 4.43906 (screening) | 107,532,018 |
| [1.1·10^8, 2·10^8) | 9,499,516 | 4.36820972726392 | 191,116,406 |
| [2·10^8, 2.54·10^8) | 5,705,420 | 4.42505873689157 | 228,357,901 |

Certified values carry Arb radii ≤ 4·10⁻¹⁵; mpmath agrees to all printed digits. Cross-check with the author's numbers: our Arb μ(g_730119) = 4.45773 against 4.46023 in 2310.14415 l. 1270, and our first-250,000 minimum 4.26161 (n = 8,569) against 4.26178 in 2610.21952 l. 351; the author's values are slightly less accurate, the indexing agrees (his corrupt g_9807962, in the Gram block 9807960–9807963, is bad and non-isolated here too, but its μ is 0.04250, not 0.0750883 as in 2310.14415 l. 733: Z'(g_n) ≈ −0.0085 is near zero, so approximate Z' values are far off). In each extreme case \|Z(g_n)\| ≈ 12–13 sits next to a neighbour with \|Z\| < 0.07 (a zero just beside g_{n±1}). The normalised minimum μ/θ'(g_n) per range (1.16, 0.80, 0.66, 0.55, 0.53, 0.54) falls nearly as fast as the threshold 4/θ'(g_n) (1.09 at n = 10^4, 0.48 at n = 2.5·10^8): the margin shrank from 0.10 (n ~ 10^6–10^7) to 0.05 (n ~ 2·10^8), so the bound looks empirical, not structural. A pair-plus-lattice model of Z'/Z = Σ 1/(t − γ) allows μ down to θ'(t)/π = (1/2π) log(t/2π), below 4 until t ≈ 5·10^11. Screening caveat: a Gram point with \|Z\| ≲ 10⁻⁵ could be misclassified in double precision; it has huge μ itself, so only its neighbours' isolation status could change.

**Novelty.** No refutation found (arXiv listings Oct–Nov 2026 for "Gram"/"viscosity"; Zenodo; the author's 2610.21952 restates the conjecture). Our exhaustive range 2.54·10^8 is 10 times the author's 2.5·10^7 (exhaustiveness not stated there); μ = 4.2155 lies below his smallest reported 4.26178.

RH relevance: N8 (Gram/Rosser counting of critical zeros) and N13 (Jerby proposes this repulsion as a mechanism behind RH); a counterexample would remove that mechanism, not touch RH.

## 2. Summatory function of (−a)^{Ω(n)} (Johnston–Leong–Tudzi; Sun)

**Statements.** 2408.04143 (Aug 2024) l. 128: Sun's conjecture \|W(x)\| < x for x ≥ 3078, verified to 2.5·10^14 by Mossinghoff–Trudgian (l. 127; \|W(x)\| < 0.979x on [1.25·10^14, 2.5·10^14], l. 649). Conjecture newcon2, l. 196: for every real a ≥ 2 and x ≥ x₀(a), \|W_a(x)\| < x^{log₂ a}, W_a(x) = Σ_{n≤x}(−a)^{Ω(n)}. Theorem 3thm, l. 217: s₃ = 0.813 ± 0.158, s_a = limsup \|W_a(x)\|/x^{log₂ a}. Conjecture, l. 1142: lim_{a→∞} s_a = 1.

**Test.** Segmented Ω-sieve, exact integers for a = 2 (`sunW.c`, n ≤ 2^32), long double for real a (`sunWa.c`, n ≤ 2^30, 17 values of a at once). Octave suprema of \|W_a(n)\|/n^{log₂ a} are stable to 4–5 digits over the last 6 octaves:

| a | 2 | 2.02 | 2.05 | 2.1 | 2.2 | 2.35 | 2.5 | 2.75 | 3 | 3.5 | 4 | 5 | 6 | 8 | 12 | 16 | 32 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| sup on [2^29, 2^30) | 0.9772 | 0.9680 | 0.9547 | 0.9338 | 0.8995 | 0.8691 | 0.8456 | 0.8227 | 0.8106 | 0.7979 | 0.7946 | 0.8005 | 0.8123 | 0.8371 | 0.8748 | 0.8992 | 0.9439 |

No a gives s_a ≥ 1; s₃ ≈ 0.8106 (inside and much sharper than 0.813 ± 0.158; numerics, not certified). For a = 2 the maximum sits at x ≈ 1.528·2^k, for large a at x = 2^k.

**lim s_a = 1 (proof).** Write n = 2^k m, m odd, g(m) = (−a)^{Ω(m)}, α = log₂ a > 1. Summing the geometric series in k, W_a(x) = (1/(1+a)) Σ_{m odd ≤ x} g(m) + (a/(1+a)) x^α Σ_{m odd ≤ x} g(m) m^{−α} h(log₂(x/m)), h(w) = (−1)^{⌊w⌋} a^{−{w}}. The first sum is O(x^{max(1, log₃ a)+ε}) = o(x^α), and Σ_{m odd} \|g(m)\| m^{−α} = Π_{p odd}(1 − a^{1−log₂ p})^{−1} < ∞ for a > 2. At x = 2^k the m-term equals (−1)^k (−1/a)^{j(m)}, j(m) = ⌈log₂ m⌉ − Ω(m) ≥ 1 for odd m > 1, so s_a ≥ (a/(1+a))(1 − Σ_{m odd > 1} a^{Ω(m) − log₂ m}); the sum tends to 0 as a → ∞ by dominated convergence (majorant 3^{Ω(m)} m^{−log₂ 3}, summable). With the upper bound of 2408.04143 Theorem limsupthm2 (tends to 1), lim s_a = 1. Check: a = 32 gives (32/33)(1 − 1/32 + 5/1024 − …) ≈ 0.944 (j = 2 for m = 5, 7, 9, 15, 27); sieve 0.9439. Reviewer agent: correct, elementary, resolves the conjecture of l. 1142. The same expansion gives s_a ≤ max((a−1)/(a+1), 2a^{−0.585}) + O(a^{−1.17}) (on [0, 0.585) the m = 3 term is −a^{−1−v}, the m = 1 term a^{−v}; every odd m ≥ 5 has Ω(m) − log₂ m ≤ −1.17), so newcon2 holds for every sufficiently large a (sketch; a₀ not computed, not reviewed).

**Novelty.** Sun's case: verified before (not ours). The latest version (v3, Sep 2024) still states both conjectures; no later paper on W_a (arXiv API, monthly listings, one GXL evidence query). The proof of lim s_a = 1 is elementary and off-tree (side result).

RH relevance: N10 (Mertens/Pólya-type size conjectures); the main term comes from the Euler factor at 2, not from zeros.

## 3. Barrett–Jarvis (Redheffer matrix)

**Statement.** 2511.13627 (Nov 2025) l. 524: "A remarkable conjecture states that the 'small' eigenvalues have negative real part; see [BJ]" (Barrett–Jarvis, LAA 162 (1992)). The original (quoted in Cardon, arXiv 0809.0076, TeX l. 287–294): the small nontrivial eigenvalues λ satisfy (i) \|λ\| < 1 and (ii) Re λ < 1.

**Test.** det(xI − R_n) = y^{n−s−1}(y^{s+1} − Σ_{k=1}^{s} F_k(n) y^{s−k}), y = x − 1, s = ⌊log₂ n⌋, F_k(n) = number of ordered factorisations into k factors ≥ 2 with product ≤ n (matrix determinant lemma on R_n = I + N + u e₁ᵀ); checked against numpy eigenvalues of R_n for n ≤ 300 (`redheffer_bj.py`). At n = 300 the small eigenvalues are 0.9649, 0.775 ± 0.263i, 0.303 ± 0.465i, −0.040 ± 0.294i, so the restatement is false; for n ≤ 262,144 the largest real part of a small eigenvalue is 0.99521 (n = 145,818), consistent with the original.

**Novelty.** Cardon disproved both parts of the original with n up to 2^36 (0809.0076, Theorem 5.3). So 2511.13627 misstates the conjecture and calls it open 17 years after its disproof. KNOWN; a misstatement only.

RH relevance: N10 (det R_n = M(n)).

## 4. Stieltjes-constant bound (Farr–Pauli–Saidak)

**Statement.** 2301.11209 (Pauli–Saidak, Jan 2023) l. 523–525, quoting Farr–Pauli–Saidak (Funct. Approx. 64, 2021): \|C_α(a)\| ≤ 2\|exp(α(log w(α) − 1/w(α)))\|, w(α) = W₀(iα/2π), C_α(a) = γ_α(a) − log^α(a)/a, 0 < a ≤ 1.

**Test.** C_n(a) = γ_n(1 + a) (no cancellation), Arb (acb.stieltjes) at 192 bits with a radius check, every n ≤ 300 and every 7th n up to 1000, a = j/16 (`fps_stieltjes.py`). A first run at 128 bits with a = 0.1 gave a spurious ratio 3.5 from cancellation in γ_n(a) − log^n(a)/a; the γ_n(1 + a) form removes it.

**Result.** Max ratio 0.6269 at n = 1, a = 1; at most 0.56 for n ≥ 2; no growth in n. Holds for integer α.

**Novelty.** No refutation found (arXiv API "Stieltjes constants" bound, listings).

RH relevance: N3 (Stieltjes constants enter λ_n).

## 5. Global PF₄ of the de Bruijn–Newman kernel (Michałowski)

**Statement.** 2602.20313 (Feb 2026) Open Problem 1, l. 820: prove or disprove K = Φ(\|u\|) ∈ PF₄ for all configurations (PF₅ fails at certified Toeplitz minors).

**Test.** (i) Coalescent limit: every r×r minor with x_i → x, y_j → y has the sign of H_r(u) = det[(−1)^j K^{(i+j)}(u)], u = x − y; mpmath, 60 digits, u ∈ [0, 0.99] step 0.01 (`pf_local.py`): H₂, H₃, H₄ > 0 throughout, H₅ < 0 exactly for u < 0.031 (their sign change at 0.0311 reproduced). (ii) 8·10^6 random configurations x, y ∈ [−0.8, 0.8]^4 with random-walk descent (`pf_search.py`): smallest det/Πdiag = 1.3·10⁻¹¹ > 0, descent drives to coalescence; 3·10^6 multiscale configurations with q = det/(V(x)V(y)Πdiag), V the Vandermonde product, finite in the coalescent limit (`pf_search2.py`): min q = 4.78 > 0.

**Result.** No negative 4×4 minor; the global question stays open.

**Novelty.** No later work on Problem 1 found (listings; Zenodo "Pólya frequency" notes concern coefficient sequences).

RH relevance: N2 (the de Bruijn–Newman kernel).

## 6. Henry's inequality

**Statement.** 2601.00631 Conjecture (l. 286) = 2601.09276 Conjecture (l. 113): for 0 < x < 1, π cot πx + x < ζ(x) − ψ(x) < π cot πx + (γ+½)x + γ − ½.

**Test.** mpmath, 30 digits, x = k/4000, plus endpoint expansions (`henry.py`): the lower gap tends to γ − ½ at 0 and 2γ − 1 at 1, minimum 0.003065 at x ≈ 0.485; the upper gap is ≈ 0.3512x at 0 and 0.6405(1 − x) at 1 (Laurent/Taylor series of ζ, ψ, cot).

**Result.** Holds; an interval-arithmetic proof looks routine (both endpoint slopes positive).

**Novelty.** 2601.09276 treats it as open; no refutation found.

RH relevance: none (ζ on (0, 1)).

## 7. Spira's RH for sections

**Statement.** 2405.12557 l. 164: all non-trivial zeros of Z_{N(t)}(t) = Σ_{k ≤ [t/2]} cos(θ(t) − t ln k)/√k are real.

**Test.** The function is discontinuous in t and "non-trivial zero" is undefined; only the natural test was run: at Lehmer's pair γ₆₇₀₉ = 7005.0629, γ₆₇₁₀ = 7005.1006 (max \|Z\| ≈ 0.004) the section Z₃₅₀₂ changes sign twice, shifted by about 0.005 (`spira.py`). The section error is about 1/√(2t) with an oscillating sign, so closer pairs of the opposite sign could defeat it; not pursued.

**Novelty.** —

RH relevance: N8.

## 8. 2511.02106 Remark (σ-analogues of Robin's criterion via LCM powers)

**Statement.** Remark after Lemma "upper bounds for x₂" (TeX l. 1285–1347): the Lemma x₂ < √(2x), proved for κ ≥ 3/2, "seems to hold for all κ > 1"; sufficient: [(x²+2)^κ − 2^κ][(x+1)^κ − 1]² − [(x²+x+1)^κ − (x+1)^κ]² increasing in x. By the proof (l. 1289–1292) the Lemma is equivalent to F(y²/2, 1) > F(y, 2) for every y > 2, F as in l. 1050.

**Test.** mpmath, 40 digits, κ ∈ {1.001, 1.01, 1.05, 1.1, 1.2, 1.3, 1.4, 1.49}, 281 log-spaced y ∈ (2, 2 + 10^8] (`lcm_lemma.py`).

**Result.** The difference stays positive (it tends to 0⁺ like log 2/(2y² log² y)); the sufficient function never decreases on x ∈ (0, 10^10). Holds.

**Novelty.** Not checked beyond the listings (no follow-up found).

RH relevance: N4 (Robin-type criteria).

## Not tested

Asymptotic or non-finite statements (HARVEST rows 9–12): Ihara's Euler–Kronecker bounds, κ(q) = o(1), 2608.24827's λ_min law (its prediction 542.6 vs the certified −ln λ_min = 556.6 at c = 47 is consistent), the effective LI conjecture with an unspecified constant, limit statements "except countably many a", RH-conditional and path-existence claims.

## Files

`gram_rep.py`, `gram_fast.c`, `gram_fast2.c`, `certify_gram.py`, `certify_minima.log`, `gram_test.log`, `gram_1e6.log`, `gf_*.log`; `sunW.c`, `sunWa.c`, `sunW_2e32.log`, `sunWa.log`; `redheffer_bj.py`; `fps_stieltjes.py`; `pf_local.py`, `pf_search.py`, `pf_search2.py`, `pf_verify.py`; `henry.py`; `spira.py`; `lcm_lemma.py`; `harvest/` (query lists, `harvest.py`, `dl_scan.py`, `scan.py`); [HARVEST.md](HARVEST.md).
