AUDIT PASSED: Λ > 3.367 (= 3367/1000) holds as written — Lemmas 3–4 re-derive and match exact5/exact6 exactly, cert5.py (J = 5, 6) and certJ5.py rerun identically, certJ5.py is certJ.py plus a valid Lemma T′, H > 0 survives a second truncation (D = 70), every tail is inside H's radius (H/rad = 2·10⁵), and the new moments ⟨S^m⟩, m ≤ 12, match an independent BM-TeX integral to 8·10⁻¹⁴.

# gapcertaudit follow-up: explore53's J = 5 and J = 6 records (Λ > 3.358, 3.365, 3.367; node N8)

Claim audited: Riemann/Deep/ExploreShape/CERT.md, sections "UPDATE … 3.358", "… 3.365", "… 3.367", "Reproduce the record",
"Code versions and reproducibility", "Reviews". Audited files (sha256 prefix): cert5.py 0271dc508c4b5aa8, certJ5.py
49b566986a1077cf, certJ.py 9023b0588d88a589 (= HEAD), exact5.py 940bc3a0840718f2, exact6.py a2da4920ed882b7e (= the files
committed in 60a0882, checked with git show HEAD:… | shasum). Runs on a copy (.scratch/gapcertaudit/run2/ES, ../GapBar symlinked
read-only; deleted after use). Derivations: FULLPROOF.md, section "Follow-up" (F1)-(F3). Unchanged from REPORT.md: the
reduction, the Cauchy argument for the mean values (it is j-free), GapBar's reduction behind the Arb main terms.

## Verdicts by item
| Item | Verdict | Evidence |
|---|---|---|
| (1) Lemma 4 (J = 6) and Lemma 3 (J = 5): per-mode form, zero-mode matrix C, sufficiency (p > 0 on [1, ∞), S_odd ⪰ 0, S_even + C/2 ⪰ 0) | holds | own derivation (FULLPROOF F1); lemma46.py (sympy, from the coefficient maps of f, …, f⁽⁶⁾): q_d, q_e, γ_dd, γ_de, γ_ee, C identical to exact6.py; N11, N12, N22, DET of the certified form identical (exact) to numer6; J = 5 numerators identical to numer5 |
| (1) "CORRECT given deg p = 6 (a6 > 0)" | holds, and is checked | cert5.py's check p(1) > 0 with no real root ≥ 1 forces a6 > 0 (a6 = 1.44·10⁻⁴) |
| (1) cert5.py JJ = 6 (code) | holds | FULLPROOF F2: exact rationals; C11 > 0, det C = 3.1·10⁻¹¹ > 0; X* = max upper end of the possibly-real roots of N11, N22, DET, leading coefficients > 0; cover starts at (π/√X*).lower() ≤ π/√X* and ends at π.upper(); k ≤ 61, dropped modes PSD (x_k ≥ 62²(π/π.upper())² ≫ X*); psd = S11 > 0, S22 > 0, det > 0 on each Arb cell |
| (1) cert5.py JJ = 6 rerun (step 2) | identical | JSON identical; X* = 6.809387, 3253 cells, min det S_odd ≥ 4.744·10⁻¹¹ (g ≈ 1.657), min det(S_even + C/2) ≥ 9.52·10⁻¹¹ (g ≈ 2.681); CERTIFIED (13 s, 18 MB) |
| (1) J = 5 step: cert5.py JJ = 5 rerun on logs/scanJ5m_r2_d80_k3.358.json | identical | JSON identical; X* = 6.492941, 3748 cells, min det S_odd ≥ 4.03·10⁻¹³, min det(S_even + C/2) ≥ 4.83·10⁻⁷, CERTIFIED (J = 5); its main terms (HEAD certJ.py, K = 10) not rerun (code and Lemma T audited in REPORT.md) |
| (2) mean values j = 5, 6 | holds | the Cauchy argument of REPORT.md (b) is j-free: (v − S)^{2j} multiplier exact, error O_j(TL^{(r+2)²+2j−1}) |
| (2) ⟨S^m⟩, m ≤ 12: indep8.py (MMAX = 12, BM's TeX only) vs certJ5.py, test amplifier d = 8, κ = 13/4, J = 6 | holds | D = 72: every m ≤ 12 agrees to ≤ 7.5·10⁻¹⁴ (m = 11, 12: 6.6·10⁻¹⁴, 7.5·10⁻¹⁴; float limit); D = 60: all values inside the (wider) balls (365 s / 251 s, 158 / 100 MB) |
| (3) certJ5.py vs certJ.py | only the tail block changed | diff: docstring + l.150-175 (base 3c⁺/2, A_M = (2r+M)^n, × V_n); arrays, s_m, c_j, H, decision byte-identical |
| (3) Lemma T′ | holds | own derivation (FULLPROOF F3): \|h_k^{(n)}(0)\| ≤ c^{n−1−k}(c+n)^k ≤ c^n(n+1)^k (z = ±ic), ≤ k! (z = 0); Leibniz, Chu–Vandermonde, sh ≤ (M+N): \|H_n\| ≤ 4V_n(M+N+1)^n(3c/2)^{M+N}CC; reviewer: CORRECT |
| (3) Lemma T′ entrywise (majcheck.py D = 30, K = 12, exact entries M, N < 30) | holds | max \|H_n[M][N]\|/T′-majorant = 0.49 (n = 0) down to 2.9·10⁻²⁰ (n = 12), all at (0, 0) (26 s, 34 MB) |
| (3) "HS bound needs M ≥ 1, M = 0 by the operator norm" | right, harmless | ‖I⁰‖ = 1 = ϑ⁰/0!; for M ≥ 1 Young's inequality gives ‖I^M‖ ≤ ϑ^M/M! directly (HS: needs 3M² ≥ 2M) |
| (4) step 4 rerun (certJ5.py, D = 72, 397 bits) | identical | all fields equal the agent's JSON except timing: H = 1.14360·10⁻³⁷ ± 5.58·10⁻⁴³, H/c0 = 1.43562·10⁻⁵ ± 5.57·10⁻¹¹, CERTIFIED (643 s on a loaded machine, 192 MB memcap, 239 MB RSS) |
| (4) second truncation/precision: certJ5.py at D = 70 (392 bits) | holds | H = 1.144·10⁻³⁷ ± 7.50·10⁻⁴¹, CERTIFIED; overlaps D = 72 (\|ΔH\| = 4.0·10⁻⁴¹ ≤ sum of radii 7.56·10⁻⁴¹); every s_m ball overlaps (s_12: Δ = 4.2·10⁻¹⁴ rel. < D = 70 radius 8.8·10⁻¹⁴) (431 s, 148 MB memcap, 165 MB RSS) |
| (5) margin: every term of H, tails included, inside the radius | holds | code (certJ5 = certJ there): sb_m = s_m + [0 ± errs_m] for every m ≤ 12 before c_j = ⟨(v−S)^{2j}⟩ and H = −Σa_jκ^{2j}c_j, all Arb with exact a, κ, v; terms −a_jκ^{2j}c_j = 7.445, −8.250, −0.0993, 2.384, −3.402, 3.162, −1.240 (·10⁻³³, j = 0..6), sum 1.14359·10⁻³⁷ = H: \|H\|/max term = 1.4·10⁻⁵; rad(H) = 5.58·10⁻⁴³ (dominated by c_6's ball, 1.8·10⁻⁴⁵ ≈ 4.5·10⁻¹⁰ rel.), H/rad(H) = 2.05·10⁵ |

## Findings (none affects the certificates)
1. Reproducibility: during the audit the J ≥ 5 code was untracked; commit 60a0882 has since added cert5.py, certJ5.py,
   exact5.py, exact6.py and the record's logs, byte-identical to the audited files (sha256 above): resolved.
2. Documentation: CERT.md "UPDATE … 3.367" says "certJ.py (K = 12, D = 72, Lemma T′)"; the run used certJ5.py (as "Reproduce
   the record" step 4 and "Code versions" say). cert5.py's docstring still describes only J = 5 (Lemma 3, C = [[a3+φ32, φ34/2], …]);
   with JJ = 6 it calls numer6 and the J = 6 matrix C (checked).
3. Lemma T′'s factor max(1, k/c)^k is never active for k ≤ 12 (k! dominates), and the (n+1)^k bound already holds with T_k = 1
   for z = ±ic: the majorant is valid as coded and more conservative than needed.
4. Unlike the J = 4 form (REPORT.md: the a0 shift was needed), the J = 6 form is certified as the SDP produced it (D0 = E12 =
   E32 = 0): its small-gap Schur LMIs were imposed in the SDP.

## Review
- reviewer (Lemma T′, exact text inline: the h-derivative bound, Leibniz, Chu–Vandermonde, the sh step, the inner-product
  bound, the tail splitting and the 700 cut): **CORRECT** ("Steps 1–4 hold for all nonnegative integer orders, not merely
  those ≤ 12"; sharper sh bound 4U_mBA^{s−j}(s!/(s−j)!)(s−j+1)^m; the 10⁻³⁰⁰ pad is safe: a_700 < 10⁻³⁷⁵, ratio < 1/40).
- Lemma 4's algebra was not sent to a reviewer: it is checked symbolically (lemma46.py) against my own derivation, and an
  earlier reviewer (CERT.md "Reviews") found it CORRECT given deg p = 6.

## What survives
Λ > 3367/1000 = 3.367 unconditionally (J = 6, degree 80): the per-gap inequality of order 6 (Lemma 4, cert5.py JJ = 6) for
every g ∈ (0, π] and H > 0 (certJ5.py, Lemma T′) at two truncations; Λ > 3.365 and the J = 5 steps (3.355, 3.358) rest on the
same lemmas and code (Lemma 3 = Lemma 4 with a6 = φ52 = φ54 = 0, checked; the J = 5 main terms use HEAD certJ.py, whose Lemma
T is the one audited in REPORT.md). Dependents to update after commit: RECORDS.md N8 line, DISCOVERIES.md l.20 and the page
discoveries/large-gaps-higher-order.md (3.34 → 3.367), TREE.md N8.

## Files
lemma46.py (sympy check of Lemmas 3–4 against exact5/exact6), indep8.py (now MMAX env: ⟨S^m⟩, m ≤ 12), margin.py (terms and
radius of H; two-truncation comparison), majcheck.py (K argument; Lemma T′ entrywise), logs/ (prefix f2_: every follow-up run
with memcap peaks). Commands, from a copy of ExploreShape with ../GapBar symlinked: JJ=6 SRC=logs/scanJ6m_r2_d80_k3.367.json
uv run … cert5.py J6m80_3.367 0 0 0 61; uv run … certJ5.py logs/amp_d80_k3.367_J6.json (and its D = 70 copy); MMAX=12 indep8.py
2 0.2499 3.25 0.5 logs/amp_test_J4_D60.json 20; certJ5.py logs/amp_test_J6_D60.json; lemma46.py logs/cert5_J6m80_3.367.json.

RH relevance: N8
