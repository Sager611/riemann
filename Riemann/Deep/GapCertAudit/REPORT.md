AUDIT PASSED: Λ > 3.34 (= 167/50) holds as written — the reduction, the order-3/4 mean values (BM's own Cauchy argument), Lemma 2's sufficiency and both Arb certificates check line by line, rerun identically, survive a second truncation/precision (D = 64), and the new moments c_3, c_4 match an independent integral from BM's TeX to 5·10⁻¹⁴.

# gapcertaudit: line-by-line referee of our discovery large-gaps-higher-order (node N8, falsification first)

Sources: BM = Bui–Milinovich arXiv 1410.3635v3 (TeX main.tex, 758 lines), BBLR = Bettin–Bui–Li–Radziwiłł arXiv 1609.02539
(Thm 1.1, l.103-136); both deleted after use. Claim audited: discoveries/large-gaps-higher-order.md, ExploreShape/CERT.md.
Full chain with line numbers: FULLPROOF.md. Code and logs: this folder (runs on a copy in .scratch, GapBar read-only).

## Verdicts by item
| Item | Verdict | Evidence |
|---|---|---|
| (1a) Reduction: Λ ≤ κ ⇒ gaps of F's real zeros ≤ (1+o(1))κπ/L; rescaling λ = (1+ε)κ/L; (V) on Re F, Im F; telescoping | holds | BM l.180-183 (gaps), l.194 (end pieces); brackets Re(F′F̄″), Re(F‴F̄″) are single smooth functions; FULLPROOF (a) |
| (1a) End and partial-gap terms | holds | |F^{(k)}|² ≪ T^{2/3+ϑ+ε}L^{2k} (Weyl + Cauchy + |M| ≪ T^{ϑ/2+ε}), length ≤ 2c/L: O(T^{11/12+ε}) = o(TL^{16}); reviewer 2: CORRECT |
| (1b) ∫|F^{(j)}|², j = 3, 4: operator identity | holds, exact for every j | F^{(j)} = (iL)^je^{ivtL}(v + L⁻¹Σ∂_α)^jG; c(α,β) (l.377-381) depends on shifts only through one exponential; own computation: L⁻¹Σ∂_α = L⁻¹Σ∂_β = −(log T/L)S on it (the (α2−α1), (β2−β1) terms cancel) ⇒ (v − S)^{2j} + O_j(1/L) (= BM l.422, l.431) |
| (1b) Error O_j(TL^{(r+2)²+2j−1}) | holds | E = I − Main holomorphic on the polydisc of radius 1/L around (α*, β*), ≪ TL^{(r+2)²−1} there (Lemma main l.383 "uniformly", complex shifts, BBLR l.104); Cauchy: (L⁻¹∂)^n costs O_n(1); times L^{2j} |
| (1b) Unsmoothing | holds (leading term) | w in [1,2], w^{(j)} ≪ T^ε (l.290): sandwich by finite partitions at scales T_i ≍ T, log T_i/L = 1 + O(1/L); BM's own h1 (l.200) uses the same step |
| (1c) Lemma 2, sufficiency: sine/cosine coefficients with free f″(0), f″(g); termwise minimisation; p > 0, S_odd ≥ 0, S_even + (a3+φ32)/2 ≥ 0 | holds | re-derived (FULLPROOF (c)): Q = (g/2)Σ[pb² + βbd + γd²] + (a3+φ32)(A1−A0)²/g, grouping per k needed (pb², γd² diverge separately); N's x⁴, x⁵ coefficients vanish identically; reviewer 1: lemma correct |
| (1d) certJ.py pipeline vs GapBar's audited reduction | holds | h_series_K/uhat_K = gcert's for any k; g_m = Σ C(m,k)conv2(u_k,w_l), H_n = ΣC(n,j)sh^j(g_{n−j}) reduce to gtw's a, b, c at n ≤ 2; ⟨S^m⟩ = ΣC(m,i)Q(H_{m−i},X^i) is the binomial expansion of commuting multiplications |
| (1d) Lemma T (truncation majorant) | holds (re-derived; even (2+M+N)^n suffices) | FULLPROOF (d); majcheck.py on the certificate's arrays (D = 72, n ≤ 8): every entry ≤ 0.49 × majorant, ≤ 1.9·10⁻⁴⁶ × majorant at max(M,N) ≥ 64 |
| (2) cert4.py | holds | exact rationals → arb; Arb-rigorous root isolation; cover starts at (π/√X_N).lower() ≤ π/√X_N, ends at π.upper(); k ≤ 61, dropped terms positive (x_k ≥ 62²·(π/π.upper())² > X_N) |
| (2) certJ.py | holds | exact decimals (GC.dec) and Fractions; float only picks β = 25 (≥ 24.74) and the precision; tails added to radii before c_j, H; a in the amplifier file = cert4's a (string-equal); φ12, φ32 do not enter H |
| (3) Reruns | identical | cert4: JSON identical, CERTIFIED; certJ D = 72: every field identical (H = 9.3142330897e-35 ± 1.92e-46) |
| (3) Independent H | holds | certJ at D = 64 (361 bits): H = 9.31e-35 ± 5.38e-38, CERTIFIED; all balls overlap D = 72's |
| (3) Independent c_3, c_4 | holds | indep8.py (BM's TeX only) vs certJ on the d = 8 test amplifier at D = 72: ⟨S^m⟩/⟨1⟩, m ≤ 8, agree to ≤ 5.5·10⁻¹⁴ |
| (3) Independent (V) at g = π, other basis (consistency: a finite basis bounds the true minimum from above; the proof is cert4) | holds | vbasis.py (sines k ≤ 30 + 8–12 polynomials, 50 digits): λ_min = 1.817086256·10⁻⁵ (stable to 10 digits); other g: 0.02407 (3.1), 0.08503 (3.0), 0.24289 (2.7707), 1.1020 (2), 7.2603 (1), 92.476 (0.3), 3406.1 (0.05) = orch_vnum's |

## The a_0 shift 2·10⁻⁵ (brief item 3), explained (shifts.py, vbasis.py; exact rationals, 40–50 digits)
| form | n3 | S_odd(π) | S_odd(10⁻³) | λ_min(π) = min Q/‖f‖² | valid? |
|---|---|---|---|---|---|
| sdp: SDP witness rationalised (a0, φ12 unshifted) | −7.8·10⁻¹⁰ | +1.2·10⁻⁹ | −7.3·10⁻¹⁶ | +1.92·10⁻⁹ | no: fails for gaps g ≲ 10⁻³ (ρ ~ n3/(4a4x) < 0 for x ≳ n2/|n3| ≈ 10⁷) |
| noshift: φ12 + 2·10⁻⁶ only | +1.98·10⁻⁸ | −1.15·10⁻⁶ | +2.5·10⁻¹³ | −1.829137·10⁻⁶ | no: fails at g = π |
| cert: φ12 + 2·10⁻⁶, a0 + 2·10⁻⁵ (certified) | +1.98·10⁻⁸ | +1.14·10⁻⁵ | +2.5·10⁻¹³ | +1.817086·10⁻⁵ | yes (cert4.py, Arb) |
- λ_min(cert) − λ_min(noshift) = 2·10⁻⁵ exactly (a0 shifts the Rayleigh quotient by a constant): the orchestrator's 1.82·10⁻⁵
  is this λ_min. So the a0 shift is not a mere safety margin: the φ12 raise needed for small gaps (n3 > 0) breaks g = π by
  1.83·10⁻⁶, and the a0 raise repairs it with margin 1.82·10⁻⁵. CERT.md/REPORT.md describe each shift correctly
  ("only helps", "makes n3 > 0") but not that both are necessary. Cost: H/c0 drops by exactly 2·10⁻⁵ (0.001383 → 0.001363).

## Findings (none affects the certificate)
1. Documentation: the page, CERT.md "Checks" and R494 (4) say crosscheck.py matches "the j = 3, 4 moments (0.0088906238,
   0.020059772)". crosscheck.py has K = 4 (l.52): these are ⟨(v−S)^3⟩, ⟨(v−S)^4⟩ (r = 1, no twist), so only c_0, c_1, c_2
   were cross-checked; c_3 = ⟨(v−S)^6⟩ and c_4 = ⟨(v−S)^8⟩ had no independent check. Now filled: indep8.py (above).
2. The 2·10⁻⁵ "safety shift" is a necessary repair (section above); the unshifted SDP witness is invalid for tiny gaps.
3. Lemma T holds entrywise on the certificate's arrays (majcheck.py: max ratio 0.49 at (M, N) = (0, 0); ≤ 1.9·10⁻⁴⁶ at
   the truncation edge max(M, N) ≥ 64) and is very conservative there; the D = 60/64 radii are this bound (unchanged with
   +160 bits, hiprec.py). At D = 72 the whole radius of H is 2·10⁻¹² of H. The sharper (2+M+N)^n also suffices.
4. Cosmetic: cert4.py l.66 comment "x_k >= 4 > X_N" is stale (K = 61 is used; the argument holds); BM l.428 typo
   "(β2−β1t1)t4" for (β2−β1)t1t4 (l.378 is right); BM l.395 divides by ie^{2ivθ(t)} (GapAudit's (a)).
6. Off the audited claim: Lemma 0 (orders ≤ 2 add nothing) spot-checked only (for J = 2 the only bracket [f′f] vanishes;
   m2² ≤ m0m4 by Cauchy–Schwarz; the moment cone of [1, ∞)); J = 5, 6 numbers are float, not audited.
5. Unsmoothing is only sketched in CERT.md ("sandwiching"); FULLPROOF (b)5 writes it out (finite partitions at scales
   T_i ≍ T; BM's w must be supported in [1, 2]). Leading term only, which is all the reduction needs.

## Review
- reviewer 1 (Lemma 2 + its use, exact formulas inline): Lemma 2 steps 1–4 valid ("no divergent series is separated", ρ =
  O(1/x), parity factors and the constant (a3+φ32)(A1−A0)²/g correct); GAP only for the end intervals E = [T, t̃_1] ∪ [t̃_N, 2T],
  which the prompt had omitted.
- reviewer 2 (that closure: pointwise bound |F^{(k)}|² ≪ T^{2/3+ϑ+ε}L^{2k}, |E| ≤ 4c/L, BM l.194): CORRECT ("every weighted
  end-interval integral and both weighted boundary products are O(T^ρ/L) = o(TL^{16})", ρ = 2/3 + ϑ + η < 1).

## Scope
- Audited code = git HEAD: my run copies of cert4.py and certJ.py are byte-identical to the committed files (git show
  HEAD:… | diff). explore53 has since edited certJ.py in the working tree (uncommitted "Lemma T′": base 3c/2, factors
  T_k, V_n, (2r+M)^n) for its J = 5 work: NOT audited, not used by the J = 4 certificate; it needs its own audit (and
  the J = 4 logs must keep pointing to the HEAD version) before any J = 5 certificate is built on it.

## What survives
Everything claimed: Λ > 3.34 unconditionally; (V) for all g ∈ (0, π] (Arb); H > 0 (Arb, two truncations). Dependents
(no change needed beyond finding 1's wording): RECORDS.md l.24, DISCOVERIES.md l.20, NOGO.md l.694, TREE.md l.530, R494.lean.

## Files
FULLPROOF.md (chain, source lines); indep8.py (BM-TeX-only ⟨S^m⟩, m ≤ 8); vbasis.py (λ_min of (V), sines + polynomials);
shifts.py (the two shifts); majcheck.py (entrywise Lemma T); hiprec.py (certJ with extra bits); logs/ (all runs, memcap peaks). Commands: in a copy of ExploreShape (.scratch/gapcertaudit/run/ES,
../GapBar symlinked; deleted after use), bash memcap.sh 290 1200 uv run --no-project --with python-flint python3 certJ.py logs/amp_d64_k3.34.json;
same with logs/amp_d64_k3.34_D64.json (D = 64) and logs/amp_test_J4_D72.json; indep8.py 2 0.2499 3.25 0.5
logs/amp_test_J4_D60.json 20; vbasis.py pi 20,30 8; shifts.py.

RH relevance: N8
