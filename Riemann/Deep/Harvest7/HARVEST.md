CANDIDATES: 3 (best: arXiv 1606.05011 = Dimitrov–Xu, TAMS 372 (2019) Thm 1.1, flagged by Zenodo 24215530; N6/LP-class RH criterion: as printed it is equivalent to not-RH; Arb certificate done here)

Scope: the arXiv announcement of 13 Nov 2026 (375 ids ≥ 2611.14518 in NT, CA, CV, SP, FA, PR, CO, math-ph, GM; 368 more in DS, NA, OA, RT, AG, cs.LO, hep-th, quant-ph by keyword), 10 arXiv API keyword queries by submission date, Zenodo records created 12–14 Nov (zsweep: 23; four extra query sets: about 150 more), 26 TeX sources scanned with rg. Harvests 1–6 and rounds R456–R467 skipped. Working log: [LOG.md](LOG.md).

Yield: a thin arXiv week. No arXiv paper of 12–14 Nov leaves a node constant unstated or rests on an input our tools supply (16550 has dropped MRH, 16828 is asymptotic, 16156 is off-tree). What matters comes from Zenodo: one false published RH criterion (rank 1, certified here) and two N8 record lines out of date (ranks 2–3; also flags 1–2).

## Ranked table

| # | id + TeX line | node | statement | what our tools give | gain | risk | pre-emption check |
|---|---|---|---|---|---|---|---|
| 1 | arXiv 1606.05011 l.198–207 (= TAMS Thm 1.1 per Koide l.67); Zenodo 24215530 (A. Koide, 14 Nov) l.42, 69, 432–451 | N6 / N2 (Laguerre–Pólya class) | "RH ⟺ for each y ∈ (−1/2,1/2)∖{0} the translates of Φ_{2,y}(t) := cosh(ty)∫(t−2s)²Φ(t−s)Φ(s)ds are dense in L¹(ℝ)"; 2nd clause: density for every y ⟺ zeros real and simple | F[Φ_{2,y}](x) = 2Re(Ξ′² − ΞΞ″)(x+iy). Arb jets of ξ (dxcert.py, 160 bits): at y = 49/100, R(304.60) > 0 > R(305.00), so there is a real zero, the translates are not dense (Wiener) and the RHS fails. Printed Thm 1.1 ⟺ ¬RH; its 2nd clause fails likewise. Jensen's corrected J = \|Ξ′\|² − Re(Ξ·conj Ξ″) > 0 at the same points (sanity check). mpmath scans: sign changes at y = 0.25, 0.40, 0.49 by height 10⁴, and at y = 0.499 already at x ≈ 110.45 (Arb: R(110.40) > 0 > R(110.50)). At y = 0.05 and 0.10 there is none in [10⁴, 10⁴ + 4] | a REFUTATIONS entry: a refereed RH criterion is false unless RH is false. Koide's corrected criterion is Jensen 1913 + Wiener (classical) | low. TAMS print vs arXiv v1: Koide quotes TAMS with the same kernel; check the TAMS PDF (AMS opens it after 5 years) | Koide proves Thm 3.2 false with a constructed LP kernel, K̂ = ρ_ε cos(Mz)cos(2Mz). He does not show that Thm 1.1 fails for Φ. Zenodo (6 queries, no date filter): nothing else. arXiv keyless search: "rejected as irrelevant" |
| 2 | Zenodo 24166562 (R. Khan, A. DeWeese, 12 Nov 17:32 UTC): q12 CandidateSpec; Q12/Simple.lean l.92–98 | N8 | dyadic and cumulative (κ − ε)·Ncount ≤ N0star, and the same with N0simple, for κ = 854046308300026565326564954594528531931648/1267673952595351605232565361802825927734375 = 0.6737113329 (q11: 0.6736952871; q10: 0.6736701854) | cannot be replayed under our cap (export 1.75 GB, nanoda and preplay ≈ 14 GB, 4–6 h). Statement audit done here (LOG.md): the trusted ChallengeDeps files are byte-identical to the arena's pinned upstream zeta-23-lean 3635e748; the three theorems are verbatim the arena template, with no hypotheses; no sorry, axiom or native_decide in the 146 Q12 files or in Zeta23/. Logs: standard axioms, nanoda rc 0, comparator preplay rc 0 | record +7.20·10⁻⁵ over the arena's 0.6736392880, for distinct and for simple zeros (records only) | the kernels were run by the authors only; no hosted check in the archive; the arena record was still 67.3639287955% at 18:36 UTC (its ledger lists only q9 = Khan, 25 Oct, 0.6736387063) | not in our files (rg of the id and of 0.67371: none) |
| 3 | anthropics/zeta-23-lean commit c77974d8 (10 Aug 2026), comparator/Challenge/XiPrime.lean l.24–52; Zenodo 24210376 (M. Iizumi, 14 Nov 00:52 UTC) README l.10–19 | N8 sibling (zeros of ξ′) | zeta23: ≥ 0.85838 of the zeros of ξ′ simple and on the line (0.86864 with the quartic window), ≥ 0.92919 distinct (0.93432), dyadic and cumulative, standard axioms, comparator "okay". Iizumi: > 0.87328 simple and on the line, using zeta23's XiPrime counts | statement audit done here (LOG.md): his theorem uses the library's Zeta23.XiPrime.Ncount and N0simple at the arena's trusted commit 3635e748; their bodies match the comparator's ChallengeDeps/XiPrime.lean (ξ = s(s−1)/2·Λ₀ + 1/2, multiplicities by analyticOrderAt). His compiler log prints only the standard axioms for all four final theorems | corrects our N8 line: R392's "ξ′ proportion 0.8137 (Conrey) unchanged, no 2025–2026 improvement" has been false since 10 Aug | Iizumi's clean build (13 Nov 20:04 JST) is his own record; no third-party kernel run. zeta23's is Anthropic's comparator run | ours: R392 missed zeta23; R425 saw only Iizumi's conditional 88.45% (24101490, Weil + Deshouillers–Iwaniec hypotheses; its HOLDS stands) |

Side results (not on a node) are listed under Routine one-liners below.

RH relevance:
- #1: low to moderate. It removes a false published RH-equivalent; the correct criterion (Jensen + Wiener) is classical. Chance 0.95, cost about 30 min (reviewer + TAMS text). Rank 1 because it is cheap, certified and a false step on a node.
- #2: low. It is a record inside a method whose ceiling is ≤ 0.6818287 (zeta23 PairCeiling, given EnclOK) and 153/227 for TEN-type windows (cglunc). Statement audit passed here (steps 1–4); the remaining risk is the unreplayed kernels. Chance 0.85; gain: a record line.
- #3: low. Zeros of ξ′ on the line do not imply RH. Chance 0.9; gain: a corrected sibling record.

## Mini-briefs
Only three candidates qualify. Ranks 4–5 would be the side one-liners below, which need no Orca agent.

### 1. Dimitrov–Xu's printed RH criterion ⟺ ¬RH (N6/N2; orchestrator + reviewer, no Orca needed)
- Goal: a REFUTATIONS entry. TAMS 372 (2019) 4107–4125, Thm 1.1 as printed cannot hold unless RH fails, and its proof is broken: Lemma 3.3 drops the factor −i of the odd-kernel transform (Koide). The general Thm 3.2 is false (Koide's counterexample).
- Inputs: arXiv 1606.05011 l.196–207 (statement; re-download https://arxiv.org/e-print/1606.05011, a single gzipped .tex) and the TAMS PDF (to confirm the kernel cosh(ty)ν₂); Koide 24215530 (.tex on Zenodo, l.42–116, 432–451, 828–846).
- Method: ν̂₂(w) = 2(Ξ′(w)² − Ξ(w)Ξ″(w)), since ∫u²Φe^{ixu} = −Ξ″ and ∫uΦe^{ixu} = −iΞ′. Then cosh(ty) shifts to (ν̂₂(x+iy) + ν̂₂(x−iy))/2 = Re ν̂₂(x+iy), and Ξ′² − ΞΞ″ = ξξ″ − ξ′² at s = 1/2 − y + ix. Opposite Arb signs at two x for one y give a real zero; by Wiener's L¹ theorem the translates are not dense.
- Certificate: [dxcert.py](dxcert.py), [dxcert.out](dxcert.out): y = 49/100, sign R(304.60) = +, R(304.70) = −, R(305.00) = −, R(305.40) = +; J/|ξ|² ∈ [9.99, 14.7]. Float scans: [dxtest.py](dxtest.py), [dxscan.out](dxscan.out).
- Verify: `uv run --no-project --with python-flint python3 dxcert.py 0.49 304.60 305.00` (14 s, 35 MB RSS, under memcap). A reviewer should check the three-line Fourier derivation and the Wiener step (Φ_{2,y} ∈ L¹: Φ decays like exp(−πe^{2|t|})).
- Novelty: Koide covers Thm 3.2 and gives corrected criteria; the Φ-instance (Thm 1.1 ⟺ ¬RH) is not in his note. Still to check: papers citing 1606.05011 that use Thm 1.1 (a scout, via arXiv listings or Google Scholar's cited-by; no contact).
- Existing scripts: CHAudit/zp_rs.py and DirLi/sanity.py use the same acb_series zeta jets; nothing in Riemann/Deep touches Dimitrov–Xu (rg).
- Smallest case: at y = 0.499 the first sign change in [14, 304] (step 0.1) is at x ≈ 110.45, next to the zero 111.0295. Certified by [dxcert2.out](dxcert2.out): R(110.40) > 0 > R(110.50). Small y needs larger heights (none at y ≤ 0.1 in [10⁴, 10⁴ + 4]); every fixed y ∈ (0, 1/2) is expected to fail high enough. Open: the infimum of the failing y.

### 2. Authenticate Khan–DeWeese κ = 0.6737113329 (N8; scout or one Orca agent, about 30 min)
- Goal: decide whether the N8 record line cites 24166562 as "author-kernel-checked, unreplayed" or waits for arena admission. Re-read www.riemannzeta.fun first.
- Inputs: Zenodo 24166562 (q9 2.3 MB, full workspace with hosted log; q10 14 MB; q11 55 MB; q12 243 MB); the arena's trusted comparator files (github josusanmartin/riemann: comparator/ChallengeDeps.lean, Challenge/Candidate.lean).
- Method, steps (1)–(4) done here (LOG.md, 19:00 UTC):
  - (1) ChallengeDeps.lean, ChallengeDeps/XiPrime.lean and Challenge/XiPrime.lean are byte-identical to zeta-23-lean 3635e748 (the arena's contract.json pin); CandidateSpec matches the arena template up to its two numerals.
  - (2) Q12/Final.lean l.2814–2826 states the three theorems verbatim as challenge/templates/Challenge.Candidate.lean.
  - (3) rg finds sorry, axiom and native_decide only in comments.
  - (4) κ is exact: 0.6737113329114, and the strict improvement is against the spec's 0.6736011515.
  - Left: (5) the EXPORT.sha256 cross-check; (6) Q12/Simple.lean's two simple-zero statements against the Multiplicity topic; (7) re-read the arena ledger. A kernel replay needs about 14 GB, far above our cap.
- Certificate: none (authentication only). Verify: sha256 and rg only, each process < 300 MB.
- Existing scripts: none. R451 and R457 authenticated arena entries through its site and GitHub (Riemann/Tree/R451.lean, R457.lean).

### 3. ξ′ record line: zeta23 0.86864 and Iizumi 0.87328 (N8 sibling; scout)
- Goal: correct the N8 line (R392's "0.8137 (Conrey) unchanged"). Cite zeta23 (Lean, comparator, 10 Aug 2026): 0.85838 / 0.86864 for simple zeros on the line, 0.92919 / 0.93432 distinct. Add Iizumi 24210376 (> 0.87328 simple on the line, unrefereed, author's build) once audited.
- Method, steps (1), (3) and (4) done here (LOG.md, 19:05 UTC):
  - (1) The statement uses the library's Zeta23.XiPrime counts at 3635e748, with definitions matching the trusted comparator file.
  - (3) His compiler log prints [propext, Classical.choice, Quot.sound] for all four final theorems; native_decide occurs only in comments of unused Distinct839/Local files.
  - (4) Clean-run JSON: success; 2553 fresh modules; Mathlib not rebuilt from source.
  - (2) `python3 tools/verify.py`: PASS (1775 byte-identical proof modules, 2,147,504 leaves, 6315 audits, four standard-axiom final audits; it checks provenance, not a fresh compilation; 126 MB RSS).
  - Left: the policy decision (cite as "statement-audited, author-kernel-checked").
- Also verify the attribution of 0.8137: is it Conrey 1983, or a later mollifier paper? Our files say Conrey.

## Routine one-liners (for COROLLARIES; side results, not on a node)
- arXiv 2611.14579 l.257 (random Chowla and Goldbach): "for any even N > 2 there is n < N with λ(n) = λ(N−n) = −1 … still open". It holds:
  - when λ(N/2) = −1 (take n = N/2);
  - when L(N−1) ≤ 0: the (N−1−L(N−1))/2 ≥ N/2 values with λ = −1 cannot all sit in distinct pairs {n, N−n}, of which there are (N−2)/2;
  - for every N ≤ 4·10¹⁸ (Goldbach verified, Oliveira e Silva–Herzog–Pardi 2014);
  - for almost all N: Σ_N |Σ_n λ(n)λ(N−n)|² ≤ sup_α|Σλ(n)e(nα)|²·x = o(x³) by Davenport's uniform bound.

  Every single N stays open. Not RH.
- arXiv 2611.16062 l.351, 599–624: the Siegel-zero term E_Siegel uses an ineffective β. The exceptional χ is quadratic, so L(s,χ) = ζ_{K′}/ζ_K, and Cureton–Thorner's Lean-checked strip (every Dedekind zeta, Re s > 31/32) gives β ≤ 31/32, hence E_Siegel ≪ X^{31/32+ε}. Their prime-ideal asymptotic is then effective (its other inputs are quantitative). N11-type, not node-central.
- Zenodo 24186091 (Mittermeier, 13 Nov): "the completed L₂₃ satisfies every Turán inequality unconditionally". Plausibly routine (Csordas–Norfolk–Varga argument: η(iy)η(23iy) > 0 and the kernel is log-concave); not checked. N6 analogue for one degree-2 L-function.

## Dropped hits (scanned, nothing for our tools)
- arXiv 2611.16550 (CLT for Frobenius, Hecke and Kloosterman angles): the MRH of earlier work is bypassed (l.377); nothing conditional is left that our tools supply.
- arXiv 2611.16828 (short exponential sums with primes, major arcs): η_k explicit (l.120–133); x ≥ x₀ asymptotic, second moment of L only; not a node.
- arXiv 2611.16156 (odd zeta values): s₀(ε) ineffective (l.141–155); off-tree.
- arXiv 2406.12152 v3 (Stieltjes error to li): fully unconditional and explicit already (section1 l.42); nothing open.
- arXiv 2611.16358, 16011, 17236 (a conjecture on permutation-matrix characteristic polynomials, l.368, off-tree), 14655, 16997 (Brown–Shields disproved; Dirichlet-space cyclicity, not Nyman–Beurling), 16506, 15505, 16713, 15534, 17218, 15558, 16254, 15379, 17189, 17210, 15696, 16536, 16923, 16157: no node statement.
- Zenodo 24210937 (Mori, Erdős #462 under RH): needs RH for intervals of length x^{1/2}(log x)²; no half-plane or zero-density input reaches that; off-node.
- Zenodo 24152440 (kusasyu36 audit of OpenAI's October math): no error in 029's zero-free strip Re s > 1 − 10⁻⁶; no Comparator mismatch in 12 formalized results; consistent with NOGO l.217 and hpaudit.
- Zenodo 24200541 (code for 24150355, R443), 24192507 (Ursachi update: no new interval), 24169192 (Munawer: d₁₄, d₂₀ only), 24197411 (heuristic HA numbers), 24159888 (prime races mod 24 to 7·10¹²; off-tree record), 24155142 (R449); the rest are AI/crank 'proofs'.
- Zenodo 24085183, 24090490 (R456); 24215600, 24215798, 24215625 (R467): already checked.

## Flags for our own files
1. N8 record line (TREE l.12): Zenodo 24166562 claims κ = 0.6737113329 for distinct and for simple zeros, above our 0.6736392880. The kernels were run by the authors only; the arena is unchanged. Mini-brief 2.
2. N8 line / R392: "ξ′ proportion 0.8137 (Conrey) unchanged, no 2025–2026 improvement" is wrong. zeta23's XiPrime theorems (anthropics/zeta-23-lean c77974d8, 10 Aug 2026; Lean, comparator okay) give 0.85838 / 0.86864 for simple zeros on the line and 0.92919 / 0.93432 distinct; Iizumi (24210376, 14 Nov) claims 0.87328. Also check that 0.8137 is really Conrey's.
3. No other pre-emption or error noticed. CGLUnc's "Alpöge–Furman PairCeiling 0.6818287" matches zeta23's README l.50.

## Drafts for the orchestrator (after a reviewer verdict)
- REFUTATIONS: Dimitrov–Xu, "Wronskians of Fourier and Laplace transforms", TAMS 372 (2019) 4107–4125 (arXiv 1606.05011).
  - Thm 3.2 (off-axis Wiener criterion for the LP class) is false: Lemma 3.3 drops the −i of an odd-kernel transform (Koide, Zenodo 24215530, 14 Nov 2026, explicit counterexample).
  - Thm 1.1 (RH ⟺ translates of cosh(ty)ν₂(Φ;t) dense for all 0 < |y| < 1/2) is equivalent to ¬RH: F = 2Re(Ξ′² − ΞΞ″)(x+iy) changes sign at y = 0.499 between x = 110.40 and 110.50, and at y = 0.49 between 304.60 and 305.00 (Arb, harvest7 dxcert).
  - The correct kernel ∫(t−2s)²cosh(y(t−2s))Φ(t−s)Φ(s)ds gives Jensen's J ≥ 0.
- N8 line: "Unrefereed, author-kernel-checked: κ = 0.6737113329 for distinct and simple zeros on the line (Khan–DeWeese, Zenodo 24166562, 12 Nov 2026; Lean + nanoda + comparator replay on the authors' machine, statement audited by harvest7, not on the arena)". Replace "xi' proportion 0.8137 unchanged" by "zeros of ξ′: ≥ 0.86864 simple on the line, ≥ 0.93432 distinct (zeta23 XiPrime, Lean, comparator, 10 Aug 2026); > 0.87328 simple on the line claimed (Iizumi, Zenodo 24210376, 14 Nov 2026, Lean, author's build, statement audited)".
