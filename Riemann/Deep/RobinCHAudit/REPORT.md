AUDIT PASSED: every step of RobinCH's proof (Lemmas 11_k, 12, 13 and the RobinHeight Lemmas 1–9 it uses) re-derived line by line, check_all2.py rerun with identical output, 30 certificate boxes recomputed independently at 256 bits, falsification hunts empty, reviewer CORRECT; no error, five harmless remarks (F1–F5); inputs inherited at statement level.

# Referee audit of RobinCH (Robin to 10^(10^46.79), Nicolas to p_k ≤ 1.78·10^48; node N4) — robinchaudit, 14 Nov 2026

Scope: Riemann/Deep/RobinCH (PROOF.md, certify*.py, check_all2.py) and the parts of Riemann/Deep/RobinHeight it uses,
line by line ([FULLPROOF.md](FULLPROOF.md), with source line numbers); Chirre–Helfgott 2512.15709v1 TeX (the only version);
Platt–Trudgian, Büthe, HSW TeX at statement level. 19:59–20:30 UTC. Nothing in RobinCH/ or RobinHeight/ edited.

## Verdict by record

| record (as written in TREE.md N4 / RobinCH REPORT) | status | rests on |
|---|---|---|
| Robin for 5041 ≤ n ≤ 10^(10^46.79) | holds | H1 RH to 3·10^12, H2 b0 = 0.87338196661320752, H3 HSW, H4/H5/H7/H8 classical, H6 Büthe, H9 Morrill–Platt, H10 CH Cor. psiexpl |
| Nicolas for p_k ≤ 1.78·10^48 | holds | same (no H9) |
| Robin 10^(10^46.44), no half-plane | holds | same without H2 (only 0 < β < 1 above T0) |
| Nicolas 6.30·10^47, no half-plane | holds | same without H2 |
| 7/8 rows (10^(10^46.79), 1.76·10^48), Level 2-CH rows | hold | rerun identical |

The chain uses RH only below T0 = 3·10^12 (and CH's proof only to 3·10^12 + 1 + π/3 < 3,000,175,332,800, PT TeX l.62).

## (1) Referee (details: FULLPROOF §§A–G)

- H10 = CH l.402–406 verbatim ("For any x ≥ 1, |ψ(x) − x| ≤ π/(3·10^12)·x + 113.67√x"); hypothesis RH to T (l.387–388,
  Prop. sagaro l.4868–4880); proof step l.4981–4991 rechecked; used only for x ≥ 10^19 e^(−10^−3). CH's analytic proof:
  inherited from CHAudit (§2.2, no error).
- Lemma 12: both directions of ψ → θ, the ψ − θ bound P, η_1 = η_CH(10^19) = 3.6263146·10^−8 on [10^19, 10^38], the (1 + η_1)
  factor beyond 10^38, monotonicity, and the window-uniform evaluation at L − λ/2: ✓.
- Lemma 11_k: decomposition, M_φ sandwich, W ≤ η H(L') E|V − c| (layer cake), MAD_k (exact), (B1'), (B2_k) (Fubini, two
  integrations by parts, MGF ((e^{sa}−1)/(sa))^k with |e^{sa} − 1| ≤ 2, H-variation term): ✓. Zero sums: low zeros with
  β_0 = 0.0461914 and S1(T0) = 119.988, high zeros via S_m(T0) and the pairing ρ ↔ 1 − ρ̄ (b0 = 1: factor (1 + 1/x')/2): ✓.
- Lemma 13: series bounds, |v| ≤ η(1+η)/L, every bracket term: ✓.
- Reductions (Lemmas 1–4, CA numbers, interval Lemma 3, Lemma 9 and its η_1 piece, Lemma 8's zones to 10^76): ✓.
- Code vs formulas (margins3/margins4, robin_pos, Mlow3/Mup3, box, run3, choose4): term by term identical (table in
  FULLPROOF §G); Lemmas 10/10_k not called (rg); λ ≤ 9·10^−4, k ≥ 1, exact binary λ.

## (2) Reruns and independent recomputation (one process at a time, memcap 290 MB, /usr/bin/time -l)

| run | result | max RSS |
|---|---|---|
| check_all2.py (unmodified copy, .scratch/robinchaudit/run) | ALL PASS; every line = RobinCH REPORT §4 (counts, minima, k ≤ 40/57, exponents 40.5656, 40.4896, 35.3789, 46.7931, 46.7900, 46.4480); 370 s | 87.6 MB |
| dump.py (RobinCH's own Level-3 search, every box recorded; 3,513 boxes) + cover.py | 0 frontier gaps; first box at 10^7 / 10^10; last box ≥ claimed end; probes fail at 10^48.252447, 10^47.156207, 10^47.799393, 10^46.810759 (= RobinCH's search ends ±1 in the 6th decimal) | 38 MB |
| indep.py (own code from the formulas, 256 bits, no RobinCH import), 30 boxes with RobinCH's (k, λ) | every RobinCH PASS passes, every probe FAIL fails | 42 MB |

Selected independent margins (normalised by √x log x; ours ≥ RobinCH's because monotone terms are taken at box ends):

| box (b0 = 0.87338…) | k, λ | ours | RobinCH |
|---|---|---|---|
| Nicolas, box containing 5·10^8 (10^8.6286–10^8.7372) | 16, 1.66·10^−11 | 1.6024 | 1.5727 |
| Nicolas, middle (10^27.9547–10^28.0633) | 39, 4.18·10^−11 | 1.6822 | 1.6590 |
| Nicolas, last box (10^48.24949–10^48.25042) | 29, 3.09·10^−11 | 0.0022612 | 0.0022593 |
| Robin, first box (10^10–10^10.0543) | 18, 1.51·10^−11 | 0.21524 | 0.21466 |
| Robin, middle (10^27.9689–10^28.0232) | 39, 4.19·10^−11 | 0.36764 | 0.36711 |
| Robin, smallest margin (10^47.14830–10^47.14914) | 30, 3.18·10^−11 | 2.04·10^−6 | 1.63·10^−6 |
| Robin, last box (10^47.15529–10^47.15534) | 30, 3.18·10^−11 | 4.497·10^−4 | 4.495·10^−4 |
| Robin, log x ∈ [20.03, 20.04] (OddRobin's use; own k, λ) | 16, 1.33·10^−11 | 0.468163 | 0.468 (oddraudit F4) |

At the Robin failure point (10^47.156207, b0 = 0.87338…) n(x_a)(Σp^−2 + L_2) = 2.7946 is used up exactly by M 1.9822, window 0.5351, D 0.2095, zeros 0.0678.
Constants recomputed: S1(T0) = 119.98765, S_2 = 2.95941·10^−12, S_3 = 4.84393·10^−25, β_0 = 0.0461914179, η_CH(10^40) = 1.04720·10^−12, C_T = 113.66888.

## (3) Falsification (FULLPROOF §I) — nothing found

- Lemma 11_k, 612 corner cases (k to 120, first zero, γ = 3·10^12 + 1 and 10^16, x = 10^7 and 1.78·10^48, β to 0.999, centred/one-sided, λ beyond 10^−3), by an MGF/Taylor method independent of test_b2k.py and cross-checked by quadrature: (B2_k) never exceeded, sharp to 1 − ratio = 6.7·10^−18 at λγ/k = π; (B1') worst 0.999995.
- Lemma 13: 144 cases, min(1 − D/bound) = 8.9·10^−15 > 0. MAD_k = quadrature to 20 digits (k = 3, 12, 40).
- Signs of all one-sided bounds, box coverage, no RH above T0: ✓.
- Endpoints: 10^(10^46.79) is below the certificate's first failure (X = 10^47.156207, i.e. 10^(10^46.7940)); the next boxes fail in both codes; a finer (k, λ) search also fails there (−6.55·10^−7): the 0.87338… records are this method's limit.

## Remarks (none affects a claim)

- F1. b0 = 1 endpoints are search-limited: at RobinCH's failure boxes, k = 55 (RobinCH tried 56/57) gives +0.0027 (Nicolas) and +0.0009 (Robin); every point box 0.002 dex further fails. Gain < 0.002 dex; claims conservative.
- F2. REPORT §1.5/§2 "λ T0 ≈ 90 (k ≈ 30–60)": λT0 ≈ 3.2k, i.e. 93–95 for the b0 < 1 rows (k = 29–30) but ≈ 182 for b0 = 1 (k = 55–56).
- F3. "404 boxes", "1161 boxes" count attempted boxes; the covering boxes are 389, 638 (the rest were halved and retried).
- F4. Platt–Trudgian's title and abstract say 3·10^12, which would not cover CH's T = 3·10^12 + 1 + π/3; their theorem (TeX l.62) states 3,000,175,332,800, as RobinCH and CH quote. Cite the theorem, not the title.
- F5. Ball evaluation of √x L·M(x) on a box loses up to 2.4% against endpoint evaluation (box widths to 0.25 in log x); only conservative.

## Reviewer

reviewer (verdict-only; Lemma 11_k (B2_k), the window bound, the high-zero pairing incl. b0 = 1, exact values pasted): CORRECT — Fubini, both integrations by parts, H'' = 2/u³ + 6/u⁴ > 0, the MGF bound, the layer cake and the pairing (fixed points on the critical line included) are valid; above T0 only the stated strip, the reflection symmetry and the S_m bounds are used, not RH.

## Inherited, not re-derived here

H2 (the 0.87338… / 7/8 half-planes; only the b0 < 1 rows), CH's analytic proof (CHAudit), H5 Ingham, H7, H9 Morrill–Platt, Robin 1984 / Nicolas 1983 (statement level via RobinHeight's quotes), python-flint/arb (E_1 = −Ei(−z), matches mpmath at two points).

## What survives

All four records as written. Weakest links: the preprint CH (above 10^19, every row) and, for the 46.79 / 1.78·10^48 rows only, the half-plane input.

Files: FULLPROOF.md; rerun.log, dump.log, cover.log, indep.log, fals.log, qcheck.log, lemma13.log, endfine.log, endreach.log; code: dump.py, cover.py, indep.py, fals.py (+ fals_lib.py), qcheck.py, lemma13.py, endfine.py, endreach.py (run dump/cover from a copy of RobinCH/; indep/fals/qcheck standalone; boxes.json regenerated by dump.py).

RH relevance: N4
