/-!
# Round 462 — node N9 (moments, Lindelöf, subconvexity), round 14, audit of two fresh unrefereed explicit S(T) bounds (Orca agent stcaudit, 12–14 Nov 2026): Casañá-Giner's |S(T)| ≤ 0.08382 log T + 4.444 and |N(T) − M(T)| ≤ 0.08382 log T + 5.326 for T ≥ 1 PASSES modulo its unrerun Arb certificates; Ismail's nine rows FAIL through a slip in Hasanalizade–Shen–Wong 2022 (log(1+η) for log ζ(1+η)) that Bellotti–Wong 2025 inherit, corrected rows computed; HSW's printed 9.3675, used by our records, is a misprint for 9.4925 (Robin2Adic unchanged on rerun); Schmidt subspace theorem, Bellman functions and moment problems transfer no tool

* Node: N9 (14th round; R438 before). Round type: audit of load-bearing explicit inputs (AGENTS.md priority): Orca agent stcaudit,
  21:27 UTC 12 Nov to 18:10 UTC 14 Nov (the Mac slept 22:31 UTC 12 Nov – 17:12 UTC 14 Nov), line by line plus falsification; texts in
  .scratch/r462/txt (deleted). Orchestrator 37's reading: Casañá-Giner (CG), Zenodo 23866074 (latest 1 Nov 2026), Thm 2.1, beats HSW's
  0.1038 log T; Ismail (IS), Zenodo 23301655 (11 Oct 2026), nine S(T) rows by the Backlund–Bellotti–Wong method.
* CG Thm 2.1, six lines of (2.1) |N − M| and (2.2) |S|, T ≥ 1 (Turing-type identity with LP-optimized kernels): AUDIT PASSED modulo the
  unrerun Arb certificates (Cζ ball sum and sup of the elementary part for all six lines; the kernel inequality (K) for ℓ150, ℓ400, ℓ6000).
  Proof re-derived (Lemmas 3.1, 3.2, Prop 4.1, (4.7)–(4.9), Lemmas 4.2–4.10); (K) certified by stcaudit's own interval code for ℓ5, ℓ20
  (smallest margin) and ℓ60; Cζ recomputed (5.76·10⁶ primes) ≥ 7·10⁻⁴ below the certified values; sup F sampled ≤ 1.6·10⁻⁶ below; all
  18 printed constants reproduce; end-to-end with Odlyzko's first 10⁵ zeros: identity residual ≤ 8·10⁻¹⁰, N − M inside the two-sided
  bound in 29/29 cases; exact S(T), N − M on [1, 74,920] (max |S| = 1.6287): every line holds, smallest slack 0.407. Thm 2.3 not audited.
* IS Thm 1.1, all nine rows, (1.3), Cor 1.4: FAILED as written. C1, C2 correct (own Bellotti–Wong transcription, ≤ 5·10⁻⁹; Ismail's
  69-check script passes); every C̃3 understated by the HSW slip below (+0.0275 … +0.0356 rows 1–6; +0.0376, +0.1401, +0.1336 rows 7–9);
  rows 7–9 also use Patel–Yang's unestablished 66.7 t^{27/164} (REFUTATIONS; 68.54 from their proof). Survives: |S(T)| ≤ 0.134791 log T
  + 0.159606 log log T + 1.534881 (T ≥ e); rows 1–6 C̃3 = 1.404841, 1.534881, 1.674012, 1.891768, 2.030092, 2.199195; rows 7–9 with
  68.54: 2.770550, 2.938158, 3.311666. Cor 1.3, Cor 1.4 and Prop 1.2 hold given the corrected rows.
* The slip (new; stcaudit's reviewer: CORRECT; the orchestrator read both TeX sources): Hasanalizade–Shen–Wong (arXiv 2107.06506) bound
  |f_N(s)| on −η ≤ σ ≤ 0 by interpolating ζ(1+η)(|Q + s|/2π)^{1/2+η} and (c₁/√(2π))|Q + s|^{1/2}(log|Q + s|)^{c₂}, then write the log
  of ζ(1+η)^{−σN/η} as −(σ/η) log(1+η); ζ(1+η) > 1/η > 1 + η and −σ/η ∈ [0, 1], so every C3 is understated: HSW headline +1.5·10⁻⁵,
  Table 2 rows +0.0078, +0.0090; Bellotti–Wong (Math. Comp. 2025, arXiv 2412.15470) inherit it (F_{c,r}, D3): headline 0.10076 log T +
  0.24460 log log T + 8.08344 needs 8.09025, other rows +0.0524, +0.0189, +1.1·10⁻⁵, +0.0593. No erratum found. → REFUTATIONS.
* Our records: HSW's printed C3 = 9.3675 is a misprint for 9.4925 (Bellotti–Wong's footnote, TeX l. 118); Robin2Adic, Ingham3, Ingham4,
  Li200k, RobinHeight and RobinCH use 9.3675. Orchestrator rerun of Robin2Adic/cert.py with 9.4926 (≥ 9.4925 + 1.5·10⁻⁵), python-flint,
  22 s, 38 MB: identical (K = 72 with KLN + RH(3·10¹²), 71 with δ ≥ 1, 75 with CS); the other reruns follow (records last).
* Exploration reruns (12 Nov ~22:00 UTC): Schmidt subspace theorem / Evertse–Schlickewei (needs algebraic inputs), Bellman functions (only
  trivial positivity of single Weil terms), Hausdorff/Stieltjes moment problems (= R296: Hankel positivity is automatic for positive atomic
  measures): no transfer.
-/

namespace Fayed.R462
end Fayed.R462
