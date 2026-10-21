NO DISCOVERY: our ring class density theorem plus Linnik's machinery gives no explicit unconditional p(Q) ≪ |D|^A below the record 455; the exponent is set by an explicit count of zeros within about 0.6/log|D| of s = 1, and the density theorem says nothing there.

# Least prime represented by a form: report (2026-10-20)
Proofs and audit: [NOTE.md](NOTE.md). Brief: `.scratch/orca/lpqf-brief.md`.

## Main question (all D): no
- Record: |D|^455 (Zaman, PhD thesis 2017, Thm 7.4.1); |D|^694 (Thorner–Zaman 2017). The arXiv and GXL sweep found nothing newer. Thorner–Zhang (arXiv 2412.01802) improve Weiss only group-theoretically, with unspecified exponents.
- Why 455: the proof needs e^{−Bλ}N(λ) < 1, where N(λ) counts characters with a zero within λ/log|D| of 1. Positivity gives N(0.569) ≤ 3365. Past 0.569 only the power-sum (Turán) bound exists, log N(0.6) ≤ 258. So B ≥ 258/0.569 = 453.
- Our density theorem at σ = 1 − λ/log|D| gives C(ε)e^{2λ}|D|^ε, which is useless for fixed λ. It is not log-free, and the |D|^ε loss is built into its proof.
- A log-free version with an unspecified constant (Kowalski–Michel type) would not help. K zeros at λ force A ≥ log K/λ, so the count has to be explicit.
- Our sharper large sieve does not move Zaman's bounds:
  - His power-sum constants do not involve the sieve.
  - The sieve's length condition (Y ≥ 2) is already met (Y ≈ 5.5).
  - Positivity is capped by the conductor term, not by the sieve.
  - The efficient mollifier method would need Möbius cancellation in Q(√D) below √|D|, which is not known.

## Side result: a new candidate (needs the reviewer)
- **Theorem B**: for every A > A* = 1.72419, where A* is the root of A³ − 3A² − 6.5A + 15 = 0 in (1,2), all fundamental D < 0 outside a density-zero set have the following property: every class contains ≫ x/(h log x) degree-one primes in (x/2, x] for every x ≥ |D|^A. So every form of such a D represents a prime ≤ |D|^{1.7242}.
- Previous record for density-one D: |D|^{20/3+ε} (Ditchen 2013, D ≢ 0 mod 8). Gaudet (2024) still cites it as the record. Ditchen's conditional bound (assuming Lindelöf) was 2+ε.
- How it works:
  - Our density theorem controls every zero left of A/2, for each single D.
  - Kowalski–Michel's log-free estimate for all weight-one theta series of level ≤ X (exponent 13/2), with Markov's inequality over D, controls the zeros right of A/2.
  - Genus characters reduce to quadratic Dirichlet L-functions; there, Kowalski–Michel, Landau–Page and the classical zero-free region exclude only o(X) discriminants.
- Theorem B needs our density theorem. With Pascadi–Thorner's density in its place the same argument gives 2.32; with Kowalski–Michel's estimates alone it gives 9. Neither is stated anywhere.
- Status: the proof is complete modulo our unrefereed density theorem and published inputs. Not refereed, no Lean. Novelty: not found (about 25 arXiv queries, 2 GXL claims; NOTE §3).

## Leads
- Any explicit count N(λ) on [0.57, 1] better than e^{258} lowers 455 directly, since A ≈ log N(0.6)/0.569.
- For fundamental D, Zaman's positivity step bounds a conductor ratio by 4/3; for these D the ratio is 1. Redoing his Table 4.14 could push positivity past 0.569 and lower 455. This is bookkeeping; not attempted.
- A better log-free density over D than 13/2 lowers A* (exponent 4 would give 1.60).
- Theorem B for non-fundamental D, and an explicit size for the exceptional set.
