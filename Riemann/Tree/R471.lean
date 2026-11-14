/-!
# Round 471 — node N3 (Li/Keiper), round 13, the unseen Zenodo N3 records (no date filter, 98 records, 14 unseen on topic): two claimed RH proofs break on the real Li coefficients — Henderson's log-convexity λ_n² ≤ λ_{n−1}λ_{n+1} fails at every 2 ≤ n ≤ 78 (λ_n is strictly log-concave), and Robertson's coefficientwise bound for the zeta block of the Li series fails at k = 1 (γ − 1) and 15 ≤ k ≤ 30 (REFUTATIONS); Ueoka–Nagi–Sui's Toeplitz–Hankel criterion, Nader's finite-height Li bound and Schuster's dichotomy hold as stated; the N3 record stands; Patterson–Sullivan, Lempel–Ziv and Lindley transfer no tool

* Node: N3 (13th round; R451 before). Round type: [falsify]/[check] of the unseen Zenodo records next to the node (lessons 216, 227):
  Zenodo without a date filter on "Li coefficients", "Li criterion", "Li's criterion", "Keiper", "Keiper-Li", "Bombieri-Lagarias",
  "Li constants", "Stieltjes constants", "Maslanka": 98 records, 14 unseen on topic; 8 texts pre-extracted (lesson 228); the
  orchestrator computed λ_1..λ_79 to 30 digits (Cauchy integral of (ξ′/ξ)(1/(1 − x)) on |x| = 0.9, M = 1024, dps 50; λ_1..λ_5 =
  0.023095708966121, 0.0923457352280467, 0.207638920554325, 0.368790479492242, 0.575542714461177) for the scouts.
* [falsify] R. Henderson, Zenodo 16990875 v1.1 (28 Aug 2025), claimed RH proof via total positivity of the Keiper–Li coefficients:
  Proposition 2.13 (l.789–794) concludes λ_n > 0 and λ_n² ≤ λ_{n−1}λ_{n+1} from the total positivity of the Stirling–Lah matrix;
  the real λ_n are strictly log-concave at every 2 ≤ n ≤ 78 (n = 2: 0.0085277 > 0.0047956): COUNTEREXAMPLE → REFUTATIONS.
* [falsify] B. H. Robertson, Zenodo 17314347 v1.0 (10 Oct 2025), claimed RH proof via an 'HB/PSD–Li identification': Corollary 6.2
  needs 0 ≤ ℓ_k ≤ [z^k]L_ζ(z) for all k, L_ζ(z) = (s² − s)(ζ′/ζ)(s) at s = 1/(1 − z); orchestrator: [z¹] = γ − 1 = −0.422784,
  [z²] = −0.0331149, [z³..z¹⁴] > 0, [z^k] < 0 for 15 ≤ k ≤ 30 (the gamma block plus this block reproduce λ_1..λ_3); Theorem 4.1's
  identification for k ≤ 40 is a residual fit with the computed λ_k as input. Scout FLAW (unproved passage to all k) plus the
  orchestrator's sign counterexample → REFUTATIONS.
* [check] Y. Ueoka, A. Nagi, Sui, Zenodo 21762261 (2 Aug 2026): Theorem 5.2, RH ⟺ every leading block of c_ij = f_|i−j| − f_{i+j+1} +
  δ_ij f_0 is positive semidefinite, f_n = λ_{n+1} − 2λ_n + λ_{n−1} (f_0 = λ_1, f_1 = 0.0461543172958); the hard direction is proved
  in the text (Lemma 5.1, a Pringsheim radius bootstrap and the identity theorem, l.328–347). Orchestrator: smallest eigenvalues
  1.0365·10⁻¹¹, 6.403·10⁻¹⁹, 1.723·10⁻²⁶, 3.493·10⁻³⁴ for N = 2..5, then below the dps-50 noise floor (about 10⁻⁷·⁷ per size step,
  so the criterion is numerically ill-conditioned): HOLDS.
* [check] J. Nader, Zenodo 22851552 v1.0 (19 Sep 2026): corrected finite-height bound λ_n ≥ L_n(T) − (n²/4)e^{n²/(2T²)}B(T),
  B(T) = (log(T/2π) + 1/3)/(6πT³) + (2 log T + 129/16)/T⁴, tested for n ≤ 79 with the first 30 zeros (T = 101.3178510057, B =
  3.2298948509·10⁻⁷): no failure; two impossibility claims withdrawn by the author: HOLDS.
* [check] L. Schuster, Zenodo 23021017 (28 Sep 2026): Bombieri–Lagarias trend/oscillation dichotomy (Theorems 6.1–6.2, finite
  resolution detects an off-line zero or leaves an unresolved tail); λ̄_79 = 84.0465373947 against λ_79 = 84.4436939806: HOLDS.
* Not read beyond titles: Iciok 23044316 (a Taylor–Stieltjes certificate that ζ ≠ 0 for 0 < σ < 1, |t| ≤ 6.53: known, the first
  zero is at 14.13), Lange 22295615 (Toeplitz–Cayley rows), Balisa 23783653, Stieltjes-constant notes 23454569, 23461957, frameworks.
  N3 record (Fabbian: λ_n > 0 for n ≤ 6.93·10²⁶) stands.
* Exploration: Patterson–Sullivan (λ₀ = δ(1 − δ) transferred to Re ρ(1 − Re ρ) ≥ 1/4 is RH itself), Lempel–Ziv/Shannon–McMillan–
  Breiman (no canonical stationary source from zeros), Lindley/Kingman (no stochastic drift in the zero gaps; queueing tried R398,
  R440): no transfer.
* Scouts: 8, all answered; the two counterexamples recomputed by the orchestrator.
-/

namespace Fayed.R471
end Fayed.R471
