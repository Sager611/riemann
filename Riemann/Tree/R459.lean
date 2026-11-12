/-!
# Round 459 — node N8 (critical-line proportion), round 26, the unseen Zenodo N8 records of R458 checked: Devine's 0.6792 simple and on the line under the vanishing-box hypothesis (Zenodo 21879591) holds as stated (conditional; was 2/3, Goldston–Suriajaya); Santibañez-Leal's positive density of simple critical zeros in windows [T, T + T^θ] from θ = 0.527 (Zenodo 23134787 v0.03; was 0.552, Steuding) has its key mean-square lemma proved (unrefereed; uses Wang's preprint); Yang–Yang's 'almost all' (22065921) exposes no file; Veltri's equivalence is Speiser's theorem; the AI monthly reports hold nothing unseen; the rest are cranks; explore44's concentration barrier goes to NOGO; learning theory, kinetic theory and Morse inequalities transfer no tool

* Node: N8 (26th round; R458 before). Round type: [check] of the unseen Zenodo records next to the node (R458's list of 20; lesson 216),
  plus the result of the exploration agent explore44 and three exploration scouts.
* [check] Devine, Zenodo 21879591 (10 Aug 2026), 'A 67.92% Bound for Simple Critical Zeros in a Vanishing Vertical Box': (27) under (1)
  max |β − 1/2| log T → 0 over T < γ ≤ 2T, liminf N0,s/N ≥ 0.6792; the box perturbation (21)–(26) is o(N) uniformly (endpoint vanishing,
  two integrations by parts, the close-pair bound of Goldston–Suriajaya 2603.28104, who proved 2/3 under (1)); the CGdL kernel and BGSTB's
  complex-difference theorem are both used on |α| ≤ 1: HOLDS (a conditional result: CGdL's RH value 0.6792 transferred to the weaker
  box hypothesis; node line).
* [check] Santibañez-Leal, Zenodo 23134787 v0.03 (4 Oct 2026), 'Levinson's method in short intervals and simple zeros of the zeta function':
  Thm 1.1 (6): for every fixed θ ∈ [0.527, 1), liminf S(T, T^θ)/N(T, T^θ) > 0, with h(527/1000; κ) > 0.0005947542001 for a degree-201
  detector (ν = 499/10000), from a short-window mollified moment (Thm 1.2, 0 < ν < min(1/2, 2θ − 1, (5θ − 2)/6)), Wang's short-interval
  pair-correlation preprint (positive for θ > 0.550193964744) and a finite Hilbert–parity inequality; Lemma 3.1 (grouped frequencies
  log(pq/n), mean square) is proved, by a Gaussian Schur bound and by direct interval Gram arguments; the paper reports Steuding's θ ≥ 0.552
  as the previous range: a third-party claim of a wider range (unrefereed): node line.
* [check] Yang–Yang, Zenodo 22065921 (23 Aug 2026), 'Almost all zeros ... are simple and on the critical line' (lim N0s/N = 1 'subject to
  the verification grading'; trace moments of a Gabor-compressed Weil matrix through order 14): no file exposed (API files [], 403):
  NOT TESTABLE; the same moment-tower method as their 79.62% (21975237), whose order-4..6 moments came from the sine/GUE model (REFUTATIONS, R458).
* [check] Veltri, Zenodo 23380594 (14 Oct 2026), 'Charges and currents on the zeros of Riemann Zeta': Prop speiser (TeX l.272) 'RH holds
  iff B_ren has no null point in 0 < σ < 1/2' restates Speiser's theorem (B_ren = i·conj(ζ'/ζ)); l.126: 2λ₁ = Σ m_ρ|ρ|⁻²;
  λ₁ = 0.0230957089661210 (mpmath): KNOWN.
* [harvest] Van Vaerenbergh's 'AI and Mathematics Research: Monthly Report' (Zenodo 23797298 and the July–October issues): every item on
  zeta, RH or L-functions is already in our files (2608.13637, 2609.02882, 2610.12234, 2610.15559, 2610.18491, 2610.29636, 2610.15483,
  2610.21099, 2610.00311).
* Dropped after one look (no checkable claim on a node): 22210004 and 20985999 (Hamzah), 22181879 (Nogueira Grossi), 21989429 (Tully),
  21699477, 21669970, 21600554, 19446164 ('Fairy Monk'), 20033800 (Bonelli), 18924002 (Morató de Dalmases), 17688310 (Perisic), 17335645
  (Kim). With R458's three checks (22066689, 22287432, 21975237) this closes R458's list of 20.
* Agent explore44 (Orca; log-Sobolev, Herbst, T₂, Weissler–Bayart → N9/N7, N8): BARRIER scoped to five transfers (Herbst lower tail;
  Prop TL, a transport form of Littlewood's lemma: zero sum ≤ e^R·W₁(μ_T, μ_∞) + A(R) + O(log T/T); Herbst for Jensen; the Weissler–Bayart
  shift ½log₂(q/2); hypercontractive Lemma M): each conversion from the Bohr torus to [T, 2T] needs a W₁ rate, mean value or large-values
  input of density or moment strength; T^∞ bounds are rotation-invariant while [T, 2T] moments are not (c7 rerun by the orchestrator:
  identical, 4th-moment ratios 0.67 and 12.01 at θ = 0.91, 149 MB); Prop TL is a standard re-encoding: NOGO, no DISCOVERY.
* Exploration (each field had been tried: learning theory NOGO, kinetic theory R369/R389/R421, Morse flux R8.6): Sauer–Shelah (N5:
  threshold classes on zero ordinates have VC dimension 1, no bound), Cercignani's entropy-production coercivity (N1: reproduces the Weil
  certificates, c = 9.99·10⁻²² on [−7/8, 7/8]), weak Morse inequalities (N1: critical points of |ζ(1/2+it)|² ≥ simple critical zeros, no
  positivity): no transfer.
* Also seen: QuasiRiemannTracker PR 137 (TheSil, 12 Nov 16:05 UTC, open): a Lean proof that every Dirichlet L-function (ζ included at
  q = 1) has no zero with Re s > 2183367/2500000 = 0.8733468 (was B_2 = 0.87338196661320752): an N7 [check] round once admitted.
-/

namespace Fayed.R459
end Fayed.R459
