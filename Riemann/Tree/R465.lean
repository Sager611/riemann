/-!
# Round 465 — node N13 (Hilbert–Pólya and spectral approaches), round 14, the unseen Zenodo N13 records (no date filter, 98 records): Alarcón's claimed 6.1σ deviation of the gap-ratio limit below GUE is a fit artifact (1.1σ with log(T/2π)); Andrews's CCM reassessment withdraws only his own claims; Martini's CCM replication and Lee's localized Weil positivity at X = 7.01 are consistent and inside the records; Elder named and first measured the prolate ladder; the rest are frameworks or cranks; Gromov–Hausdorff spectral convergence, matrix completion and Bruhat–Tits buildings transfer no tool

* Node: N13 (14th round; R432 before). Round type: [check]/[falsify] of the unseen Zenodo records next to the node (lesson 227; Zenodo
  without a date filter on "Hilbert-Polya", "Hilbert–Pólya", "Berry-Keating", "self-adjoint operator" "Riemann", "spectral interpretation"
  "zeta", "Connes" "trace formula", "prolate" "Riemann zeros", "Connes-Consani": 98 records, 78 unseen); texts pre-extracted (lesson 228).
* [falsify] D. Alarcón, Zenodo 19297695 (28 Mar 2026), 'Gap ratio statistics of Riemann zeros': fits <r>(T) = 0.59891(13) + 1.245/log²T on
  21 points (Table 1, log T = 9.736 ... 24.145, Odlyzko and Platt zeros) and claims the limit lies 6.1σ below GUE (0.59971). Refit
  (orchestrator, weighted least squares on his table): a + b/log²T: 0.598908 ± 0.000132 (6.1σ, χ²/dof 0.497, reproduced); a + b/log²T +
  c/log³T: 0.598695 ± 0.000517 (2.0σ, 0.515); GUE fixed with b/log²T + c/log³T: χ²/dof 0.690; a + b/log²(T/2π) (the natural variable of
  the zero density): 0.599582 ± 0.000114 (1.1σ, 0.686). FLAW: the deviation is a model artifact, not a measurement.
* [check] R. Andrews Jr., Zenodo 23937987 v2.1 (4 Nov 2026), 'Convergence Claims for the CCM Construction: An Empirical Reassessment':
  withdraws his prolate/Lemma 7.2 falsification, the Śliwiński mean-error argument, a Mellin crossing comparison and a super-exponential
  claim (txt l. 42–45, 249–254, 280, 311); keeps 'accurate first roots do not give a uniform rate over an expanding window'. Our N13 line
  cites ACLG only for truncated-matrix numerics: HOLDS.
* [check] A. F. Martini, Zenodo 23393038 v4 (15 Oct 2026), 'Independent Replication of the CCM Spectral Construction': λ_min(c = 13) =
  5.999e-59, λ_min(c = 29) = 1.013e-97, zero-extraction errors 4.2e-55, 1.4e-93; his v1–v3 Theorem 5.3 withdrawn: HOLDS (consistent with
  our certified reports). Martini, Zenodo 23393801 v2: central conclusion withdrawn by the author.
* [check] B. Lee, Zenodo 22740489 v3.1r1 (14 Sep 2026), Thm 1.1: rigorous computer-assisted λ_{a_7} > 0 for the closed localization
  Q_W^{a_7}, a_7 = (1/2) log 7.01: inside Botha's L <= 2 and below Chekmachev's c <= 47: HOLDS, no record.
* [check] S. Elder, Zenodo 23113186 v1 (3 Oct 2026), 'The Prolate Ladder of Weil's Windowed Quadratic Form': exact low eigenvalues 1–5
  times one prolate leakage (n = 4, 8, 12, ... even; 6, 10, 14, ... odd); our ACLGFals Galerkin factors 1.67–12 are upper bounds,
  consistent; the name and first measurement are Elder's: node line now cites him.
* One-liners: Aksman 21813624, 21893282 (X0(163) and Fredholm–scattering 'Hilbert–Pólya' frameworks claiming RH), Warren 24068307 (claimed
  proof via Berry–Keating self-adjointness), Liu 23633249 (position bands from stated zero-free and proportion inputs): frameworks or cranks,
  no checkable new bound.
* Exploration: Gromov–Hausdorff / Kuwae–Shioya spectral convergence (N13): no metric measure space whose Laplacian carries Weil's signed
  prime kernel: no transfer. Matrix completion (N8): deterministic zero-Hankel samples have a rank-2 diagonal null space, so null-space/RIP
  recovery fails: no transfer. Bruhat–Tits buildings (N13): Ramanujan bounds for Hecke spectra have no trace map to zeta zeros: no transfer.
-/

namespace Fayed.R465
end Fayed.R465
