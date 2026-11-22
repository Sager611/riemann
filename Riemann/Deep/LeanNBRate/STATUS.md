WORKING
UTC: 2026-11-22 13:45
sorry count: 0 in finished modules (MLBound/… all sorry-free); open obligations listed in PLAN.md §2–4
files: Riemann/Deep/LeanNBRate/: Defs (statement), LitInputs, ExpSum, ExpSumOps, Kernel, Multiplier, MLBound,
  Gevrey (port of zeta23 Taper/Gevrey.lean, Apache-2.0), Windows, Nikolskii, Bernstein, LemmaB, FourierDeriv,
  LemmaLambda (part 1); PLAN.md; audit/ (memory probes of 22 Nov, superseded)
main theorem: (planned) Riemann.Deep.LeanNBRate.nb_rate : LitInputs → MainClaim (Defs.lean)
literature inputs (LitInputs.lean, `def … : Prop`, bundled in `structure LitInputs`): Burnol2002_Thm (math/0103058);
  Nazarov1993_Turan (as stated in Friedland–Yomdin 1107.0039 Thm 1.1); GoldstonGonek2007_Thm1 (math/0511092, RH);
  Titchmarsh_Thm9_2, Titchmarsh_Thm9_7, Titchmarsh_2_1_5; Rudin1987_Thm19_3 (Paley–Wiener); Davenport_Ch12
proved so far (ours, sorry-free): Turán–Nazarov extended to jets (forward differences); Fourier kernels and the
  multiplier identity/bound (replaces [N] Lemma A); Nikolskii; Bernstein on [0,ℓ]; [N] Lemma B for one group;
  Λ-representation F(ζ) = Λ_r(Θ_ζ)
NOT-FORMALIZED: none decided
last reviewer verdict: none
next step: Lemma Λ bound, Lemma D (Schur form), groups + RH cluster bound, Lemma E, Claim F, dual form, U2; then U1
note for the orchestrator: Gevrey.lean is proof code copied from anthropics/zeta-23-lean (Alpöge–Furman's repo,
  Apache-2.0): machine-checked, not an input; it proves only generic Gevrey bounds for exp(−1/x) bumps.
