WORKING
UTC: 2026-11-22 15:06
sorry count: 0 (all committed modules sorry-free; remaining obligations are unwritten modules, see PLAN.md)
files: Riemann/Deep/LeanNBRate/: Defs (statement), LitInputs, ExpSum, ExpSumOps, Kernel, Multiplier, MLBound,
  Gevrey (port, Apache-2.0), Windows, Nikolskii, Bernstein, LemmaB, FourierDeriv, LemmaLambda, LemmaD, Groups, Config,
  HighR, Legendre, LowPart; PLAN.md (revised DAG); audit/ (22 Nov memory probes, superseded)
main theorem: (planned) Riemann.Deep.LeanNBRate.nb_rate : LitInputs → MainClaim (Defs.lean)
literature inputs (LitInputs.lean, def … : Prop bundled in structure LitInputs): Burnol2002_Thm (math/0103058);
  Nazarov1993_Turan (Friedland–Yomdin 1107.0039 Thm 1.1); GoldstonGonek2007_Thm1 (math/0511092, RH);
  Titchmarsh_Thm9_2, Titchmarsh_Thm9_7, Titchmarsh_2_1_5; Rudin1987_Thm19_3; Davenport_Ch12
proved (ours, sorry-free): Turán for jets; Nikolskii/Bernstein/Markov for exponential sums; [N] Lemma B, Lemma Λ, Lemma D
  (almost orthogonality, Plancherel); grouping at scale s0; Theorem R'' high part; Christoffel bound (Legendre);
  oscillatory IBP and cross-term bounds
NOT-FORMALIZED: none decided
last reviewer verdict: none
next step: LowPart Claim F core, Lemma E', U2 assembly, Dual (representer, OU kernel), then U1 (PLAN.md §3)
