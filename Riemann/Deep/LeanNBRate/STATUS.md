WORKING
UTC: 2026-11-22 18:10
sorry count: 0 (all committed modules sorry-free; remaining obligations are unwritten modules, see PLAN.md)
files: Riemann/Deep/LeanNBRate/: Defs (statement), LitInputs, ExpSum, ExpSumOps, Kernel, Multiplier, MLBound,
  Gevrey (port, Apache-2.0), Windows, Nikolskii, Bernstein, LemmaB, FourierDeriv, LemmaLambda, LemmaD, Groups, Config,
  HighR, Legendre, LowPart, LemmaE, Counting, U2parts, ZeroConfig, GroupSums, U2grp, U2aux, U2asym, U2, U2props,
  U2main, U2final, LitInputsU1, Mellin, Closure, LaplaceClosure; PLAN.md (lemma DAG); audit/ (22 Nov memory probes, superseded)
main theorem: (planned) Riemann.Deep.LeanNBRate.nb_rate : LitInputs → MainClaim (Defs.lean)
literature inputs (LitInputs.lean, def … : Prop bundled in structure LitInputs): Burnol2002_Thm (math/0103058);
  Nazarov1993_Turan (Friedland–Yomdin 1107.0039 Thm 1.1); GoldstonGonek2007_Thm1 (math/0511092, RH);
  Titchmarsh_Thm9_2, Titchmarsh_Thm9_7, Titchmarsh_2_1_5; Rudin1987_Thm19_3; Davenport_Ch12;
  Rudin1987_Thm9_13 (Plancherel on L1 ∩ L2, LitInputsU1.lean)
proved (ours, sorry-free): Turán for jets; Nikolskii/Bernstein/Markov for exponential sums; [N] Lemma B, Lemma Λ, Lemma D;
  grouping at scale s0; Theorem R'' high part; Christoffel bound (Legendre); [R] Claim F core; Lemma E'; counting and
  RH cluster bounds; summability of Σ m²/|ρ|²; **U2 = Theorem R' + Corollary** (U2final.u2: RH ⇒ ∀ε>0, ∀ℓ large, every
  zero configuration of height ≤ ℓ⁷ has |p(0)|² ≤ (S+ε)/ℓ·Q_ℓ(p)); Mellin transform of {θ/x};
  Mellin–Plancherel; closure: nymanDist(e^{−L})² ≤ (2π)⁻¹∫|(1+ζ(s)(c₀+∫_0^L g e^{−us}))/s|² for continuous g with
  c₀+∫g e^{−u} = 0 (LaplaceClosure.nymanDist_sq_le_laplace)
NOT-FORMALIZED: none decided
last reviewer verdict: none
next step: Dual (representer p*, E_low, OU kernel identity), window Φ, contours/Cauchy, K/J, trivial zeros, entire G,
  growth, Paley–Wiener (PLAN §3), lower-bound bridge, value, Main
