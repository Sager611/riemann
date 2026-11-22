WORKING
UTC: 2026-11-22 17:09
sorry count: 0 (all committed modules sorry-free; remaining obligations are unwritten modules, see PLAN.md)
files: Riemann/Deep/LeanNBRate/: Defs (statement), LitInputs, ExpSum, ExpSumOps, Kernel, Multiplier, MLBound,
  Gevrey (port, Apache-2.0), Windows, Nikolskii, Bernstein, LemmaB, FourierDeriv, LemmaLambda, LemmaD, Groups, Config,
  HighR, Legendre, LowPart, LemmaE, Counting, U2parts, ZeroConfig, GroupSums, U2grp, U2aux, U2asym, U2, U2props,
  U2main, U2final; PLAN.md (lemma DAG); audit/ (22 Nov memory probes, superseded)
main theorem: (planned) Riemann.Deep.LeanNBRate.nb_rate : LitInputs → MainClaim (Defs.lean)
literature inputs (LitInputs.lean, def … : Prop bundled in structure LitInputs): Burnol2002_Thm (math/0103058);
  Nazarov1993_Turan (Friedland–Yomdin 1107.0039 Thm 1.1); GoldstonGonek2007_Thm1 (math/0511092, RH);
  Titchmarsh_Thm9_2, Titchmarsh_Thm9_7, Titchmarsh_2_1_5; Rudin1987_Thm19_3; Davenport_Ch12
proved (ours, sorry-free): Turán for jets; Nikolskii/Bernstein/Markov for exponential sums; [N] Lemma B, Lemma Λ, Lemma D;
  grouping at scale s0; Theorem R'' high part; Christoffel bound (Legendre); [R] Claim F core; Lemma E'; counting and
  RH cluster bounds; summability of Σ m²/|ρ|²; **U2 = Theorem R' + Corollary** (U2final.u2: RH ⇒ ∀ε>0, ∀ℓ large, every
  zero configuration of height ≤ ℓ⁷ has |p(0)|² ≤ (S+ε)/ℓ·Q_ℓ(p))
NOT-FORMALIZED: none decided
last reviewer verdict: none
next step: Dual (representer p*, E_low, OU kernel identity), then U1 (window Φ, rectangles/Cauchy, K_n/J_n, trivial
  zeros, entire G, growth, Paley–Wiener, closure), lower-bound bridge, value, Main
