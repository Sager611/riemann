REVIEW
UTC: 2026-11-22 23:25
sorry count: 0 (rg 'sorry|admit|native_decide|axiom|implemented_by|extern|unsafe': one hit, Main.lean:37 '#print axioms nb_rate', the check command itself)
main theorem: Riemann.Deep.LeanNBRate.nb_rate (Main.lean) : NBInputs → (hRH : RiemannHypothesis) →
  Summable (m_ρ²/|ρ|² over NTZ) ∧ L·D(λ)² → zeroSum (λ → 0⁺) ∧ (all m_ρ = 1 → zeroSum = 2 + γ − log 4π) ∧ L·D_c(λ)² → zeroSum
build: bash .scratch/lean4/leanrun.sh -- lake build Riemann.Deep.LeanNBRate.Main → "Build completed successfully (3768 jobs)", peak 3118 MB
#print axioms nb_rate: [propext, Classical.choice, Quot.sound]
literature inputs (structure NBInputs, Main.lean; statements and verbatim sources in LitInputs.lean, LitInputsU1.lean):
  Burnol2002_Thm (Adv. Math. 170 (2002), arXiv:math/0103058, Thm, TeX l.169–206); Nazarov1993_Turan (as in Friedland–Yomdin,
  Studia Math. 218 (2013), arXiv:1107.0039 Thm 1.1); GoldstonGonek2007_Thm1 (Bull. LMS 39 (2007), arXiv:math/0511092, RH);
  Titchmarsh_Thm9_2, Titchmarsh_Thm9_7, Titchmarsh_2_1_5 (Titchmarsh, 2nd ed. 1986); Davenport_Ch12 (GTM 74, ch. 12 (10)–(11));
  Rudin1987_Thm9_13 (Plancherel on L¹ ∩ L², Real and Complex Analysis 3rd ed.). (Rudin1987_Thm19_3, Paley–Wiener, is defined in
  LitInputs.lean but not used.)
files: Riemann/Deep/LeanNBRate/: Defs (statement objects), LitInputs, LitInputsU1; U2 (dual bound): ExpSum, ExpSumOps, Kernel,
  Multiplier, MLBound, Gevrey (port, Apache-2.0), Windows, Nikolskii, Bernstein, LemmaB, FourierDeriv, LemmaLambda, LemmaD,
  Groups, Config, HighR, Legendre, LowPart, LemmaE, Counting, U2parts, ZeroConfig, GroupSums, U2grp, U2aux, U2asym, U2,
  U2props, U2main, U2final; closure: Mellin, Closure, LaplaceClosure, Bridge, Assembly; U1 (construction): Dual1–5, GammaBd,
  ZetaBd, Window, Rect (part ported from PNT+, Apache-2.0), Strip, Contour1–8, HoloInt, KJ1–5, ContourBd, KBd, RtrBd, GF1,
  LineE, HConc1–3; final: Final1–5, Main. PLAN.md (lemma DAG); audit/ (22 Nov memory probes, superseded).
proof outline: lower bound from Burnol (Bridge); D ≤ D_c; upper bound for D_c (Final5.upper_bound): for large L, the
  representer p* of the zero configuration up to a Titchmarsh-9.7 height T₂ ∈ [L⁶, L⁶+1] has N(p*) ≤ (S+ε)/ℓ (U2,
  ℓ = L − √L); H = (1 − E_low)Φ is entire, vanishes to order m_ρ at the low zeros (HConc1–2); g = inverse Laplace of H/ζ;
  contour shifts give ζG = H − K − ζR_tr on the line (KJ5); closure (LaplaceClosure) with c₀ = G(1); line bounds:
  ∫|E_low/s|² = 2πN(p*) (LineE, Plancherel), 1 − Φ, K, R_tr, G(1) tiny (KBd, RtrBd, GF1, Final1–4).
NOT-FORMALIZED: the decimal "= 0.0461…" (numerical value of 2 + γ − log 4π; needs a 5-digit enclosure of γ, absent from Mathlib)
last reviewer verdict: none yet (round 1)
next step: reviewer rounds (statement faithfulness + vacuity; literature inputs), then DONE

