WORKING
UTC: 2026-11-22 21:31
sorry count: 0 (all committed modules sorry-free; remaining obligations are unwritten modules, see PLAN.md)
files: Riemann/Deep/LeanNBRate/: Defs (statement), LitInputs, ExpSum, ExpSumOps, Kernel, Multiplier, MLBound,
  Gevrey (port, Apache-2.0), Windows, Nikolskii, Bernstein, LemmaB, FourierDeriv, LemmaLambda, LemmaD, Groups, Config,
  HighR, Legendre, LowPart, LemmaE, Counting, U2parts, ZeroConfig, GroupSums, U2grp, U2aux, U2asym, U2, U2props,
  U2main, U2final, LitInputsU1, Mellin, Closure, LaplaceClosure, Bridge, Assembly, Dual1–5, GammaBd, ZetaBd,
  Window, Rect (part ported from PNT+, Apache-2.0), Strip, Contour1–8, HoloInt, KJ1–5;
  PLAN.md (lemma DAG); audit/ (22 Nov memory probes, superseded)
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
  c₀+∫g e^{−u} = 0 (LaplaceClosure.nymanDist_sq_le_laplace); Bridge: Burnol lower bound, D ≤ D_c, value;
  Assembly.mainClaim_of_upper: MainClaim from Burnol + T9.2 + Davenport + UpperBoundNyman (the only open part);
  Dual1–5: OU form as inner product, Riesz representer p*, N(p*) ≤ (S+ε)/ℓ from U2, |E_low| ≤ √N(e_s)√N(p*),
  Green form, E_low entire; GammaBd: |Γ(z)| ≤ √(2π)|z+1|e^{−π|Im z|/2} on ½ ≤ Re ≤ 3/2 (Phragmén–Lindelöf);
  ZetaBd: trivial zeros, ζ′(−2k), |1/ζ| left/right, |ζ| from (2.1.5), ζ′ Cauchy bound; Window Φ = M_a(s−½)^n;
  Rect/Strip: residue theorem on rectangles (finitely many simple poles), contour shifts across strips,
  half-strip Cauchy; Contour1–3: g = inverse Laplace of H/ζ on Re = 3, G = ∫_0^L g e^{−us} (Fubini),
  Cauchy formula on the right half-plane; Contour4–8: strips, trivial-zero residues, far-left limit,
  middle rectangle ⇒ G = F − R_tr − J on Re s > 3 (generic in H with HBound/HZeros);
  HoloInt: holomorphic parametric integrals; KJ1–5: K (divided differences) analytic off 1, K = ζJ on Re s > 3,
  R_tr and G analytic, identity theorem ⇒ ζG = H − K − ζR_tr on −1 < Re s < 1 (the critical line)
NOT-FORMALIZED: none decided
last reviewer verdict: none
next step: concrete H = (1−E_low)Φ (HBound; HZeros via interpolation/vanishing order), OU identity
  ‖E_low‖²_w = N(p*), zero configuration up to T₂, line bounds (1−Φ, K, ζR_tr, G(1)), assembly of UpperBoundNyman, Main
