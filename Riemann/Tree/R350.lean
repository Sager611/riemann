/-!
# Round 350 — node N5 (Nyman-Beurling), round 7: the Moebius series converges in H^2(C_alpha) for every alpha > B (answers 2611.03924's question beyond alpha > 1; routine)

* Node: N5 Nyman-Beurling / Baez-Duarte (seventh round; R340 before).
* [check] Zenodo (since 28 Oct) and arXiv: no new d_N bound, criterion or claimed proof (2511.13496, 2507.13780 concern Beurling
  zeta functions and zero-free regions, not d_N): HOLDS.
* Harvest settled: 2611.03924 (TeX l.168-170, 199-206) asks for which alpha the Moebius series E = sum mu(k) G_k,
  G_k(s) = (k^{-s} - k^{-1}) zeta(s)/s, converges in H^2(C_alpha); it proves alpha > 1 and says alpha < 1 would need a zero-free
  half-plane beyond the known ones (it does not cite the 7/8 theorem). With the Lean-checked half-plane Re s > B,
  B = 0.87338196661320752: convergence for every alpha > B, so its threshold alpha_c lies in [Theta, B]. Proof (orchestrator):
  S_n - E = -(zeta(s)/s)(T_n(s) - T_n(1)), T_n(s) = sum_{k>n} mu(k) k^{-s}; Perron's formula shifted to Re(s+w) = B + delta gives
  T_n(s) << n^{B+2 delta-sigma}(|t|+2)^eps uniformly, and int |zeta(alpha+it)|^2 |t|^{2eps-2} dt < infinity. A first scout sketch
  via Abel summation alone loses a factor |t| (it gives only alpha > (1 + B)/2). COROLLARIES (routine).
* Computation: the paper's exact series for s_n(alpha) (Prop. exact) recomputed for n <= 200: log-slopes 0.436, 0.064, -0.119,
  -0.465, -0.616 at alpha = 0.3, 0.5, 0.6, 0.8, 0.9, consistent with its n^{1-2 alpha} envelope (TeX l.867-875); the law changes at
  1/2, not at B (KNOWN). An independent Gram-matrix recomputation of d_N timed out.
* Wild: a Nyman-Beurling criterion for number fields exists only for p = 2 (de Roton, JNT 129 (2009), Selberg class); the L^p
  form (1608.07887) is proved for Dirichlet L only, so density for p <= 32/31 for every Dedekind zeta (Cureton-Thorner's
  Lean-checked Re s > 31/32) needs a proof of the L^p criterion (open, small).
* Exploration: wavelets (dilation completeness is N5's own toolbox, R311, R340: KNOWN), Pisot/Salem numbers (Wiener-Tauberian
  reach only Re s = 1), Anderson localization (no Schrodinger cocycle realizes Xi; R61): no idea.
-/

namespace Fayed.R350
end Fayed.R350
