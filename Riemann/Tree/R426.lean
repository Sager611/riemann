/-!
# Round 426 — node N11 (GRH families), round 15: harvest6's quasi-RH consumers checked: Gauss's class-number mean values (2611.10606) and the sphere problem (2610.29636) are unconditional by OpenAI's 7/8 (routine, COROLLARIES); super-positivity for every primitive real character of conductor up to 400,000 (Goldfeld–Huang criterion and Platt's verification, routine); 2611.09031's exponents do not depend on the strip; explicit short character sums from the half-plane still open; tribology, phase-field metallurgy and lithography give no route

* Node: N11 (fifteenth round; R414 before). Round type: harvest (harvest6's candidates 3-5 and its routine one-liners).
* 2611.10606 (Jia-Luo-Song-Xiong), Thm thm:cond (TeX l.105-112): only hypothesis 'no Dirichlet L-function has a zero in
  Re s > 7/8', used only in Lemma Lapprox for quadratic characters (l.810-846); l.105 cites OpenAI's announcement but keeps
  the theorem conditional: with OpenAI's Lean-verified 7/8 the error O(N^(2274811/3487744+eps)) is unconditional: HOLDS,
  routine (COROLLARIES).
* 2610.29636 (Bellotti-Kerr), Thm thm:qrh (l.68-75; hypothesis stated in the proof, l.1001-1004): sphere-problem exponent
  2274811/1743872 under Dirichlet zero-freeness in Re s > 7/8; not called unconditional by the paper: HOLDS, routine
  (COROLLARIES, joint entry).
* 2611.09031 (small defining polynomials II), Thms thm:qrh, thm:cubic (l.176-207): 'Only the implied constant depends on N';
  the exponents (below 1/2 without index-two subfields, (n-1)/(2n) primitive, 1/3 cubic) are the same for N = 32
  (Cureton-Thorner) and N = 8, so a wider strip (Liu's B_2 for Hecke L) gains nothing: HOLDS, no candidate.
* 2611.12152 (Sheth-Voo) l.242-255: super-positivity follows from Goldfeld-Huang 2017 Thm 1.2 if Lambda(s, chi) has no zero with
  1/2 < sigma < 1, |t| <= sigma - 1/2; the paper states positivity of theta(t, chi) only for 6,165,194 of 30,396,324 positive
  fundamental discriminants up to 1e8 (l.247) and no all-conductor range; Platt's verification (all primitive chi mod q <= 4e5,
  rigorous zero counts that exclude real zeros in (1/2, 1)) covers the region: super-positivity for every primitive real chi of
  conductor <= 400,000 (routine, COROLLARIES).
* Pre-emption check for the Orca candidate charsum78 (2611.08816 l.163: explicit constants for S_chi(N) << N^(7/8) q^eps 'remain to
  be seen'): arXiv and Zenodo give only GRH versions (2311.17231 Prop charsum) and long Polya-Vinogradov sums (Zenodo 23881837):
  still open: HOLDS.
* Exploration: tribology (Greenwood-Williamson moment positivity, rate-and-state Hopf stability), phase-field metallurgy
  (Mullins-Sekerka gradient flow with a positive Dirichlet-to-Neumann operator) and lithography (Hopkins' positive semidefinite
  transmission cross-coefficient and its coherent-mode decomposition): generic positivity, no prime structure: no route.
-/

namespace Fayed.R426
end Fayed.R426
