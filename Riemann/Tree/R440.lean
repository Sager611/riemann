/-!
# Round 440 — node N7 (zero-free regions), round 13, audit of the strongest input of Liu's number-field B_2: Remark 12.18 of 2610.26835v2 (TeX l. 7836–7866) derives a Lindelöf-on-average fourth moment Σ_{q_k ≤ U} |X^{−1/2} Σ_n ψ_k(n) W(q_n/X)|⁴ ≪ U^{1+ε} at X = U^{1/2} for the sextic family over a cyclotomic F ⊃ μ12 (and up to X = U^{1/2+r_F/100}) from the z = 0 clause of Theorem M2, while for cubic characters the best published input (Blomer–Goldmakher–Louvel 2014, Thm 1.3, M + N + (MN)^{2/3}) gives only U^{4/3+ε}; no v3, no dependent paper; the clause becomes the next audit candidate; Kingman coalescent, portfolio theory and seismic inversion transfer no tool

* Node: N7 zero-free regions (13th round; R424 before). Round type: audit the load-bearing unrefereed input of a frontier claim (LiuAudit's open
  item G1, the z = 0 clause of Thm 12.16/M2, and its strongest consequence, Remark 12.18).
* [falsify] 2610.26835v2 (e-print number_fields_v4.tex), Remark 12.18, l. 7836–7866: for z = 0, ψ_k(n) = ν(n)χ_n(k) (sextic residue symbol,
  F ⊃ μ12 cyclotomic, exterior ideals prime to S), q_k ≤ U, ψ_k nonprincipal, X = U^{1/2}: Σ |X^{−1/2} Σ_n ψ_k(n) W(q_n/X)|⁴ ≪ U^{1+ε}; the same
  clause gives U^{1+ε}(1 + H)^A for X = U^m, m ≤ 1/2 + r_F/100 + ε, as used in Prop. refined-numerator; stated as a consequence of the M2
  argument (exact correlations, ordinary reflection, centering), not of the cubic large sieve; the cubic case 'requires a separate
  computation and is not decided here' (the sextic exceptional congruence 2 + v_p(h') ≡ 0 mod 6 forces v_p(h') ≥ 4, mod 3 it permits 1;
  residual power exponent 1/6 vs 1/3). Quoted exactly: the statement is as our audit cited; its proof is the unrefereed prose family construction.
* [check] cubic state of the art: Blomer–Goldmakher–Louvel 2014, Thm 1.3, M + N + (MN)^{2/3} gives U^{4/3+ε} at X = U^{1/2}: a
  Lindelöf-on-average fourth moment is open for cubic characters, so Remark 12.18 is a large-sieve-strength sextic estimate if correct.
* [check] OpenAI's Q(√−3) model, Lemma 18.1 (Lean-checked), is stated for 3/4 ≤ κ ≤ 1 with n1 + n2 + 6κz ≤ M when z > 0
  (Riemann/Deep/Kappa/NOTE.md l. 14); whether its z = 0 case contains the balanced fourth moment is not settled here.
* [check] versions: 2610.26835 is still v2 ('Comments are welcome'), Chang–Wang 2611.00924 v1; QuasiRiemannTracker PR 117 (0.9672) does not
  improve Liu, PR 114 pending; issues 128–136 report no problem with Liu's bounds. Dependents: 2611.00924 and 2610.24900 (both seen);
  nobody uses Remark 12.18.
* Exploration (tool transfers, tested): Kingman's coalescent and Ewens' formula need exchangeable sampled partitions, Markowitz/Kelly needs a
  covariance or return model (Kelly tried R143), least-squares moment-tensor inversion certifies only finite Gram positivity (statistical
  seismology tried R412): no transfer.
-/

namespace Fayed.R440
end Fayed.R440
