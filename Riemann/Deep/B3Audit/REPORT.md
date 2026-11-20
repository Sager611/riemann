AUDIT PASSED — every step of Theorem E1′ (Lemma B, Corollary B′, Lemma R and the residues, the §4 branches, the compact-window step and the θ-free chain) re-derived from SL's TeX and Wang's TeX, falsified numerically by an independent route to τ = T/H² = 251, Arb rerun identical and κ, k*, h recomputed independently to 30 digits; no error found.

# B3Audit — line-by-line referee of B3Cert Theorem E1′ (simple critical zeros in every (T, T + T^θ], θ > 371/750) — node N8
b3audit (Orca), 20 Nov 2026, 05:40–06:20 UTC. Object: Riemann/Deep/B3Cert/REPORT.md, FULLPROOF.md (B3F). Sources: SL v0.03 main.tex (Zenodo 23135354, the
line numbers B3Cert/E1Audit use), Wang arXiv 2609.07918 v2 TeX. Chain with line numbers: FULLPROOF.md (§ numbers below). No Lean, no edits outside B3Audit/.

## Verdict per item of the brief
| Item | Check | Result |
|---|---|---|
| (1) Lemma B (i) bulk | coordinates Λ = H(e^{2y/H} − 1), ξ = Tu/2π, joint phase K[ℓ − Λu + σc_r√u], Jacobian f = ½e^{−y²−(1+2β)y/H}, Step b cases, completed square, (c1)/(c2), ψ″ signs and c_r ≤ H/46 ⟸ 32√2 < 46, Titchmarsh 4.4 with boundary term, (c_r/K)^{1/2}X0^{−1/2} = K^{−1} (§1.1–1.5) | ✓ re-derived; reviewer subagent: CORRECT |
| (1) Lemma B (ii), Cor B′ | pointwise ∂^kA ≪ H^{−k} (|λ| ≤ (4π\|y\|/H)e^{2\|y\|/H}), IBP gain 1/c_r; n-split at N0T^ε, N0 ≤ y²T ≤ T², tail exponent, dual total y²T^{1/2+ε}H^{−1/2} (§1.6–1.7) | ✓ |
| (1) uniformity | constants depend on χ, C0, N only; no τ or θ in any step; shifts only in amplitudes (holomorphic, Cauchy for derivatives); centre T′ ≍ T | ✓ |
| (1) falsification | own direct quadrature (y- and δ-trapezoid, mpmath-validated to 10⁻¹³), 27 blocks, 216 rows, T = 10⁴–10⁶, θ = 0.30–0.45, τ = 2.5–251, three shift sets, r = 0.003–8; one term recomputed fully in mpmath (25 digits, = float to 8 digits) | max \|𝒯\|/(HX0^{−1/2}) = 0.668; each branch = e^{−πr/2}\|(T/2π)^{−d/2}\|/(2√2) (own 2D stationary phase, §1.8) to 10⁻³; no failing case |
| (2) Lemma R | v-substitution, Parseval, ∫ω̂e^{2πik} = ω(1) = 1, O(H/T) uniform on \|w − 1\| ≤ 1/4 (§2.1); F4: \|E\| ≈ 0.035/T ≪ H/T | ✓ |
| (2) residues | SL l. 441–447 with M_lead and outer factor 2π/(gq) reproduce SL's displayed R_lead (l. 713–719) exactly (own computation); pole cancellation, max modulus, o(H/L) (§2.2) | ✓ |
| (3) §4 branches | off-band (own derivative bound, tail e^{−H²/400}), SL (11) Taylor part absent, ε = +1 ≪ (H/T)^A for all y, K-branch with sup\|A\| only, E1 Steps 1–2 exact for any A ∈ C_c^∞ (E1.md l. 23–31) (§3) | ✓ no θ > 1/2, no L¹ gain, no Taylor step |
| (4) compact windows | SL l. 731–758: profile-free; "θ − η > 1/2" is SL's own Gaussian range, replaced by ν < (3(θ − η) − 1)/4 (§4.1) | ✓ |
| (4) chain | lem:count / Thm 1.3 (l. 823–925: RvM, Backlund, Littlewood, (Emean) l. 812–815), Lemma diagonal (ν < 1/2 only, rem:uses l. 336–345), Thm D (l. 932–954): parity product + Wang thm:local-pair "0 < λ < θ < 1" (TeX l. 232–240), error O(H + T^λL²) (§4.2–4.4) | ✓ θ-free on (0, 1) |
| (5) Arb | rerun of b3c_arb.py: output byte-identical; rerun of b3c_fals.py (10⁶: θ = 0.45, 0.40): 8 rows identical | ✓ |
| (5) independent κ, k*, h | mpmath 50 digits, trig-form Q, Gauss–Legendre, cross term integrated (Q(1) = 0 not assumed): κ(121/1000, 44/5) = 0.0850387603917657602445070238217, k*(371/750) = 0.08499141434975944309603663, h(371/750; κ) = 3.354573606042305490·10⁻⁵, θ_det = 0.494636567599178 | ✓ = Arb in every digit |

## Remarks (not errors)
- Lemma B is sharp termwise: the true size is HX0^{−1/2}e^{−πr/2}/(2√2) per branch, uniformly in τ; the gain over the trivial bound ‖A‖₁X0^{−1/2} ≍ √T X0^{−1/2} is τ^{−1/2}.
- Sanity: at ν → 0 the range θ > 1/3 + 4ν/3 tends to θ > 1/3, the range of the Gaussian-smoothed mean square of ζ from Atkinson's formula (smoothed error
  ≍ T^{1/2}H^{−1/2}, the same dual length N0 = T/H² and per-term size). E1-M′ is the twisted analogue; its exponent is the natural one.
- The hypothesis c_r ≤ H/46 of Lemma B(i) is a proof artefact (rows beyond it obey the same bound); Cor B′ only needs it for n ≤ N0T^ε.
- Inputs not re-audited here (audited earlier, unchanged for θ < 1/2): SL (10), (12), (14), E1 Steps 1, 2, 5 (E1Audit), Wang thm:local-pair (WangAudit),
  parity product (SantiAudit), Young's Lemma 7; for each, the θ-dependence was re-checked in the source (§3.5, §4).

## Files
FULLPROOF.md (§1 Lemma B, §2 Lemma R + residues, §3 branches, §4 windows + chain, §5 Arb, §6 numerics, §7 review); code/b3a_lemB.py (F3, F5),
b3a_lemR.py (F4), b3a_kappa.py, b3a_mpT.py (F6), run copies b3c_arb.py, b3c_fals.py; logs/ (memcap peaks ≤ 149 MB in *.err).

RH relevance: N8
