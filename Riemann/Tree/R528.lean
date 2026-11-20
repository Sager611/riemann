/-!
# Round 528 — node N1 (Weil/Connes positivity), round 17, [check] of Ursachi's fifth version (Zenodo 24400244, 20 Nov 2026 15:10 UTC): ten new statements (§§55–61, short and long excursions at ordinary primes; Thm 60.1 a cofinal RH-equivalent) hold as checked; no RH claim; three explorations, no transfer

* Node: N1 (17th round; R518 before). Round type: [check] of a new version of a paper on the node line (o44_zsweep.py 15:41 UTC: 24400244 = v5 of concept 23401731, file CP8_26_Short_and_Long_Prime_Excursions.pdf, 15,050 pdftotext lines against 13,914 for v4 24395180).
* Orchestrator diff (statement labels, lesson 343): 73 labels against 63, none removed; new: Lemma 55.1, Lemma 56.1, Prop. 56.2, Prop. 57.1, Lemma 57.2, Cor. 57.3, Prop. 58.1, Prop. 59.1, Thm 60.1, Prop. 61.1 (§§55–61, v5 l. 13380–14536), split into five chunks for the scouts.
* (1) c1: Lemma 55.1 (an exact normalisation of window coefficients) HOLDS; numerical check at D = 10⁴: equality after the fractional final coefficient.
* (2) c2: Lemma 56.1 (a cofactor mask with two exact Möbius moments) and Prop. 56.2 (native-tent truncation outside O(X^{1−5η+ε}), D = ⌊X^{4/5−η}⌋) HOLD; finite check at X = 10⁶.
* (3) c3: Prop. 57.1 (U_X and V_X have growth exponent 2Θ − 1), Lemma 57.2 (diagonal energy Σ d_j ~ (πa/32)c_g ℓ_X) and Cor. 57.3 (an RH-equivalent for the fully centred distinct-prime contribution; its arithmetic inequality is stated as open) HOLD.
* (4) c4: Prop. 58.1 (grid error o(Δ)) HOLDS. Its FLAW on Prop. 59.1 (V₀ ≪ H X^{3 − (14/13)(1−Θ)+o(1)} ≤ X^{901/260+o(1)}; 'the (U/K)^{15/13+ε} density growth is asserted, not derived') is wrong (orchestrator, c4.txt l. 136–146): it is Guth–Maynard's N(σ, U) ≪ U^{(30/13)(1−σ)+ε} with (30/13)(1 − σ) ≤ 15/13 for σ ≥ 1/2, uniform in σ, and the factor (K/U)² makes the dyadic sum converge; the exponent 901/260 = 3/5 + 3 − 7/52 uses the 7/8 theorem (Θ ≤ 7/8). Not counted.
* (5) c5: Thm 60.1 (with s = (14/13)(1 − Θ), sharp short-interval mean squares of ϑ(y) − ϑ(y − M) − M give a cofinal RH-equivalent; arithmetic 1 + 2·7/8 + 3/5 = 67/20 checked) and Prop. 61.1 (short-excursion prime-count bounds; 4052235189/4082225020 = 0.99265 < 139/140) HOLD.
* Explorations (3 scouts; one returned no answer): Foster–Vohra calibration (game theory, N7: calibration averages paired real parts to 1/2 and cannot place individual zeros; a symmetric off-line pair passes); Kronrod positivity and Stieltjes interlacing (numerical quadrature, N13: automatic for positive atomic zero data, no operator, no control of off-line zeros); Kramers–Kronig sum rules: no answer.
* Verdict: v5 adds correct statements (one more cofinal RH-equivalent) but no claim on RH and no node record.
-/

namespace Fayed.R528
end Fayed.R528
