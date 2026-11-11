RECORD: D(lambda) > 0 for every lambda >= 0.7851, unconditionally (was 0.7892, Lamzouri 2611.14486)

# SmallGap: small gaps between zeta zeros (node N8)

## Claim
A positive proportion of consecutive ordinates satisfy γ_{n+1} − γ_n ≤ 2πλ/log T for every λ ≥ 0.7851, unconditionally: liminf_T (1/N(T))·#{0 < γ_n ≤ T : γ_{n+1} − γ_n ≤ 2πλ/log T} > 0.

## Before
- 0.7892: Lamzouri 2611.14486 Thm Main (l.195-197), Selberg minorant λR^(λu) in his Thm Main2 (l.392-400) with θ = 4/7 − 10^{−5} (l.452-469).
- 0.7864: stated in his remark (l.199-201) with Chirre-Gonçalves-de Laat's degree-162 polynomial (optimised for θ = 1), no details.
- Earlier: 1 − 2.36·10^{−6} (Zhao 2603.17334), 1 − 10^{−3·10^{13}} (Simonič-Trudgian-Turnage-Butterbaugh 2010.10675).
- Under RH: λ > 1/2 (Durkan-Pearce-Crump), 0.6039 (BGMM). Not the target.

## Method
1. Reproduced Lamzouri's 1.000100147 at λ = 0.7892 (Arb: 1.00010014764303636; Selberg family crosses 1 at 0.7891108).
2. Solved his extremal problem (maximise J = f(0) + 2∫_0^θ f(α)α dα over even f ≥ 0, f^(0) = 1, f^ ≤ 0 on |ξ| ≥ λ) by SDP: f = polynomial × Gaussian in the Hermite basis, both sign conditions as sums of squares, Gram blocks matched exactly in the orthonormal basis (sdp.py). Validation: at θ = 1 it gives λ* = 0.60400 (d = 41), BGMM's 0.6039.
3. Certificate f: degree 122 polynomial × e^{−πx²/s²}, s = 11/10, 62 rational Hermite coefficients (cert.json): 0.98 × optimum (d = 61, λ = 0.7851) + 0.02 × a strictly interior point (d = 41), so both signs hold strictly.
4. Lemma T (DERIVATION §3): any even nonnegative Schwartz f with these properties yields an F ∈ C_c² satisfying Thm Main2 exactly as stated: F = f·h_R + c·(λR^(λx) + (λ/φ)R^(λx/φ)), φ the golden ratio; the two Selberg minorants give a strictly negative tail |F^| ≥ κ(λ/ξ)^6 (Diophantine bound ‖nφ‖ ≥ 1/(3n)) that absorbs the cut-off error.
5. Inputs (brief item 2): θ < 4/7 comes only from Conrey's density slope 8/7 (l.221-239, l.381-389). A slope c > 8/7 would substitute verbatim (θ < min(c/2, 1)), but none is known. 2610.08965's proportion of simple critical zeros cannot replace it: the proof needs the weighted moment Σ_{β≥1/2} T^{θp(2β−1)}, not a count. Details: DERIVATION §2.

## Certificate
`cd Riemann/Deep/SmallGap && uv run --no-project --with python-flint python3 cert.py cert.json` (8 s, 37 MB) → ALL CHECKS PASS (logs/cert.log):
- f > 0 on R: p(w) = P(√w) has no root in [0, ∞): exact Sturm count 0 and Arb root isolation (only real root ≈ −0.441).
- f^ < 0 on |ξ| ≥ 0.7851: q has no root in [w0, ∞), w0 ∈ Q, w0 < 2πs²λ² (Arb); Sturm count 0; its only real root 4.6861382 < w0 ≈ 4.6861393.
- J(f)/f^(0) = 1.0000219019375662111 ± 7·10^{−21} > 1, θ = 399993/700000 = 4/7 − 10^{−5} exactly.
Independent cross-check: `uv run --no-project --with mpmath python3 xcheck.py cert.json` (10 min, logs/xcheck.log): f^ by direct quadrature equals the Hermite formula to 1.6·10^{−61}; J/f^(0) by quadrature agrees to 20 digits.
Regenerate cert.json (bit-identical): `OMP_NUM_THREADS=1 uv run --no-project --with cvxpy --with numpy --with scipy python3 make_cert.py 61 11/10 0.7851 41 1e-6 0.02 cert.json` (16 s, 287 MB).

## Checks
- Reviewer 1 (Lemma T and reduction): PASS on G ≥ 0 and C_c², R bounds, Diophantine step, cut-off tail, F meeting Thm Main2 as stated. Its item on the input chain needed lines 221-239 and 322-389; reviewer 3 read them: Conrey enters only via Lemma ZeroDensity in the Hölder step, a slope c gives θ < min(c/2, 1), and a proportion cannot replace it (PASS).
- Reviewer 2 (certificate): PASS on the Fourier convention, the sign logic, the Sturm and Arb root-isolation code, the functional formula and θ; it reran cert.py: ALL CHECKS PASS.
- Float SDP table and the θ = 1 sanity run: DERIVATION §5, logs/sanity_theta1.log.

## Novelty
Not pre-empted as of 11 Nov 2026 (NOVELTY.md): arXiv math.NT listings Oct-Nov 2026, arXiv API, OpenAlex citations of CGdL, Zenodo, HAL, our files. Lamzouri's paper is one day old.

## Caveat
- The gain is small (0.7892 → 0.7851) and anticipated by Lamzouri's remark; a referee may call it routine optimisation.
- Rests on 2611.14486v1 (unrefereed) and its inputs: Conrey's density estimate (announced; proof in Baluyot's thesis), BGST's unconditional Montgomery theorem, Fujii's moments.
- Certified λ is 0.7851. Float thresholds keep falling with the degree (0.785901 at d = 21 to 0.785072 at d = 61, ≈ 0.7849 at d = 100), and we have no dual bound, so the method's floor is unknown; 0.7850 needs d ≳ 80 in better than double precision (beyond the 300 MB cap with cvxpy).
- Lemma T is ours (reviewed, not refereed). The sign checks trust python-flint (exact rationals; Arb root isolation), with Sturm as a second method.

RH relevance: N8

## Files
Riemann/Deep/SmallGap/: REPORT.md, DERIVATION.md, NOVELTY.md, cert.json (certificate), cert.py (exact/Arb checker), xcheck.py (mpmath cross-check), make_cert.py (SDP → certificate), sdp.py (SDP), scan.py and thresholds.py (λ* bisection), sanity_theta1.py (θ = 1), repro_selberg.py (Lamzouri's value); logs/: cert.log, xcheck.log, make_cert.log, repro_selberg.log, sanity_theta1.log, thresholds.log.
