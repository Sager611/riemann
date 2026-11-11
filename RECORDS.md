# Records

Constants, ranges, exponents and proportions (AGENTS.md: records, not discoveries); "half-plane" marks dependence on the Lean-checked zero-free half-plane. Pages as in [DISCOVERIES.md](DISCOVERIES.md).

## On the RH tree

- (N1) Superseded record: Weil's quadratic form certified positive on test functions supported in [−1.55, 1.55] (was 1.5; superseded by Chekmachev's 1.925, Zenodo 23848597, audited sound). [details](discoveries/weil-positivity.md)
- (N10 (barrier)) Record: the integer Chebyshev constant of [0, 1] is at least 0.4222355 (was 0.4222286). [details](discoveries/integer-chebyshev.md)
- (N4, half-plane) Record: Robin's inequality verified for all n ≤ 10^(10^46.79), Nicolas' for primorials with p_k ≤ 1.78·10^48 (was ours 10^(10^38.71); Polak 10^(10^22.85)). [details](discoveries/robin-nicolas-ranges.md)
- (N10, half-plane) Explicit power-saving Mertens bound: |M(x)| ≤ 9.77·10^4 x^0.9382 and ≤ 1.54·10^4 x^0.9409 for all x ≥ 1 (11 Nov: + RH to 3·10^12 on a hybrid contour; was ours 1.84·10^6 x^0.9382, the first, 5 Nov; Chirre–Helfgott input refereed, [CHAudit2](Riemann/Deep/CHAudit2/REPORT.md): robust to its square-free chain and to Platt's sums up to ≈ 7,690). [details](discoveries/mertens-power-saving.md), [MertensCH](Riemann/Deep/MertensCH/REPORT.md)
- (N9) Record: |ζ(1/2+it)| ≤ 0.568 t^{1/6} log t for all t ≥ 3, and ≤ 0.4821 t^{1/6} log t for t ≥ 8.97·10^17 (was 0.611, 0.566); holds with Hiary–Patel–Yang's refereed lemmas, 0.56797; Revers's lemmas audited correct, [RevAudit](Riemann/Deep/RevAudit/REPORT.md)). [details](discoveries/zeta-critical-line-constant.md)
- (N9) Record: |ζ(1/2+it)| ≤ 54.6 t^{27/164} for all t ≥ 3 (was 66.7, Patel–Yang 2024, whose proof drops a factor ≈ 1.10 and gives 68.54 as computed). [details](discoveries/sub-weyl-constant.md)
- (N10) Conditional record: lim sup M(x)/√x > 2.07390 and lim inf < −2.07641, given Weisz's published zero data, as his 2.06175 / −2.06471 is. [details](discoveries/mertens-limsup.md)
- (N11, half-plane) Record: all imaginary quadratic fields with class number h ≤ 10,000 determined without RH or GRH: 410,810,155 fields, largest |d| 28,154,696,707 (was h ≤ 5000, ours; 2000, Oct 2026). [details](discoveries/class-number-10000.md)
- (N6) Superseded record (by the explicit Kim–Lee wedge below): Riemann's Jensen polynomials J^{d,n} are hyperbolic whenever n² log(n+2) ≥ 24 d³ (Savva's two-thirds wedge had no number); supersedes our n³ log²(n+2) ≥ 7000 d⁵. [details](discoveries/jensen-two-thirds-wedge.md)
- (N6, half-plane) Explicit Kim–Lee wedge: Riemann's Jensen polynomials J^{d,n} are hyperbolic for every n ≥ 9.99·10^11 and d ≤ 6.45 n²/log² n (was n² log(n+2) ≥ 24 d³). [details](discoveries/jensen-kim-lee-wedge.md)
- (N8) At least 0.6724 of the zeta zeros are simple, on the line and isolated from both neighbours (every c ≤ 0.00476; was 0.6552, c unspecified). [details](discoveries/isolated-zeros.md)
- (N8) Large gaps between consecutive zeros on the critical line: Λ > 3.2648 unconditionally (Bui–Milinovich's criterion with a degree-14 amplifier, Arb certificate; was 3.18, Bui–Milinovich 2018). [GapLambda](Riemann/Deep/GapLambda/REPORT.md)
- (N8) Small gaps: a positive proportion of consecutive ordinates have γ_{n+1} − γ_n ≤ 2πλ/log T for every λ ≥ 0.7851, unconditionally (Lamzouri's criterion with an SDP-optimised test function, exact certificate; was 0.7892, Lamzouri 2611.14486, audited: [LamzAudit](Riemann/Deep/LamzAudit/REPORT.md); its density input is proved only in Baluyot's 2017 thesis, refereed by us: [Baluyot](Riemann/Deep/Baluyot/REPORT.md)). [SmallGap](Riemann/Deep/SmallGap/REPORT.md)
- (N8) Jerby's Gram-point repulsion conjecture (|Z′(g_n)| > 4|Z(g_n)| at every isolated bad Gram point; arXiv 2406.06548, 2310.14415) holds for every n < 2.54·10^8 (was 2.5·10^7, Jerby 2610.21952); smallest ratio 4.21553 at n = 69,784,845 (Arb-certified; the range is double-precision screening). [ConjFals](Riemann/Deep/ConjFals/REPORT.md)
- (N8, N12) At least 2.06% of the zeros of a combination of three Dirichlet L-functions lie on the critical line (was 4.91·10⁻⁹); 34.71% for two (was ours 34.49%). [details](discoveries/linear-combinations-three.md)
- (N12, N8) At least 6.8·10⁻⁴ of the zeros of a combination of four Dirichlet L-functions lie on the critical line (was 3.7·10⁻⁹), and c/N for N summands. [details](discoveries/selberg-linear-combinations.md)
- (N8, N12) At least 7.005% of the zeros of every GL(2) newform L-function lie on the critical line (was 6.96%); arXiv 2508.11108's Question 5 settled. [details](discoveries/gl2-critical-line.md)
- (N5) Báez-Duarte distance certified: d_N² for every N ≤ 1280 (Vasyunin's Gram summed exactly, ball width < 10⁻²³); d_N² log N − C (C = 2 + γ − log 4π) changes sign 15 times on [2, 1280], first below C at N = 41, below on [584, 1280] (was N ≤ 128 rigorous, Zenodo 24088197; d_300² interval, Zenodo 21898513; 2611.03924 in double precision). [NBRate](Riemann/Deep/NBRate/REPORT.md)

## Side records

- Record: counterexample to Gale's conjecture in R³ with ratio 0.8906211 (previous 0.8892812). [details](discoveries/gale-r3.md)
- Records: stronger counterexamples to Gale's conjecture in R⁴ (0.9164674, was 0.9131485), R⁵ (0.9283700, was 0.9257636) and R⁶ (0.9373384, was 0.9369672). [details](discoveries/gale-r4-r5.md)
- Record: non-separable homothetic triangles with covering factor 1.0633023 (was 1.0631883), and a proof that λ_simplex ≥ 1.0633083. [details](discoveries/tricover.md)
- Records: quaternary covering codes K₄(10,3) ≤ 768 (was 832) and K₄(11,3) ≤ 1920 (was 2048). [details](discoveries/covering-code-k4-10-3.md)
- Record: a body of volume 1.8755015 moves through the 3D snaky corridor K₄ (was 1.8558404). [details](discoveries/sofa3d.md)
- Record: a convex cover of area 0.2453479 contains a copy of every planar arc of length 1 (Moser's worm problem; was 0.2463322). [details](discoveries/moser-worm.md)
- Record: kissing number in dimension 38 is at least 624020 (was 611522): a known 3-modular lattice gives a far larger equator. [details](discoveries/kissing-38.md)
- Record: ξ(K₆) ≥ 33 — a 3-colouring of K₃₃ in which every K₆ sees all three colours (was 25). [details](discoveries/xi-k6.md)
- Record: the ambidextrous moving sofa constant is at most 19/10 (was 2.2195, from an unverified announcement). [details](discoveries/ambidextrous-sofa.md)
