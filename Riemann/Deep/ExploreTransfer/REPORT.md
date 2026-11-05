NO IDEA: transfer operators reach ζ only through the modular surface's scattering poles s = ρ/2, which sit behind the cusp-form line Re s = 1/2; positivity and Dolgopyat/Baladi–Vallée cancellation give contraction only in front of that line (Re s > 1/2, i.e. ζ(w) ≠ 0 for Re w > 1, already known), while on Re s = 1/4 the Mayer operator is expanding (at heights 30 and 40, 14 of its 15 non-negligible eigenvalues have modulus 1.06–2.07); a flow whose closed orbits are the primes would need the conclusion as input.

# Exploration explore3: transfer operators and Dolgopyat cancellation → RH tree (2026-11-05)

Field: thermodynamic formalism; Ruelle–Pollicott resonances; Dolgopyat's oscillatory cancellation; Baladi–Vallée (Euclidean algorithms); Mayer's operator for PSL(2,Z) and its twists (Chang–Mayer, Fraczek, Pohl); Naud, Oh–Winter, Bourgain–Kontorovich/Zaremba; fractal uncertainty (Dyatlov–Zahl, Bourgain–Dyatlov).
Sources: arxiv.org site search and abs pages (the export API returned HTTP 429), one TeX e-print (2511.06513, theorem statements only), listing scans of math.DS 2026-10 and math.NT 2026-10/11. No Evidence-API call (explore2 reported today's quota exhausted). Our TREE/NOGO (R7.8, R48, R58, R302).
Three transfers; T1 and T2 tested numerically (Arb via python-flint, mpmath), T3 argued. None survives.

## Survey (2015–2026)
- **Mayer's identity (background).** Z(s) = det(1 − L_s²), L_s f(z) = Σ_m (z+m)^{−2s} f(1/(z+m)) (Mayer 1991). Efrat 1993: det(1 − L_s) holds the even Maass forms and the ζ-zeros s = ρ/2; det(1 + L_s) holds the odd ones. Lewis–Zagier 2001: the Eisenstein period function Σ_{m,n≥1}(mz+n)^{−2s} solves the three-term equation up to ζ(2s)(z+1)^{−2s}, so ζ(2s) = 0 is the eigenvalue-one condition.
- **Baumgartner 2511.06513 (Nov 2025).** Mayer's operator on Hölder spaces C^α([0,1]): essential radius ≤ λ₁(Re s + α); the eigenvalues outside are the analytic ones; ζ-zero eigenfunctions need α > 1 − Re s (3/4 on the line). Functional analysis only, no positivity in Re s < 1/2.
- **Bertozzi–Bonanno–Varandas 2603.11852 (Mar 2026).** Dolgopyat's method for the geodesic flow of the modular surface (non-compact); recovers Ratner's exponential mixing. The rate is set by the spectrum on Re s = 1/2; the ζ resonances lie behind it.
- **Bettin–Drappeau–Lee 2610.41973 and 2610.42008 (30 Oct 2026).** Baladi–Vallée machinery for cuspidal orbits of any cofinite Fuchsian group gives Gaussian laws for central values L(x, 1/2) of additive twists (Maass forms, Dirichlet-theta series). Before: Bettin–Drappeau 1903.00457 (Estermann central values, modular symbols, Dedekind sums), 1905.01954 (cotangent sums). Distributions of values, not zeros.
- **Farey/Gauss operators.** Bonanno 2211.11664 (eigenvalue one of generalised Farey operators as an infinite-matrix problem); Bonanno–Isola 0907.1471, 1607.03414; Alkauskas 1210.4083 (structure of the Gauss–Kuzmin–Wirsing spectrum).
- **Twists and deformations.** Fraczek–Mayer 1011.4441 and Bruggeman–Fraczek–Mayer 1201.2324: under the character deformation of Γ₀(4) the λ = −1 zeros stay on the line (Phillips–Sarnak), the λ = +1 zeros leave it. Pohl et al.: 1503.00525 (unitary twists), 1606.09109, 1709.00760 (non-expanding cusp monodromy), 2209.05927 (strict transfer operators for all geometrically finite orbisurfaces).
- **Dolgopyat strips for thin groups.** Naud 2005; Oh–Winter 1410.4401 (uniform in square-free level); Magee–Oh–Winter 1412.4284, 1601.03705; Calderón–Magee–Naud 2407.21506 (random Schottky, optimal gap); Soares 2610.17276 (Fourier decay of Patterson–Sullivan measures, Oct 2026). Strips sit at Re s > δ − ε with δ < 1; no L-function.
- **Fractal uncertainty.** Dyatlov–Zahl 1504.06589, Bourgain–Dyatlov 1612.09040: an essential gap beyond the unitarity axis for every convex co-compact surface, needing δ < 1. Lattices have δ = 1 (already on the tree as R7.8).
- **Zaremba.** Bourgain–Kontorovich 1107.3776 (density one); Shkredov 2603.14116 (prime denominators); X. Zhang 2605.02518 (claims the full conjecture, May 2026). No RH link.
- **Euclidean Dirichlet series.** Baladi–Vallée cs/0307062: Dolgopyat bounds for (1 − H_s)^{−1} in Re s > 1 − ε, |t| ≥ t₀, giving Gaussian laws with power savings; the plain series is Σ_{u/v} v^{−2s} = ζ(2s−1)/ζ(2s).

## Transfers and tests

### T1 (N12, N13, N7): a Dolgopyat strip for Mayer's operator and the scattering zeros s = ρ/2
- **Object.** The matrix of L_b in the basis (z−1)^k on |z−1| < 3/2: M[n,k] = (−1)^{n+k} Σ_{j≤k} C(k,j)(−1)^j (2b+j)_n/n! · ζ(2b+j+n, 2). This is Mayer's continuation to every b ≠ (1−n)/2.
- **Sanity.** At b = 1 the eigenvalues are 1, −0.303663002898594, 0.10088451 (Wirsing: −0.3036630028987).
  - At b = ρ/2: |det(1 − L_b)| = 5.5·10⁻¹⁶, 8.9·10⁻¹⁴, 7.8·10⁻¹¹ for γ₁, γ₂, γ₃ (K = 50), against 0.07–0.17 at +0.05i. |det(1 + L_b)| = 1.68–2.05 there.
  - Cusp forms: det(1 − L) = 1.4·10⁻⁹ at the even r = 13.7798; det(1 + L) = 8.7·10⁻¹³ at the odd r = 9.5337. So the ζ-zeros sit in det(1 − L) only (Efrat).
- **The Eisenstein matrix element is the scattering coefficient.** Each rational in (0,1) has two continued fractions, so ⟨(1 − L_b)^{−1}1, δ₀⟩ = 2ζ(2b−1)/ζ(2b).
  - It matches to ≤ 10⁻¹¹ at nine points, including b = 0.3 + 5i and 0.25 + 9i in the continued region, and to 9·10⁻⁸ at b = 1/2 + 17.74i, next to a cusp form (K = 90).
  - It stays finite at the even cusp form (the pole cancels in this matrix element) and grows like 1/distance at ρ₁/2.
  - Any bound on it is a bound on ζ(2b−1)/ζ(2b) in the critical strip; it gives no operator handle.
- **Contraction map.** Spectral radius r(σ,t) on σ ∈ [0.2, 1.1] × t ∈ [0.5, 27] (grid 0.05 × 0.5, file t1_grid.tsv). Multiprecision QR at t = 25 agrees within about 2%.
  - r(1,t) ∈ [0.38, 0.56] for 9 ≤ t ≤ 27: Dolgopyat's regime, flat in t.
  - At the same heights r(1/2,t) ∈ [0.77, 1.05] and r(1/4,t) ∈ [0.99, 2.11].
  - The contraction boundary σ*(t) (r < 1 for all σ ≥ σ*) lies in [0.45, 0.55] at 38 of 39 heights in [8, 27]. The exception is t = 10.5, just below the ζ event at 10.511, where it dips to 0.25.
  - For 1 ≤ t ≤ 6.5 the operator contracts down to σ = 0.2, but there are no zeros there anyway.
- **Branches.** Through every eigenvalue-one point tested, |λ(σ+it)| decreases strictly in σ and crosses 1 transversally:
  - ρ₁/2: 1.248, 1.000, 0.582 at σ = 0.05, 0.25, 0.75;
  - odd cusp 9.53: 1.871, 1.000, 0.728 at σ = 0.05, 0.50, 0.75.
  - At ρ₃/2 another eigenvalue (modulus 1.63) dominates.
- **No branch is unimodular along its line.** For Δt ∈ [−1, 1]:
  - the cusp branches vary by up to 0.068 (even) and 0.236 (odd) on σ = 1/2;
  - the ζ branches vary by up to 0.32 (ρ₁) and 0.57 (ρ₂) on σ = 1/4.
  - Even Selberg's theorem (cusp zeros on the line) is invisible branch by branch. It comes from the Laplacian, not from the operator. The line σ = 1/4 is not a symmetry of the family either: the Selberg functional equation is b ↔ 1 − b.
- **The bulk** (multiprecision QR; converged K = 180/220 at t = 30 and 220/260 at t = 40).
  - The eigenvalues of modulus ≥ 0.2 number about t/5: 1, 2, 3, 4, 5, 6, 8–9 at t = 5, 10, 15, 20, 25, 30, 40.
  - On σ = 1/2 they cluster at the unit circle (7 of 9 in [0.96, 1.023] at t = 40).
  - On σ = 1/4 they exceed 1: all 6 at t = 30 (1.06–1.95) and 8 of 9 at t = 40 (1.17–2.07).
  - On σ = 1 they lie in 0.21–0.47.
  - The scale fits damped-quantum-map behaviour, e^{Λ(1/2−σ)} with Λ = π²/(6 log 2) = 2.373: 1.81 at σ = 1/4 and 0.305 at σ = 1.
- **Verdict: DEAD.**
  - No region {Re b > σ₀, |Im b| ≥ t₀} with spectral radius below 1 can have σ₀ < 1/2, because the cusp forms put eigenvalues ±1 on Re b = 1/2 at every height. Numerically σ*(t) ≈ 0.5 for t ≥ 8.
  - So every uniform contraction or pole-free statement for L_b lives in Re b > 1/2. There Z(b) ≠ 0 for b ≠ 1 is already Selberg's theorem; translated, ζ(w) ≠ 0 for Re w > 1.
  - The ζ-zeros are eigenvalue-one crossings of a block of about t/5 eigenvalues that is expanding on σ = 1/4. They occur at rate about (1/π) log(t/π) per unit height, against t/6 for the cusp forms.
  - So no norm, Lasota–Yorke or Dolgopyat bound can exclude eigenvalue one in 1/4 < Re b < 1/2.
  - This quantifies R48 (sign-changing weights) and R58 (no positivity at complex s).

### T2 (N7): Baladi–Vallée strips for Dirichlet series with arithmetic coefficients
- **Identities.** Checked by direct summation over rationals with v ≤ 8000 (B and C to 8 digits, A to 10⁻⁵ because of its slow tail):
  - (A) Σ v^{−2s} = ζ(2s−1)/ζ(2s);
  - (B) Σ χ₋₄(v) v^{−2s} = L(2s−1, χ₋₄)/L(2s, χ₋₄);
  - (C) Σ χ₋₄(u) v^{−2s} = D(2s)/L(2s, χ₋₄), with D(w) = Σ_{m ≡ 1,2 (4)} m^{−w}. It has no pole at s = 1 (value 1.5102) and a pole at 0.25 + 3.0105i, from the first zero 1/2 + 6.0209i of L(·, χ₋₄).
- **General case.**
  - Twists through (u, v) mod N give finite sums of L(2s−1,·)/L(2s,·) and L(2s,·)/L(2s,·).
  - Twists through the whole matrix mod N give Kloosterman–Selberg series Σ_c S(m,n;c) c^{−2s}. Their poles are the Maass spectrum of Γ(N) (Selberg 1965), which is the cusp wall again.
  - Von Mangoldt weights appear only on the periodic-orbit side: prime geodesics of norm ε_d² with class-number multiplicities, not primes.
- **Why it stops.**
  - Möbius enters only through coprimality of (u, v), and scaling (u, v) ↦ (du, dv) multiplies v^{−2s} by d^{−2s}. So the L-function sits at 2s, and its zeros sit at Re s ≤ 1/2, half a unit behind the leading pole s = 1.
  - Strip width: Dolgopyat's proven ε is tiny and non-explicit. Even the best possible strip (numerically σ* ≈ 0.5) only restates Re w > 1.
  - Height: r(1,t) stays at 0.38–0.56 for 9 ≤ t ≤ 40. Meanwhile r(1/4,t) grows from about 1.2 (t = 10) to about 2 (t = 25–40), and the number of eigenvalues above 1 at σ = 1/4 grows like t/5.
- **Verdict: DEAD.** Arithmetic coefficients are reached, but only at the doubled argument; nothing beyond the classical half-plane.

### T3 (N7, N10): a flow whose closed orbits are the primes (argued)
- **Set-up.** Suppose a suspension over a real-analytic expanding Markov map T (countably many full branches, as for Gauss) has roof r and primitive closed orbits in bijection with the primes, of length log p.
  - Then its Ruelle zeta is the Euler product ζ.
  - Ruelle's trace formula gives ζ(s) = det(1 − L_{s,1})/det(1 − L_s), with L_s f = Σ_h e^{−s r∘h} f∘h and L_{s,1} f = Σ_h e^{−s r∘h} h′ f∘h.
  - So ζ-zeros are eigenvalue-one points of L_{s,1}.
- **Positivity suffices.** |L_{s,1}ⁿ f| ≤ (positive operator)ⁿ|f|, so L_{s,1} contracts for Re s > σ₁, where P(−σ₁r − log|T′|) = 0. If r ≤ A log|T′|, then σ₁ ≤ 1 − 1/A.
  - Such a realization would give ζ ≠ 0 on Re s > 1 − 1/A without Dolgopyat, on the half-plane where both operators are defined.
  - A ≤ 2 would give RH. Conversely, the zeros on the line force A ≥ 2 if the operators reach Re s = 1/2.
- **Surfaces are excluded.** For compact hyperbolic surfaces ζ_R(s) = Z(s+1)/Z(s). Then ζ = ζ_R would force Z(s) = Π_{k≥0} ζ(s+k)^{−1}, which has poles at ρ − k, but Z is entire. This is the sign obstruction again.
- **Why it stops.** Building r with exact lengths log p and Hölder or analytic control at depth k needs primes in intervals of relative size θ^k around prescribed targets, at every depth, together with a global bijection. That is information about primes at all scales, at least as strong as the strip it would yield. No construction is known; the primon gas and Bost–Connes models (R48) have no expanding base.
- **Verdict: CIRCULAR.** The reformulation is "RH ⟸ such a realization with A ≤ 2".

### Also checked (argued)
- FUP (R7.8) needs δ < 1. As Hecke triangle groups Γ_w → PSL(2,Z) (w → 2), δ(w) → 1 and the gap closes.
- Non-unitary twists (Fedosova–Pohl) raise the pressure, which moves contraction to the right, not the left. Unitary twists keep the cusp wall.
- Character deformation of Γ₀(4): the λ = +1 zeros leave the line. A deformation that keeps a functional equation but loses the Euler product moves zeros off the line (our Davenport–Heilbronn filter, R2), so it gives no deformation argument.
- Bettin–Drappeau–Lee central-value laws (N9, N5): they give distributions of L(x, 1/2) over cusp orbits, which neither bound moments nor locate zeros. The Nyman–Beurling distance d_N is extremal over coefficient vectors, so laws for Vasyunin cotangent sums cannot decide d_N → 0.

## Ranking (RH relevance × chance × gain)
| Transfer | Node | RH relevance | Chance | Gain if it worked | Verdict |
|---|---|---|---|---|---|
| T1 Mayer/Dolgopyat strip → zeros ρ/2 | N12, N13, N7 | high (exact ζ link) | ≈ 0: contraction only for Re b > 1/2; expanding at 1/4 | nothing beyond Re w > 1 | DEAD |
| T2 Baladi–Vallée → arithmetic Dirichlet series | N7 | medium | ≈ 0: L-function at 2s | nothing beyond Re w > 1 | DEAD |
| T3 prime-orbit flow + positivity/Dolgopyat | N7, N10 | high | ≈ 0: needs prime data at all scales | quasi-RH, RH if A ≤ 2 | CIRCULAR |
| FUP strip (R7.8) | N7 | high | 0 for lattices (δ = 1) | — | known |
| Bettin–Drappeau–Lee laws | N9, N5 | low | — | distributions only | DEAD |

## Barrier (proposed NOGO line, not written by explore3)
- explore3 transfer operators (Mayer, Dolgopyat, Baladi–Vallée, twists, prime-orbit flows): ζ enters only at s = ρ/2 behind the cusp-form line; contraction only for Re s > 1/2 (= Re w > 1); on Re s = 1/4 the Mayer operator is expanding (14 of 15 non-negligible eigenvalues of modulus 1.06–2.07 at t = 30, 40). [report](Riemann/Deep/ExploreTransfer/REPORT.md)

## Files (Riemann/Deep/ExploreTransfer/)
Run any script with `uv run --no-project --with python-flint --with numpy python <file>`. t1_mayer.py and t2_twist.py need mpmath instead of python-flint (`--with mpmath --with numpy`). All runs were single-threaded under alarms of 1800 s or less.
- t1_mayer.py, t1_mayer_out.txt: mpmath Mayer matrix; Gauss–Kuzmin–Wirsing sanity check; det(1 ∓ L) at ρ/2 and at cusp forms (K = 30, 40, 50).
- t1_flint.py, t1_grid.tsv: Arb version; spectral-radius grid r(σ,t) (`grid` mode, 875 s). Grid values with t ≥ 28 are unreliable (double-precision eigenvalues): r(1,32) = 1.6 is impossible, since |L_{1+it}ⁿf| ≤ L₁ⁿ|f| gives r ≤ 1. t1e_highprec.py redoes t = 30 and 40 in multiprecision.
- t1b_branches.py, t1b_out.txt: Eisenstein matrix element against 2ζ(2b−1)/ζ(2b); branch moduli along horizontal lines; the ρ₁ branch along σ = 1/4.
- t1c_unitarity.py, t1c_out.txt: branch moduli along σ = 1/2 (cusp forms) and σ = 1/4 (ζ zeros).
- t1d_bulk.py, t1d_out.txt: eigenvalue cloud at t = 5–25.
- t1e_highprec.py, t1e_out.txt: multiprecision QR check at t = 25 and the bulk at t = 30 and 40.
- t2_twist.py, t2_out.txt: identities A, B, C; the pole of (C) at the first zero of L(·, χ₋₄).
- Working notes: .scratch/explore3/STATE.md, .scratch/explore3/survey.md.

RH relevance: N7, N12, N13 — no advance; barrier recorded (transfer-operator contraction is capped at Re w > 1).
