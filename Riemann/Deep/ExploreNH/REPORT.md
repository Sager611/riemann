NO IDEA: non-Hermitian physics makes spectra real in only two ways, a positive (Krein-definite) metric or continuation from a Hermitian limit that meets no exceptional point; for ξ the first is Weil/de Branges positivity and the second is the existing de Bruijn–Newman collision-exclusion machinery, whose exceptional points pile up at t = 0⁻ (Rodgers–Tao).

# explore24: non-Hermitian and PT-symmetric physics → RH tree
Orca agent explore24. Worked 6 Nov 2026 14:51–15:00 UTC and 11 Nov 11:28–12:05 UTC (login outage in between). Every number below comes from the scripts in this folder; each runs in under 30 s and under 40 MB, single-threaded.

## Short answer
- **Exceptional point (EP).** An EP is a parameter value where two eigenvalues and their eigenvectors coalesce. For the de Bruijn–Newman flow H_t it is a collision of two real zeros. Its normal form is universal: gap² = 8(t − t_c), with no dependence on height, amplitude or the other zeros. The EP time of a close ("Lehmer") pair is exactly the Csordas–Smith–Varga bound. For Lehmer's pair at 7005, our 81-zero ODE gives t_c = −7.11203·10⁻⁴ and the CSV formula gives −7.11206·10⁻⁴. These are lower bounds only, moot since Λ ≥ 0.
- **Krein signature.** The Krein signature is the sign of the indefinite norm of an eigenvector; Krein's rule says only eigenvalues of opposite signature can collide and leave the real axis. For a zero of H_t it is the sign of the de Branges phase velocity. A signature that forbids collisions for t > t₀ exists ⟺ H_{t₀} is real-rooted ⟺ Λ ≤ t₀. Weil's form is a Pontryagin space (a space with finitely many negative directions) whose negative index κ equals the number of off-line zero pairs. "Unbroken PT symmetry" (a positive metric) is therefore Weil positivity (N1). The Bender–Brody–Müller metric Δ†Δ is positive for free and carries no RH content (2607.19067, R308).
- **Adversarial tests.** Every transfer fails for Davenport–Heilbronn and for a moved pair, as it must: off-line pairs show up as EPs at t > 0, for example t_e = 0.1507 for Lehmer's pair moved to σ = 1/2 ± 0.3. But every transfer only restates a known equivalent (Weil, de Branges/Lagarias, CSV, de Bruijn strip, Lax–Phillips) or is pre-empted in this repo (explore1, 5, 7, 15, 18, 22; R62, R71, R278, R327).

## Survey (arXiv listings, abstracts and two TeX sources; queries in .scratch/explore24/survey*.sh)
- **PT symmetry.** Bender–Boettcher 1998: x²(ix)^ε has real spectrum for ε ≥ 0. Dorey–Dunning–Tateo 2001 and Shin 2002 prove it with ODE/IM functional relations. That proof is a Hermite–Biehler (HB) argument: equal moduli at rotated points plus real zeros of a self-adjoint spectral determinant. R62 records ODE/IM as dead for Ξ.
- **Pseudo-Hermiticity.** Mostafazadeh 2002: a diagonalizable spectrum is real ⟺ the operator is quasi-Hermitian with a positive metric η. Krein-space PT: Langer–Tretter 2004 and Caliceti–Graffi–Sjöstrand 2005 show the spectrum stays real while a PT perturbation stays below half the gap. Survey 1610.08101.
- **Krein collisions.** Krein 1950, Gelfand–Lidskii 1955, MacKay 1986: in Hamiltonian–Hopf collisions only opposite signatures leave the axis. Kollár–Miller 2014 give the graphical signature (minus the slope of the eigenvalue curve). Kirillov's review 1907.04144. Generalized Nevanlinna classes N_κ: Krein–Langer 1977; Pontryagin spaces of entire functions: Kaltenbäck–Woracek.
- **EP perturbation theory.** Kato 1966 (Puiseux series), Berry 2004, Heiss 2012, Miri–Alù 2019 (EP sensing), the Petermann excess-noise factor, pseudospectra (Trefethen–Embree 2005; Krejčiřík–Siegl–Tater–Viola 2015).
- **Lasers and absorbers.** Fox–Li 1961 (open-resonator modes); Chong–Ge–Cao–Stone 2010 (coherent perfect absorbers are time-reversed lasers); scattering-matrix pole–zero pairs. Modular-surface scattering: Pavlov–Faddeev 1972, Lax–Phillips 1976 (R5.5, R58).
- **Non-Hermitian band theory.** Hatano–Nelson 1996; Yao–Wang 2018 (skin effect, generalized Brillouin zone); Schmidt–Spitzer 1960 (Toeplitz spectra with open boundaries); Beraha–Kahane–Weiss (zeros on equimodular curves); 2303.15186 (necessary and sufficient condition for a real open-boundary spectrum); 2302.13309.
- **Links to zeta.**
  - Bender–Brody–Müller 1608.03679; Bellissard's comment 1704.02644; the BBM note 1704.04705 and reply 1705.06767.
  - The metric completion 2607.19067 (SEEN, R308): η = Δ†Δ is positive, not coercive, and gives the dilation generator with spectrum ℝ, deficiency indices (∞, ∞).
  - Ahmed–Jain quant-ph/0407167: 2×2 pseudo-unitary ensembles with GUE spacing; no zeta content.
  - Pei 2504.07928: a KKR determinant via Krein's spectral-shift theorem, not Krein signatures.
  - 2408.01095: computing zeros by light scattering.
  - Suzuki's screw functions 2206.03682 and 2606.09096 (positivity of a Krein screw kernel ≡ RH; the TeX states no count of negative squares).
- **Links to de Bruijn–Newman.**
  - Csordas–Smith–Varga 1994 (Lehmer pairs; formula checked in Stopple's TeX 1508.05870).
  - Simonič 1612.08627 (Lehmer pairs via derivatives of Z).
  - Rodgers–Tao 1801.05914; Polymath15 1904.12438.
  - Planat 2609.37164 (the top non-real zero obeys d(y²)/dτ ≤ −2 − 4y²G).
  - Planat–Raymond 2611.06245 ("collision exclusion", Λ < 0.156821).
  - 2610.32694 (genus-two spectral curve).
  - Backward heat flow of polynomials as the Brown measure of A + i·semicircle (non-Hermitian free probability): Hall–Ho 2202.09660, Hall–Ho–Jalowy–Kabluchko 2308.11685, Kabluchko 2203.05533. This area is covered by R278 and R327.
- **Not found (1998–2026):** any paper that uses exceptional points, Krein collisions or PT thresholds on the de Bruijn–Newman flow or on Jensen polynomials.

## Transfers, tests, verdicts (derivations: [derivations.md](derivations.md))

### T1 Krein signature / Pontryagin index (N13 → N1, N2)
- **Claim.**
  - Weil's form W(f) = Σ_ρ f̂(γ_ρ)·conj f̂(γ̄_ρ) has one 2×2 block of signature (1, 1) per off-line pair γ, γ̄. So (test functions, W) is a Pontryagin space Π_κ with κ = number of off-line pairs (D3).
  - Likewise −ξ′/ξ ∈ N_κ, with κ = number of zeros of H_0 in the upper half-plane.
  - For H_t the Krein signature of a zero is the sign of the de Branges phase velocity. E_t = H_t + iH_t′ is Hermite–Biehler ⟺ H_t is real-rooted, and the heat flow keeps this only forward in t (D4).
  - Collisions need opposite signatures, or a real zero of E_t, which is where the metric degenerates.
- **Test.** We take the "Riemann–Siegel half" E = e^{iθ}ζ₃₃(1/2 + it), the main term of the approximate functional equation, and compute its phase velocity on [7003, 7009]. θ′ = 3.508.
  - The phase velocity is negative only on [7004.74, 7005.08].
  - At Lehmer's zeros it is −0.712 and +0.0905: opposite signatures, as the collision rule requires. The main term indeed loses this pair (T3).
  - Isolated zeros have +2.1 … +8.13.
  - The pair 7002.353 / 7002.692 also has opposite signs (−1.46 / +1.38) but survives (explore1). Opposite signature is necessary, not sufficient; the EP distance decides.
- **(a) DH and moved pair.** κ ≥ 1 for both, and the HB phase velocity turns negative at DH's zero 0.8085 + 85.70i (explore22's Re Λ′/Λ scans). No false proof.
- **(b) Restatement.** κ = 0 is Weil positivity. A definite phase is Lagarias' Re ξ′/ξ > 0 / de Branges (explore18 "definite Krein metric", explore22). A scalar function has no protection, because consecutive zeros always carry opposite graphical signature. BBM's η is positive for free. **Dead (restatement).**

### T2 EP normal form of de Bruijn–Newman collisions (N2): the numeric test
Method: the zero ODE ż_j = 2Σ_{k≠j} 1/(z_j − z_k) on zeros 6670–6750 (z = 2γ), with a static smooth-density tail, RK4, run backward to the collision (t2_ep.py, t2_csv.py).

| Quantity (Lehmer pair 6709/6710, z = 14010.126 / 14010.201) | Value |
|---|---|
| gap δ₀ / mean gap; bath S₀ = Σ(x̄ − z_k)⁻² | 0.075397 / 1.791; 0.6073 |
| EP time, full 81-zero ODE | **−7.112025·10⁻⁴** |
| EP time, pair-plus-bath closed form (D1) | −7.112023·10⁻⁴ |
| CSV λ = ((1 − 5Δ²g/4)^{4/5} − 1)/(8g) | −7.112064·10⁻⁴ (g = 1.21971) |
| −δ₀²/8 (bare EP) | −7.105884·10⁻⁴ |
| slope d(δ²)/dt at t = 0 / at the EP | 7.9863 (predicted 8 − 4S₀δ₀² = 7.9862) / **8.00000** |
| other EPs in the window (two-body model) | 34 real (next: −0.0507, −0.0613, −0.0652), 36 complex (e.g. \|Im t\| = 1.13–1.52) |

- **(a) Moved pair.** Replace the pair by x̄ ± iη and run the flow forward. Re-entry happens at an EP with t_e > 0:

| η (z-units; σ − 1/2 = η/2) | t_e (ODE) | de Bruijn η²/2 | screened (1/4S)log(1 + 2Sη²) |
|---|---|---|---|
| 0.05 | 0.001248 | 0.001250 | 0.001248 |
| 0.20 | 0.019542 | 0.020000 | 0.019529 |
| 0.60 | 0.150703 | 0.180000 | 0.149324 |

  DH: its pair 0.8085 + 85.70i is an EP at t ≈ +0.1764 (explore7). The EP map reports the broken phase correctly. It never forces t_e ≤ 0, so it cannot bound Λ from above.
- **(b) Restatement.** Real EPs are CSV 1994 (it agrees to 5 digits). Re-entry is de Bruijn's strip plus screening: Planat 2609.37164 and explore15's factor 2Sδ²/log(1 + 2Sδ²). Continuation without EPs is collision exclusion (Polymath15, Raymond, 2611.06245). EPs accumulate at 0⁻ (Rodgers–Tao): at height T the nearest has |t_c| ≈ δ²/8 with δ ≈ 4π/log(T/2π)·s and arbitrarily small GUE gaps s. So the continuation must stop at some t₀ > 0, behind a height cut-off that grows as t₀ ↓ 0. **Dead for Λ ≤ 0; known.**

### T3 EP sensitivity of approximate real-rooted models (N13, N8, N2 barrier)
- **Claim.** A real perturbation of size ε moves an isolated zero by ε/|Z′| along the line. A pair whose EP distance |Z(t_m)| ≈ |Z″|δ²/8 is below ε leaves the line, by about √(2(ε − |Z(t_m)|)/|Z″|).
- **Test (t3_afe.py).** At Lehmer's pair:
  - Z(t_m) = 0.003967, and |Z″|δ²/8 = 0.003970: the parabola is good to 0.06%.
  - The Riemann–Siegel remainder is 0.0700, so the N = 33 main term has no real zeros there. It has a zero at 7005.084141 + 0.076484i, i.e. σ = 1/2 ± 0.0765. The Puiseux prediction is 0.07686 (0.5%).
  - Zero 6712, for contrast, shifts linearly: −0.00856 against a predicted −0.00868.
- **(a)** This is a property of approximations: it applies to DH alike and proves nothing about either function.
- **(b)** Pre-empted: explore1 found the same complex zero (1/2 ± 0.07648 + 7005.08414i, √(2|R|/|Z″|) = 0.078) and 1.2% of main-term zeros off the line on [10³, 3·10⁴]. **Barrier only.** Any model whose error exceeds the smallest EP distance (→ 0 along Lehmer pairs) loses reality. This is the "barely true" Λ = 0 phenomenon.

### T4 Lasers, coherent perfect absorbers, scattering matrix (N7, N8)
- **Claim (D6).** For the modular surface, S(s) = Λ(2s − 1)/Λ(2s) is unimodular on Re s = 1/2 with resonances at ρ/2.
  - A zero or pole on the unitary axis (an absorber or a laser at threshold) ⟺ ζ(1 + 2ir) = 0. Its absence is the prime number theorem (Sarnak 2004, Eisenstein series).
  - RH ⟺ every resonance has the same decay rate 1/4, "one lasing threshold".
- **(a)** DH has no automorphic scattering realization, and the formalism puts no constraint on resonance positions.
- **(b)** Lax–Phillips / Pavlov–Faddeev (R5.5, R58). In de Branges form, a real zero of S = E^#/E is a real zero of E, i.e. an EP of Re E (D4). **Dead (restatement).**

### T5 Skin effect / generalized Brillouin zone (N8, N13) and PT thresholds for Jensen polynomials (N6)
- **Skin effect.**
  - An open-boundary Toeplitz spectrum is real when an imaginary gauge transformation makes the symbol real on a circle. The matrix is then similar to a Hermitian one (2303.15186).
  - ξ is real on the critical line, but so is DH. Reality of zeros needs a determinant/operator realization, i.e. Hilbert–Pólya. **Fails (a).**
  - The equimodular "two dominant terms" principle is Polymath15's large-height mechanism for t > 0 (a locally HB dominant half). It has no analogue at t = 0 (D7).
- **Jensen polynomials.** Krein-space threshold estimates (spectrum real while the perturbation is below half the gap) are weaker than the pointwise sign test at Hermite extrema used by GORZ, Savva, Holland and our wedges. The sharp PT threshold is the EP distance, i.e. that same sign test. No new range; not pursued.

## Ranking (RH relevance × chance × gain)
| Transfer | Node | Relevance | Chance | Gain | Verdict |
|---|---|---|---|---|---|
| T1 Krein index / signature | N13 → N1, N2 | high | ≈ 0 (it is RH) | none | Weil, de Branges, Lagarias |
| T2 EP normal form of collisions | N2 | medium | 0 for upper bounds | lower bounds only | CSV 1994 |
| T3 EP sensitivity of models | N13, N8 | low | – | barrier, pre-empted | explore1 |
| T4 laser / absorber / scattering matrix | N7, N8 | medium | 0 | prime number theorem only | Lax–Phillips, Sarnak |
| T5 skin effect / PT thresholds for Jensen | N8, N13, N6 | low | 0 | none | DH passes / IVT |

## Filters learned
- Every non-Hermitian reality theorem is one of two kinds:
  - a positive metric (→ N1 or de Branges; κ = number of off-line pairs);
  - continuation without EPs from a Hermitian limit (→ N2 collision exclusion).
- Exceptional points are invisible to holomorphic spectral invariants (D5). Taylor coefficients of H_t are entire in t; power sums and log-potentials are real-analytic through EPs. Λ is a signature threshold (the Hermite–Grommer Hankel index N₊(t)), so Lee–Yang/Fisher analyticity arguments cannot locate it.
- A scalar real function carries no Krein protection: adjacent real zeros always have opposite graphical signature. Protection needs a fixed definite metric.

## Suggested NOGO line
- explore24 (Orca, exploration: non-Hermitian/PT physics: exceptional points, Krein signatures and collisions, lasers and coherent perfect absorbers, skin effect → N2, N13, N1, N6, N8) NO IDEA: non-Hermitian reality comes only from a positive metric (for ξ, Weil's form is a Pontryagin space Π_κ with κ = number of off-line pairs, or the de Branges phase; BBM's metric is positive for free, 2607.19067) or from continuation without exceptional points (the dBN collision-exclusion machinery); a dBN collision is a square-root EP with universal gap² = 8(t − t_c), and a Lehmer pair's EP time is the Csordas–Smith–Varga bound (pair 7005: −7.11203·10⁻⁴ by an 81-zero ODE, CSV −7.11206·10⁻⁴); off-line pairs (DH, moved pair) are EPs at t > 0 (Lehmer's pair moved to σ = 1/2 ± 0.3: 0.1507), EPs accumulate at t = 0⁻ (Rodgers–Tao) and are invisible to every holomorphic spectral invariant; absorber/laser points of the modular scattering matrix are the prime number theorem (Sarnak), and the skin-effect reality criterion needs an operator (DH is real on the line too) [ExploreNH](Riemann/Deep/ExploreNH/REPORT.md)

## Files
- REPORT.md (this file), [derivations.md](derivations.md) (D1–D7), [outputs.txt](outputs.txt) (all printed results).
- t2_ep.py: zero ODE, EP time, slope, moved-pair re-entry. t2_csv.py: CSV formula. get_zeros.py + zeros_6670_6750.txt: zeros 6670–6750 (mpmath).
- t3_afe.py: main term vs Puiseux law at 7005. t1_phase.py, t1_phase2.py: phase velocity (graphical Krein signature) of the Riemann–Siegel half.
- Run: uv run --no-project --with mpmath --with numpy python3 <script> (in this folder).

RH relevance: N2 (primary), N13, N1.
