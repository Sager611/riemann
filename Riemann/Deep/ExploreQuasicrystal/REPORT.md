NO IDEA: Fourier-quasicrystal / Lee–Yang theory classifies constant-density zero sets of polynomials on a finite torus, and ζ's zero measure is outside that class on every count (density log(T/2π)/2π, Fourier masses growing like e^{X/2}, a continuous archimedean part, and a prime-power spectrum that on any finite prime torus forces arithmetic progressions); the Lee–Yang / Hermite–Biehler mechanisms that do reach ζ-type functions are either crystalline toys (Euler-product torus models are real-rooted only for small X, where from height 10⁴ on their zeros are measurably more rigid than ζ's) or exactly as strong as the zero-free half-plane (Lagarias' differenced ξ, threshold h = sup Re ρ − 1/2).

# Exploration explore5: Fourier quasicrystals, crystalline measures, Lee–Yang varieties → RH tree (2026-11-05)

Brief: `.scratch/orca/explore5-brief.md`. Sources: arXiv abs pages and e-prints (export API throttled, HTTP 429, from 12:45), arXiv listings 2026-10 and 2026-11 for math.CA, math.NT, math.CV, math.SP, math.MP, and our NOGO/TREE/ExploreStatMech. Not repeated: NOGO R2.5, R2.6, R7.9, R34, R13.7, explore1 transfer B (Lee–Yang lift of the Riemann–Siegel main term). Derivations: `derivations.md`.

## Survey (one page)
- **Classification (Dyson's "classify 1-D quasicrystals" is solved for integer masses, and ζ is not in the class).**
  - Olevskii–Ulanovskii 2009.12810: a Fourier quasicrystal (FQ) with unit masses is the zero set of an exponential polynomial with finitely many frequencies.
  - Alon–Cohen–Vinzant 2303.03201: every real-rooted exponential polynomial is a Lee–Yang (LY) polynomial restricted to a positive Kronecker line. So all ℕ-mass FQs come from the Kurasov–Sarnak construction (2004.05678).
  - Wei 2610.10652 (Oct 2026): positive-integer-mass FQs are exactly the real sections of effective LY cycles. The support is the common zero set of finitely many trigonometric polynomials, and a finite-dimensional torus is recovered from the spectrum.
- **Statistics and geometry.**
  - Alon–Vinzant 2307.13498: gap distributions of ℕ-FQs exist and are absolutely continuous. There are sequences tending to Poisson, and random LY polynomials tending to CUE.
  - Higher dimensions: Alon–Kummer–Kurasov–Vinzant 2407.11184; Alon–Kummer 2507.16029 (periodic "lighthouse" hypersurfaces are LY); Lawton–Tsikh 2302.07464.
- **Beyond FQs.**
  - Favorov: a crystalline measure that is not an FQ (2401.01121); uniqueness theorems (2006.07610).
  - Strip analogues (2408.09563; Favorov–Değer 2605.10766, May 2026): for measures on a strip with pure-point Fourier transform, Σ|b_γ|²δ_γ "has exponential growth" (their wording).
  - Favorov 2411.07190: an absolutely convergent Dirichlet series with only real zeros is a finite product of sines under a simple condition.
- **Direct links to ζ.**
  - Gonçalves–Kerr–Klurman–Vedana 2610.26705 (Oct 2026): a real-rooted trigonometric polynomial shares at most 2/3 of its zeros with ζ on the line (crystallinity via KS + ACV). This is a separation theorem, not an approach.
  - Fourier uniqueness pairs: Bondarenko–Radchenko–Seip 2005.02996 (zeros of ζ and log n), Kulikov–Nazarov–Sodin 2306.14013, Lysen 2509.17600 (critical pairs beyond BRS).
- **One-variable Lee–Yang = Hermite–Biehler for ξ.**
  - Taylor 1945 and Lagarias math/0601653: A_h = ½(ξ(s+h) + ξ(s−h)) and B_h have only simple critical-line zeros for |h| ≥ 1/2, with trivial (equal) spacing.
  - Gonek 0704.3448: finite Euler products with a functional equation. His RH-by-construction model is P_X(s) + χ(s)P_X(s̄); the analytic P_X(s) + χ(s)P_X(1−s) he calls the "wrong AFE", because |P_X(1−s)| reaches exp(X^{σ−1/2}/log X).
- **Listings 2026-10 and 2026-11, other relevant items:**
  - Lamzouri 2610.26841: at least 31 % of the zeros of DH-type combinations lie on the line with odd multiplicity.
  - Hofgard 2610.09485: LY tensors.
  - 2611.02838: boundary zeros of stable polynomials in 𝔹₂.
  - 2610.36606: arithmetic Poisson formulas.

## T-A (N13, N1): which FQ hypotheses fail for the zeta zero measure, and by how much (`t1_audit.py`, `t1_out.txt`)
Under RH, the explicit formula says the Fourier transform of μ = Σδ_γ consists of:
- atoms −Λ(n)n^{−1/2} at ±log n,
- the pole term 2cosh(u/2),
- a continuous archimedean part.

| hypothesis (OU / KS / ACV / Wei) | ζ | gap |
|---|---|---|
| constant density Σd_jℓ_j/2π | log(T/2π)/2π | matching with one degree per prime needs θ(P) = log(T/2π): P = 7, 17, 31, 67, 251 at T = 10³, 10⁶, 10¹², 10²⁴, 10¹⁰⁰ — the torus must grow with T |
| |μ̂| translation-bounded | atomic mass in (X, X+1] = 2345.5 at X = 15, equal to the pole term 2(e^{(X+1)/2} − e^{X/2}) to 2·10⁻⁴; growth rate 1/2 | exponential; temperedness needs cancellation against 2cosh(u/2), and that cancellation is RH (residual max 3.44 for X ≤ 17; unconditionally only O(e^{(0.3734+ε)X})) |
| pure point Fourier side | continuous Γ-part and the pole term | not crystalline, even after subtracting the mean density |
| spectrum in a rank-r semigroup ⟨m,ℓ⟩ | {log n}: count e^X | rank-11 (31-smooth): 4.4·10⁶ frequencies ≤ 30 against e³⁰ = 1.1·10¹³; infinitely many variables are needed |

**Rigidity lemma (ours, elementary; `derivations.md` D1).** In KS the coefficient of μ̂_P at ⟨m,ℓ⟩ is ∝ ⟨m,ℓ⟩c_m, where log P = Σc_m z^m on 𝔻ⁿ. If the spectrum lies on prime powers {k log p} (as in ζ's explicit formula), then P = c·Π_j Q_j(z_j). LY then forces unimodular roots, so the zero set is a finite union of arithmetic progressions.
- On any finite prime torus, ζ's prime-power-only spectrum therefore gives periodic zero sets.
- Non-periodic structure needs composite or reflected lines, which ζ does not have, or the archimedean factor, whose phase θ(t) is not a torus flow. Linearising θ is explore1's dead Riemann–Siegel model.

**Toy (argued, not computed; D2).** A finite Euler product with exact conductor, F_P = (P#)^{s/2}/ζ_P(s) + (s ↔ 1−s), is Lee–Yang for every P. It satisfies RH by Hermite–Biehler, since each factor |p^{s/2} − p^{−s/2}| grows across the line.
- Its prime lines carry the 1/ζ sign, +Λ(n)/(2π√n). ζ has −Λ(n)/(2π√n).
- The ζ sign would need the pole (1 − p^{−1/2}z)^{−1} on every prime axis, i.e. no polynomial. The ζ_P version fails Hermite–Biehler at p ≤ 7 and at the poles on Re s = 1.

**DH test.** The lemma uses the Euler product (DH's −f'/f has coefficients on all n), but it kills the approach; it proves nothing.

**Verdict:** dies (no FQ mechanism on the prime torus).

## T-B (N13, N12): Euler-product torus models with the true Γ-factor — where the zeros sit (`t2_models.py`, `t2_out.txt`, `t2b_euler_part.py`, `t2b_out.txt`)
**Models.** P_X = exp Σ_{n≤X} Λ(n)n^{−s}/log n.
- ζ: G_X = E(s) + E(1−s), with E = π^{−s/2}Γ(s/2)P_X. E never vanishes: it is the LY / Hermite–Biehler candidate.
- DH: H_X = Φ_X(s) + \overline{Φ_X(1−s̄)}, with Φ_X = (5/π)^{s/2}Γ((s+1)/2)[c₁P_X(s,χ) + c₂P_X(s,χ̄)], χ mod 5, χ(2) = i, Davenport–Heilbronn's κ.

**Method.**
- On Re s = 1/2 the real zeros coincide with Gonek's model.
- Windows [T, T+200].
- Real zeros come from phase crossings. All zeros come from the argument principle on 0.02 ≤ σ ≤ 0.98 ("full") and on 0.3 ≤ σ ≤ 0.7 ("inner").
- Zeta zeros are from Odlyzko's table. Gap variance is of unfolded gaps.

| T | ζ: gap var. | X | ζ-model off-line full/inner | ζ-model gap var. | median distance to ζ zeros (spacings) | DH-model off-line full/inner | DH Euler-part zeros, Re s > 1/2 |
|---|---|---|---|---|---|---|---|
| 10³ | 0.153 | 10 | 0/0 | 0.132 | 0.066 | 12/12 | 8 |
| | | 20 | 0/0 | 0.177 | 0.043 | 14/6 | — |
| | | 30 | ≈16/0 | 0.185 | 0.046 | 22/10 | 12 |
| 10⁴ | 0.168 | 10 | 0/0 | 0.078 | 0.116 | 14/≈12 | 8 |
| | | 30 | ≈8/0 | 0.146 | 0.064 | 18/10 | 14 |
| | | 50 | ≈20/0 | 0.173 | 0.060 | ≈30/14 | — |
| 5·10⁴ | 0.158 | 10 | 0/0 | 0.051 | 0.132 | 14/≈14 | 8 |
| | | 30 | 0/0 | 0.117 | 0.082 | 22/≈12 | 11 |
| | | 100 | ≈66/0 | 0.165 | 0.057 | ≈86/≈16 | — |
| 7·10⁴ | 0.167 | 30 | 0/0 | 0.111 | 0.089 | 16/8 | 10 |
| | | 100 | ≈56/0 | 0.158 | 0.064 | ≈86/≈22 | — |

Further rows:
- X = 2: ζ-model gap variance 0.005–0.017, i.e. a crystal.
- X = 5: both models have no off-line zeros.
- "≈": the unwrapped phase made a step > 1 rad (a zero close to the contour), so the count is less certain. X ≥ 300 is unreliable throughout. Full table in `t2_out.txt`.

**Findings.**
1. The ζ-model's zeros in 0.3 ≤ σ ≤ 0.7 are all on the line for X ≤ 100 at every height. Off-line zeros appear only near the strip edges, from X ≈ 30–50, and fewer at larger heights (X = 30: 16 → 8 → 0 → 0). This matches Gonek's "finitely many".
2. **Rigidity versus real-rootedness.** The largest X with no off-line zero anywhere in the strip is 20, 20, 30, 30 at T = 10³, 10⁴, 5·10⁴, 7·10⁴.
   - There the model's gap variance is 0.177, 0.133, 0.117, 0.111, against ζ's 0.153, 0.168, 0.158, 0.167 (statistical error about ±0.015).
   - The X at which the model first reaches ζ's variance moves up with height (between 10–20, 30–50, 50–100, 100–300). The onset of off-line zeros moves up more slowly (between 20–30, 20–30, 30–50, 30–50).
   - So from T ≈ 10⁴ on, the real-rooted models are measurably more rigid than ζ; at T = 10³ they are not.
   - Heuristic, unproved: the real-rooted range is limited by the Euler fluctuation against the conductor log(T/2π) and grows slowly with T. Matching ζ's local statistics needs much larger X. Lee–Yang-like real-rootedness is then a crystalline phase, with ζ-like statistics beyond it. The heights tested are small.
3. **DH test passed.** The DH model has off-line zeros inside 0.3 ≤ σ ≤ 0.7 from X = 10 at every height (8–14 per window, 2–5 % of zeros), not thinning with height. Their number matches the zeros of its Euler part c₁P_X(χ) + c₂P_X(χ̄) right of the line (5–7 model zeros against 6–8 Euler-part zeros per window, X = 10; not matched point by point). For ζ, P_X never vanishes; what DH lacks is a single Euler product.
4. No route to ζ: as X → ∞ the model is the wrong AFE (|P_X(1−s)| ~ exp(X^{σ−1/2}/log X)). Taking X → ∞ inside the real-rooted phase is RH (Gonek's theorem: approximation by P_X near the line ⟺ no zeros there).

**Verdict:** dies; known in kind (Gonek 2007). New here: the rigidity-versus-real-rootedness numbers and the DH contrast.

## T-C (N13, N7): Lee–Yang / Hermite–Biehler threshold of differenced ξ (`t3_lagarias.py`, `t3_out.txt`, `t3b_bisect.py`, `t3b_out.txt`)
**Exact statement** (`derivations.md` D3; ⇐ by pairing ρ with 1−ρ̄, ⇒ by the HB theorem): A_h and B_h have only real, simple, interlacing critical-line zeros ⟺ ξ ≠ 0 on Re s > 1/2 + h. So the Lee–Yang threshold equals the N7 half-plane constant (R47-type), and RH ⟺ threshold 0.
- The de Branges spaces H(E_h), E_h = ξ(s+h), give self-adjoint operators (canonical systems) whose spectra are the zeros of A_h, for every h above the threshold. This is a Hilbert–Pólya ladder whose missing rung is N7.

**Routine corollary.** For h > B₂ − 1/2 = 0.37338196661320752 (Bebel's Lean-checked half-plane; 3/8 from 7/8), A_h and B_h have only simple, interlacing critical-line zeros. Lagarias proved this for |h| ≥ 1/2.
- His "trivial spacing" theorem very likely extends unconditionally to h > 0.3734, since ζ'/ζ is o(log t) in zero-free half-planes by Borel–Carathéodory; this is unchecked.
- Light novelty check: arXiv search ("differenced xi", "xi(s+h)+xi(s−h)") and the 2026-10/11 listings of math.CA, math.NT, math.CV and math.SP show no later improvement of |h| ≥ 1/2; older literature is not scanned. Candidate for COROLLARIES after a reviewer verdict, not for DISCOVERIES.

**DH test (numerical, true DH function, off-line zero ρ₀ = 0.80851718246 + 85.69934848538 i, Re ρ₀ − 1/2 = 0.30852).**
- A_h has a zero right of the line in [80, 92] for h ≤ 0.29805 and none for h ≥ 0.29812.
- So the real-rootedness threshold sits 0.0104 below Re ρ₀ − 1/2: the region {|E_h| < |E_h^#|} around ρ₀ − h attaches to the line.
- Control: ξ gives none for h = 0.05, 0.15, 0.25.

The mechanism transfers to DH with DH's own half-plane, so on its own it is DH-blind; its only arithmetic input is the zero-free half-plane.

**Verdict:** circular (≡ N7). One routine corollary.

## Ranking (RH relevance × chance × gain)
| rank | transfer | node | relevance | chance | gain | result |
|---|---|---|---|---|---|---|
| 1 | T-C: HB / LY threshold of A_h, B_h | N13, N7 | high (exact RH-equivalent family) | 0 (≡ zero-free half-plane) | routine: h > 0.37338 (was 1/2) | circular; corollary candidate |
| 2 | T-B: Euler-product torus models with the Γ-factor | N13, N12 | medium | ≈ 0 (real-rooted only in the crystalline phase; X → ∞ is the wrong AFE) | DH-filter illustration | dies; known in kind (Gonek) |
| 3 | T-A: FQ hypotheses audit + prime-power rigidity lemma | N13, N1 | medium | 0 | NOGO line | dies |

## Suggested NOGO line (orchestrator)
- explore5 (Orca exploration: Fourier quasicrystals / Lee–Yang) NO IDEA. ζ's zero measure fails every FQ hypothesis: density log(T/2π)/2π, Fourier masses ~ e^{X/2} cancelled only by the pole (that is RH), a continuous Γ-part, and a prime-power spectrum, which on any finite prime torus forces arithmetic progressions (lemma). Euler-product torus models with the Γ-factor are real-rooted only for small X (from T ≈ 10⁴ on with zeros measurably more rigid than ζ's), while the DH model is off-line from X = 10 at every height. Lagarias' differenced ξ has a Lee–Yang threshold equal to the zero-free half-plane (now h > 0.37338; DH numerically 0.2981 against 0.3085). [ExploreQuasicrystal](Riemann/Deep/ExploreQuasicrystal/REPORT.md)

## Files
| file | content | run time |
|---|---|---|
| `derivations.md` | lemma D1, toy D2, threshold D3 | — |
| `t1_audit.py`, `t1_out.txt` | hypotheses audit | 14 s |
| `t2_models.py`, `t2_out.txt` | ζ and DH torus models; needs Odlyzko's `zeros1` (re-download to `.scratch/explore5/dl/`) | 1034 s |
| `t2b_euler_part.py`, `t2b_out.txt` | zeros of the DH Euler part (imports t2, so it also needs `zeros1`) | 400 s |
| `t3_lagarias.py`, `t3_out.txt` | A_h counts for DH and ξ; numpy Hurwitz–Euler–Maclaurin checked against mpmath to 5·10⁻¹² | 79 s |
| `t3b_bisect.py`, `t3b_out.txt` | DH threshold bisection | 217 s |

All runs were single-threaded at nice 19 under a hard alarm (≤ 1700 s), one process at a time. An earlier mpmath-only t3 hit its alarm and was replaced.

RH relevance: N13 (with N7; N1 and N12 touched)
