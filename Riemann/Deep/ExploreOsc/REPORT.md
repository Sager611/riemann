BARRIER: Eremenko–Novikov and Sturm–Hurwitz–Kellogg need an exact spectral gap, but Z (like Davenport–Heilbronn, χ₈, χ₋₈) has a flat local spectral measure reaching frequency 0, and any argument whose only input on Z is that measure (all untwisted mean values of Z^{(i)}(t)Z^{(j)}(t+h)) proves at most 1/√3 = 57.735% of zeros on the line, below the 67.36% record, since a variance-matched Gaussian process has the same measure and exactly that zero density (Rice) (node N8)

# explore60 — oscillation theory (spectral gaps force sign changes) → N8 (Orca exploration, 20 Nov 2026)

## 0. Repo check
- `rg -i 'Sturm|Hurwitz|Eremenko|Novikov|Poltoratski|Mitkovski|Kellogg|zero.crossing|Kedem|Logan|sign change|spectral gap|Beurling.Malliavin|Kac.Rice'`
  over NOGO.md, Riemann/Deep/*/REPORT.md: Eremenko–Novikov and Mitkovski–Poltoratski only as pointers in GapJInf §Novelty
  (reverse question); Sturm only as nodal counts (ExploreGFF) and Prüfer angles (ExploreBS: Prüfer angle of Z = Levinson–
  Montgomery's arg ζ′ count); Kac–Rice in R343 ("random-model counts do not transfer"). No hit: Kedem, Logan, higher-order
  crossings, Rice's formula as a method ceiling.
- Repeats avoided: explore37/41/42/58/59, GapJInf, R343, R34/R56 (Laguerre–Pólya toolbox), explore57 (Prüfer = Levinson).

## 1. Transfer
Tools, exact statements:
- (b) Eremenko–Novikov (math/0301060, abstract): if f̂ = 0 on [−a, a], f a locally integrable temperate distribution, the
  lower density of sign changes of f is at least a/π.
- (a) Sturm–Hurwitz–Kellogg, finite form: if f is orthogonal on I (positive weight) to a k-dimensional Chebyshev system,
  f has ≥ k sign changes in I; equivalently, f has ≥ k sign changes in I iff no nonzero P of the system has fP ≥ 0 on I.
- (d) Kedem: the counts D_k of sign changes of k-th differences are nondecreasing in k and (Gaussian case) determine the
  spectrum; Logan: a band-pass signal of less than one octave is fixed by its zero crossings.
- (c) Mitkovski–Poltoratski: G_X = 2πD_*(X) depends only on the Beurling–Malliavin density of X.
Node N8. Sign changes of Z ≤ zeros on the line ≤ N(T); N(2T) − N(T) = (θ(2T) − θ(T))/π + O(log T), so the maximal
sign-change density is θ′/π, and the record (≥ 67.36% of all zeros simple and on the line) is ≥ 0.6736·θ′/π sign changes.
- (C1) E–N transferred: if W ≥ 0 and (ZW)^ = 0 on (−a, a), then N₀(T) ≥ (a/π)T(1 − o(1)) (W ≥ 0 keeps Z's sign changes).
  A fixed a gives ≍ T, i.e. 0% of N(T) ~ (T/2π) log T; a proportion c needs a local gap (−cθ′(t), cθ′(t)).
- (C2) Z has no gap. (i) Original variable: with φ̂ supported in (−a, a), φ_X(t) = φ(t/X)/X, X = √t_n, t_n = 2πn²,
  stationary phase at the n-th Riemann–Siegel term (phase −πn² − π/8 + (t − t_n)²/(4t_n)) gives
  (Z * φ_X)(t_n) = (−1)^n 2n^{−1/2} Re(e^{−iπ/8}∫e^{iw²/4}φ(w)dw)(1 + o(1)), nonzero for generic φ, so Ẑ ≠ 0 on (−a, a)
  for every a > 0 (standard, not written out; the stationary-zone means in 2d are ±(0.1–1.7)t_n^{−1/4} for ζ); for χ₈,
  ∫_0^T Z ~ CT^{3/4} (2d). So the hypotheses of (a) and (b) (exact gap, exact orthogonality) fail for Z, Z_χ, Z_DH.
  (ii) Unfolded variable: F_Z(v) = Z(t)/√(2θ′(t)), θ(t) = πv. The diagonal of Z = 2Σ_{n≤N} n^{−1/2}cos(θ − t log n) +
  O(t^{−1/4}) puts mass 1/(2θ′n) at ±π(1 − log n/θ′): the unfolded spectral measure → uniform on [−π, π] weakly, error
  O(1/log t), the n = 1 (Gram) atom at the top carrying ≈ 1/(θ′ + γ) (30% at t = 1500). Its Rice ratio, from
  Σ_{n≤N} log^k n/n = log^{k+1}N/(k+1) + γ_k + o(1): ρ(T)² = Σ n^{−1}(θ′ − log n)²/(θ′²Σ n^{−1}) =
  (θ′³/3 + γθ′² − 2γ₁θ′ + γ₂)/(θ′²(θ′ + γ)) = 1/3 + 2γ/(3θ′) + O(θ′^{−2}): .682 at T = 1500 (time domain .681, acf.py),
  .6290 at 10⁶ (rmodel .6290), .601 at 3·10¹², → 1/√3. Same picture for every primitive Z_χ and, |a_n|² averaged, DH.
- (C3) Ceiling. 𝒮 = {f : F_f(v) = f(t)/√(2θ′(t)) satisfies (1/V)∫_V^{2V} F_f^{(i)}(v)F_f^{(j)}(v + u) dv →
  c_f(−1)^i (d/du)^{i+j} sinc(u), sinc(u) = sin(πu)/(πu), for all u and i, j ≤ J, some c_f > 0}. Z ∈ 𝒮 by the diagonal
  computation (ii), granted the off-diagonal and differentiated Riemann–Siegel error estimates (standard Montgomery–Vaughan
  type, not proved here; numerics consistent). Claim: a theorem "every f ∈ 𝒮 has ≥ c(θ(2T) − θ(T))/π real zeros in
  [T, 2T] for large T" needs c ≤ 1/√3. Proof: Y stationary Gaussian, spectral density uniform on [−π, π], covariance
  sinc: absolutely continuous spectrum, so Y and (Y, Y′, …, Y^{(J)}) are mixing (Maruyama), hence ergodic, and
  E[Y^{(i)}(v)Y^{(j)}(v + u)] = (−1)^i sinc^{(i+j)}(u); so X(t) = √(2θ′(t))·Y(θ(t)/π) has F_X = Y and lies in 𝒮 a.s.
  Its real zeros are a.s. simple, and those in [T, 2T] are Y's in [θ(T)/π, θ(2T)/π]: by Rice (λ₂ = π²/3 < ∞) and the
  ergodic theorem, a.s. ~ (1/π)√(λ₂/λ₀)·(θ(2T) − θ(T))/π = (1/√3)(θ(2T) − θ(T))/π. Exact, no sampling. So no argument
  whose only input on Z is membership of 𝒮 (local spectral measure, all untwisted mean values ∫ Z^{(i)}(t)Z^{(j)}(t + h) dt)
  proves more than 57.735% of zeros on the line, however counted; the record (simple on-line zeros, which are sign
  changes) needs inputs outside 𝒮 (twisted moments, argument principle).
- Evidence the true infimum over 𝒮 is lower (not proved; 2e): amplitude-modulated chirps have 9.7% of the maximal density
  with flat windowed spectra; repeated periodically they have line spectra with a flat envelope, so they lie in 𝒮 only
  after smoothing (e.g. random periods), which is not checked.

## 2. Tests (scripts in this folder, outputs in logs/)
Functions: ζ; Davenport–Heilbronn f (mod 5, odd); χ₈ (mod 8, even, S(χ₈) = 8e(−1/16) ≠ 0: Z-biased, 2610.31744); χ₋₈ (odd,
S = 0). oscz.py: Euler–Maclaurin Hurwitz ζ with 12 Bernoulli terms; Z against mpmath ≤ 7.9·10⁻¹³, discarded imaginary part
≤ 2·10⁻¹² of max|Z| (logs/zgrid.log). Grid t ∈ [0, 2000], step 0.01. Sign changes to 2000: ζ 1517 (= N(2000)), χ₈ and χ₋₈
2178 (= N), DH 1861 against ≈ 2029 zeros (168 off the line). Peaks: 226 MB (zgrid, mpmath checks), 70 MB (rmodel), else < 30 MB.

### 2a. Smallest case: first 100 zeros (t ≤ 236.6)
- Power fraction below aθ′, a = .1/.25/.5/.75/1, at t = 150 (window 60, Hann): ζ .058/.152/.311/.530/.778; DH
  .075/.143/.350/.390/.689; χ₈ .121/.211/.351/.513/.753; χ₋₈ .021/.125/.281/.461/.727 (flat = a; top deficit = leakage).
  No gap at any a: (C1) gives nothing.
- Rice density of the windowed spectrum /(θ′/π): .750, .813, .751, .789 (window-inflated; exact flat value .577); actual
  zero density /(θ′/π): .924, 1.007, .997, 1.037.
- Kedem D_0..D_5 on [20, 236.6]: ζ 99, 98, 99, 98, 99, 98 (θ-count 98.5); DH 146, then 154 for k ≥ 1 (θ-count 154.0;
  8 zeros off the line); χ₈ 170/171; χ₋₈ 170 (θ-count 170.2).

### 2b. Bounded scale (t ≤ 2000: 1517, ≈ 2029, 2178, 2178 zeros)
- F(a) at t = 1800 (window 80): ζ .065/.235/.438/.637/.868; DH .083/.146/.334/.515/.780; χ₈ .141/.191/.382/.611/.818;
  χ₋₈ .039/.206/.400/.621/.823. Window 200, a = .25/.5/.75 (analyze2.log): Z .17–.21/.33–.41/.49–.68 for all four.
- sign(Z) = ZW with W = 1/|Z| (the ideal E–N multiplier): .056–.077/.115–.185/.186–.251 for all four: damped, not gapped.
- Z|M|², M = Σ_{m≤30} μ(m)(1 − log m/log 30)m^{−s} (ζ's mollifier): ζ .13–.14/.19–.23/.24–.28, the others
  .14–.22/.25–.41/.48–.66. Only the arithmetic mollifier, and only for ζ, moves power up: Selberg's mechanism, outside 𝒮.
- Unfolded autocorrelation on [1000, 2000] (acf.py, logs/acf.log) at u = .25/.5/.75/1/2, sinc = .900/.637/.300/0/0:
  ζ .863/.514/.112/−.174/.184, DH .834/.421/−.029/−.301/.344, χ₈ .844/.455/.020/−.261/.274 (top atom: excess at u = 2);
  time-domain Rice ratio ζ .681, DH .751, χ₈ .726, χ₋₈ .724; sign changes/θ-count ζ .9996, χ₈ .9997, DH .9008.
- Rice /(θ′/π) at t = 1000, 1800 (Hann windows): .665–.762 for all four; actual: ζ 1.022, 1.000; χ₈ 1.011, 1.015; χ₋₈ 1.011, 1.015;
  DH .847, .843 (15% of its zeros off the line there).
- Kedem on [20, 2000]: ζ D_0..D_5 = 1516 (constant; D_6 is rounding noise), χ₈ 2170–2171, χ₋₈ 2171 (θ-count 2170.9);
  DH D_0 = 1855, D_1..D_6 = 2023 (θ-count 2022.8). For DH the first derivative already reaches the θ-count; the jump
  D_1 − D_0 = 168 = θ-count − D_0, its off-line deficit up to S(t). For ζ the sequence is saturated at k = 0 (RH: Z′/Z decreases between zeros).

### 2c. Same spectrum, other phases (rmodel.py; Riemann–Siegel main sums; window 300 at T = 10⁴, 150 at T = 10⁶)
| phases (sign changes / θ-count) | ζ, 10⁴ | ζ, 10⁶ | DH, 10⁴ | DH, 10⁶ |
|---|---|---|---|---|
| true: θ(t) − t log n | .975 | .997 | .906 | .910 |
| random multiplicative φ_n = Σ v_p(n)φ_p (6 / 4 seeds) | .783 ± .026 | .748 ± .058 | .846 ± .016 | .818 ± .016 |
| independent φ_n | .700 ± .026 | .665 ± .015 | .782 ± .015 | .715 ± .021 |
| Rice, discrete spectrum (→ 1/√3) | .659 | .629 | .721 | .680 |
Each column has one spectrum; the count runs from ≈ Rice to ≈ 1. Independent prime phases (Bohr–Jessen without the θ
coupling) do not restore it; maximality needs the functional-equation phases. True-phase deficits at 10⁴: main-sum truncation.

### 2d. Frequency 0 (χ₈; 2610.31744 Thm GL1)
I(T) = ∫_0^T Z. χ₈: I/T^{3/4} = 1.2575, 1.2713, 1.2226, 1.1964, 1.2705 at T = 250, 500, 1000, 1500, 2000, against
C = 2^{3/2}π^{1/4}/3 = 1.255195; |I − CT^{3/4}| ≤ 0.453 T^{1/4} log T on [100, 2000]. ζ, DH, χ₋₈: |I| ≤ 6.71, 6.53, 3.30
times T^{1/4}. ζ's stationary-zone means alternate as (−1)^n (t_n = 2πn²); χ₈'s add up. The frequency-0 term is arithmetic
(S(χ)) but carries T^{−1/2}/log T of the power: invisible to every spectral-measure argument, and Hardy–Littlewood's
detector (|∫_I Z| < ∫_I |Z|) works the same for χ₈ (coherent part H T^{−1/4} ≪ H).


### 2e. Deterministic witnesses in 𝒮 with few zeros (witness.py, logs/witness.log; unfolded units, maximal density 1)
Amplitude-modulated chirp: instantaneous frequency w(u) spends time τ(w) ∝ 1/(w + δ) at w ∈ [0, π], amplitude √(w + δ),
so the power a²τ is flat; zeros are exactly those of cos(phase), density 1/log((π + δ)/δ) − δ/π.
| δ | F(.1/.25/.5/.75/1), flat = a | zeros per unit (predicted) |
|---|---|---|
| 10⁻² | .090/.261/.469/.736/.982 | .1738 (.1706) |
| 10⁻³ | .100/.268/.480/.729/.982 | .1242 (.1239) |
| 10⁻⁴ | .101/.260/.476/.725/.982 | .0966 (.0965) |
| ∞ (linear chirp) | .103/.256/.502/.740/.993 | .5008 (.4997) |
Windowed spectra are flat to ±.035 in F (Z itself only to ±.15, 2a–2b), but a periodic chirp has a line spectrum
(flat envelope), so this is evidence, not proof, that the infimum over 𝒮 is below 1/√3, possibly 0. Classical analogue:
Σ_{a∈A} cos(ax), |A| = n, can have as few as n^{5/6+o(1)} zeros (Borwein–Erdélyi–Ferguson–Lockhart 2008) and
C(n log n)^{2/3} (Juškevičius–Sahasrabudhe, 2005.01695). Certified ceiling: 1/√3 (Gaussian, Rice); numerical: 0.0966.

## 3. Adversarial test
- DH: same flat spectrum (2a, 2b), Kedem-saturated from k = 1, 8–15% of its zeros off the line for t ≤ 2000. Every
  statement proved for all of 𝒮 holds for DH and for the Gaussian X (a positive-definite model: covariance = Fourier
  transform of the flat measure). So (C1)–(C3) are not about ζ: that is the barrier.
- χ₈: same spectrum plus a coherent frequency-0 term (2d) of relative power T^{−1/2}/log T: no difference for (a)–(d).
- What does tell ζ apart in these data: its true phases (2c: .997 against .75 for random Euler-product phases) and its
  mollifier (2b: low band halved for ζ only). Both are arithmetic and outside 𝒮.

## 4. Where the transfer breaks
- E–N needs an exact gap; Z, Z_χ, Z_DH have none (flat down to 0; the frequency-0 content is forced by the stationary
  points t = 2πn²/q). A fixed gap gives ≍ T zeros = 0% of N(T).
- Robust forms are false: f = (1 + cos x)^{2K} ≥ 0 has no sign change, yet its power below frequency 1 (the DC term) is the
  fraction C(4K, 2K)²/C(8K, 4K) ~ 1/√(πK) → 0 of its total (exact: .2476, .0796, .0252 at K = 5, 50, 500); Z is as peaky (E|Z|/√(E Z²) → 0), so L¹ input is needed.
- Repairs leave 𝒮: E–N on Z|M|² with L¹ control (|∫_I ZW| < ∫_I |Z|W) is Selberg's method (k = 1); Kellogg's k-fold version
  needs Z ⊥ all frequencies below ≈ cθ′ on each window, which the flat spectrum forbids; the argument form of E–N (sign
  changes ≥ phase winding/π, zeros of the analytic signal add winding) is Levinson's method; Kedem's HOC = the Prüfer/
  Levinson–Montgomery count (explore57) and Bui–Hall's pairing of zeros of Z^{(k)} (2304.05178); Logan needs a band-pass
  signal (Z's band reaches 0; its top octave n ≤ N^{1/2} is a short Dirichlet polynomial whose zeros are not ζ's);
  Mitkovski–Poltoratski's G_X = 2πD_*(X) sees only the density of the zero set (equal for ζ, DH, χ after unfolding).
- Scope of the ceiling: inputs = Z's local spectral measure, i.e. all untwisted second-order mean values
  ∫ Z^{(i)}(t)Z^{(j)}(t + u) dt. Twisted moments (mollifiers), fourth moments and the argument principle are outside;
  every N8 record uses them.


## 5. Novelty
- arXiv API: "higher order crossings" AND zeta; Eremenko AND Novikov AND zeta; "spectral gap" AND "sign changes" (AND zeta);
  "Rice formula" AND zeta; "Kac-Rice" AND "Riemann zeta"; "zero crossings" AND Riemann; random AND "Hardy Z-function" AND
  zeros; "Hardy Z-function" AND derivatives AND zeros (hits: 2304.05178 Bui–Hall, paired zeros of Z^{(k)}; 1612.08627
  Lehmer pairs; neither states a spectral ceiling); Littlewood AND cosine AND zeros (2005.01695, 2407.16075, 1610.07680:
  flat-weight cosine sums). Listings 2026-10/11 in math.NT/CA/CV/SP/PR/FA/CO (regex sign change|spectral gap|zero.crossing|
  Eremenko|Rice formula|Kac.Rice|Sturm|higher.order crossing): 27 title hits, none on ζ. Zenodo (three queries): nothing.
- Classics behind each step: Rice (1944–45), Kac (1943), Dunnage (1966: random trigonometric polynomials of degree n have
  ≈ 2n/√3 zeros, the same 1/√3), Kedem (higher-order crossings), Logan (1977), Eremenko–Novikov (2003), BEFL (2008).
  Verdict: routine consequences of classical facts; new only as a stated ceiling for spectral-only N8 arguments.
- Review: two verdict-only reviewer passes on (C3). First (recap after its time limit): Rice/ergodicity and the Gaussian
  ceiling sound; fixes applied: precise class 𝒮 (normalised unfolded autocorrelation, σ² = 2θ′), E–N/Kellogg excluded by their
  hypotheses rather than by "reading only the measure", sign-change vs zero counts, chirps downgraded to evidence,
  robust-form counterexample made exact. Second: "after these repairs, the Gaussian gives the claimed 1/√3 universal
  sign-change ceiling"; its fixes applied: (−1)^i (not (−1)^j) in the mixed covariances, Z ∈ 𝒮 stated as granted (not
  proved), the no-gap claim for E–N argued in the original variable (stationary phase), simple ⇒ sign change only.

## Suggested NOGO line
- explore60 (Orca, exploration: oscillation theory: spectral gaps and sign changes — Eremenko–Novikov, Sturm–Hurwitz–Kellogg, Kedem's higher-order crossings, Logan, Mitkovski–Poltoratski → N8) BARRIER: Eremenko–Novikov and Sturm–Hurwitz–Kellogg need an exact spectral gap, but Z's local spectral measure is flat on [−θ′, θ′] down to frequency 0 (so are Davenport–Heilbronn's, χ₈'s, χ₋₈'s); a variance-matched Gaussian with that measure has exactly 1/√3 = 57.735% of the maximal zero density (Rice), so arguments using only that measure (all untwisted mean values of Z^{(i)}Z^{(j)}) cannot pass 57.735% < 67.36% (reviewer: ceiling sound after fixes; chirps with 9.7% are numerical evidence only); same spectrum with independent or random Euler-product phases gives .665 and .748 against ζ's .997 (T = 10⁶); Kedem's D_k is constant for ζ, χ₈, χ₋₈ and jumps at k = 1 by DH's 168 off-line zeros; χ₈'s frequency-0 term (2610.31744) confirmed (∫Z/T^{3/4} = 1.20–1.27 vs 1.2552); the repairs (mollifier with L¹ control, phase winding) are Selberg's and Levinson's methods [ExploreOsc](Riemann/Deep/ExploreOsc/REPORT.md)

RH relevance: N8 (critical-line proportion): an exact ceiling (1/√3) for every argument that uses only Z's local spectral
measure; unconditional for the class 𝒮, applied to Z through the standard diagonal mean value (granted, not re-proved);
it does not touch mollifier, pair-correlation or TEN methods.

## Files
- oscz.py (Z for ζ, DH, χ₈, χ₋₈), zgrid.py (grids; logs/zgrid.log), analyze.py (E1 spectra, E2 ∫Z, E3 multipliers, E4 Kedem;
  logs/analyze1.log, analyze2.log), rmodel.py (phases; logs/rmodel.log), witness.py (chirps; logs/witness.log), acf.py (unfolded
  autocorrelation, time-domain Rice ratio; logs/acf.log).
- Reproduce: `uv run --no-project --with mpmath --with numpy --with scipy python3 Riemann/Deep/ExploreOsc/zgrid.py 2000 0.01
  zeta,dh,chi8,chim8` (90 s; writes logs/grid_*.npz, deleted after use, 8 MB), then `analyze.py E1,E2,E3,E4`,
  `rmodel.py 10000 300 6`, `rmodel.py 1000000 150 4`, `witness.py <delta> 10 20`; every run under memcap 290 MB.

