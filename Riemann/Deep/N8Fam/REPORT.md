CONDITIONAL RECORD: Johns's weighted family of primitive Dirichlet L-functions, zeros in (T, 2T] with T ≤ Q^κ: at least 79.784% (κ = 3), 76.528% (κ = 5), 72.879% (κ = 10) simple and on the critical line (Johns: 79.721%, 76.414%, 72.748%), and at least 89.963%, 88.402%, 86.624% distinct (Johns: 89.860%, 88.207%, 86.374%); for one prime modulus (Hua–Yang) and for each fixed χ, 67.353% and 83.900% (before 67.250%, 83.625%). Conditional on a routine, unwritten transfer lemma (§4) and on Johns's unrefereed (Lean-checked) energy bound. At polylogarithmic height no gain: 92.645% < Johns's 93.228% (node N8)

# N8Fam — transferring arXiv 2610.08965's certificates to families of Dirichlet L-functions

Run 2026-11-05, 12:02–13:40 UTC, one single-threaded process at a time. Agent n8fam.

## Short answer
1. **The brief's baseline is out of date.** Conrey–Iwaniec–Soundararajan's 56% is not the family record. Unconditional
   family results: Sono (arXiv 2105.07422, 2021) 61.07% on the line, 60.44% simple and on the line; Alpöge–Furman
   (arXiv 2608.13637v2, Remark 7.2, TeX `rem:Davg`) 0.811, announced without details; Hua–Yang (arXiv 2608.16034,
   2609.27808) 0.6725 for the characters modulo one prime; **Johns (Zenodo 23070760 v1.0.1, 30 Sep 2026; headline
   theorems Lean-checked, unrefereed): 0.9322 − ε at polylogarithmic height (Thm 1.1) and p(β̄_T; F₁) − ε at heights
   T ≤ Q^κ (Thm 1.4(a)), e.g. 0.8657 for T ≤ Q.** Johns proves an unconditional *upper* bound for the family energy with
   the GUE form factor F₁(α) = min(α, 1) for |α| < 2, by a sharp Toeplitz large sieve. That is the strongest
   unconditional family pair-correlation input found.
2. **Every input of 2610.08965 transfers** (table §1). Its proof uses the energy only as an upper bound,
   tr A² ≤ (C(f) + o(1))N, and Johns supplies exactly that, character by character, summed with weights ω_χ ≥ 0.
3. **Dilation lemma (§3, exact; reviewer: sound).** If (x, b) is a critical-line certificate for f with separation θ,
   then (x, λb) is one for f_λ(u) = f(u/λ)/λ with separation θ/λ, same Λ and β. So the paper's certificate
   (x = 4361/500000, b = 7·627·10⁻⁶) gives at bandwidth λ: s(λ) = (2 − C_fam(f_λ) − λb)/(1 − x). The distinct-zero
   certificate (first window, c = 67/20) scales the same way.
4. **Result.** The spread charge λb grows with λ, the reward x does not. Near bandwidth 2 (polylogarithmic height) the
   scaled certificate loses: 0.926455 against Johns's 0.932282. For bandwidth ≤ 5/4 (heights from Q³ on) it beats
   Johns's optimal plain constants by 0.0006–0.0013 (simple, on the line) and 0.0010–0.0025 (distinct), tending to
   0.673538 − 0.672501 and 0.839005 − 0.836250 as κ → ∞.
5. **First input that fails at polylogarithmic height:** the local certificate. A record there needs a new 8-point
   certificate for a window near Johns's optimiser; trial states cap it at ≤ 0.93307 (gain ≤ 0.0008), §5.

## 1. Inputs: zeta versus family
| Input of 2610.08965 (TeX line) | ζ | Family (Johns: primitive χ mod q, ηQ ≤ q ≤ Q, weight ω_χ = w(q/Q)q/φ(q)) | Transfers? |
|---|---|---|---|
| Energy (Prop. `prop:energy`, l.780) | BGST Lemma 5 (arXiv 2306.04799): tr A² = (C(f)+o(1))N for supp f*f ⊂ [−1,1]; form factor δ + \|α\| | Johns Prop. `prop:second` / `H:prop:second` with Prop. `prop:TIsharp` / `H:prop:TIsharp`: 𝔐 = Σ ω_χ‖Ĝ_χ‖² ≤ (1−2θ_J)Q_F(f_v)N + o(N), Q_F(f) = ∫f² + ∬f(x)f(y)F(\|x−y\|), F → min(α,1); windows on [−λ/2, λ/2], λ < β̄_T = log(Q²T)/log(QT) ∈ (1, 2]; zero density ℓ*/2π, ℓ* = log(QT) | yes: only "≤" is used (proof of Thm `thm:line-general`, l.1103) |
| Zero symmetry ρ ↦ 1 − ρ̄, same multiplicity (l.766) | functional equation | per χ: functional equation + L(s, χ̄) = conj L(s̄, χ) (Johns Lemma `lem:blocks`) | yes |
| Threshold lemma, level c = 2 (l.517, l.1072) | linear algebra | same, per χ, applied to Â_χ = P_χ + R_χ with n₊(R_χ) ≤ s₂ + p (Johns Lemma `lem:blocks`) | yes |
| Pinching, separated points (Λ < 2), close pairs (l.537–589) | window only | window only; dilation-invariant (§3) | yes |
| Local certificate (x, b), 8-point frames (l.1120–1236) | 535,332,163-box search (Lean, native_decide) | dilated window: (x, λb, θ/λ) | yes (§3) |
| Spread ≤ N + o(N) (l.1103) | heights in an interval of length TL/2π | spn(Y_χ) = (1+o(1))N_χ, since log q = log Q + O(log 1/η); Σ ω_χ spn ≤ (1+o(1))N | yes |
| Counting N, N₀ˢ | N(T) ~ (T/2π) log T | N = Σ ω_χ N_χ = (1+o(1))HTℓ*/2π (Johns Lemma `lem:RvM`); N₀ˢ = Σ ω_χ N₀ˢ_χ | yes (weights ≥ 0) |
| Gram matrix of simple on-line zeros = (K(y_i − y_j)) | exact (Lamzouri's operator) | Gabor form: Poisson-exact for interior zeros; a finite Gabor system is PSD-dominated by the full one, so U_SS ⪯ ΛI survives | routine (to write, §4) |
| Exterior zeros, bad characters | none | Johns §`sec:zero`: trace-norm tails; Montgomery's zero-density theorem drops characters of weight o(1) | yes (Johns) |

Families named in the brief:
- **q ~ Q average:** Johns's weighted family. The unweighted dyadic family has only the classical large sieve
  (Johns Cor. `cor:dyadic`, 0.7224); not pursued.
- **Single modulus q → ∞:** Hua–Yang for prime q: bandwidth 1 at polylogarithmic and sub-polynomial heights,
  1/(1+θ) at T = q^θ; form factor \|α\|, so C_fam = C for the paper's window at λ = 1.
- **Quadratic characters:** no unconditional family energy theorem found (arXiv, Zenodo). First missing input: the energy bound.
- **Fixed χ, T → ∞:** Alpöge–Furman Theorem B (TeX `thm:L`), bandwidth 1, the same C(f) as ζ.

## 2. Computations (scripts and logs in this directory)
Checks: C(f) of the paper's second window 1.32794727284456 (mpmath; paper 1.32794727284); first window
1.33270495828 (paper same); flat 4/3; Montgomery–Taylor 1.327499296321; flat at λ = 2 with F₁: 13/12; s = 0.6735383 and
0.8390046 (paper 1669159/2478195, 1645064/1960733); our 400-cell QP reproduces Johns's table of p(β̄; F₁) to 10⁻⁶.

**C_fam(f_λ) = Q_{F₁}(f_λ), paper's second window:** λ = 1: 1.327947; 12/11: 1.272769; 7/6: 1.236269; 5/4: 1.203633;
3/2: 1.136994; 2: 1.072848. For λ ≤ 1 it is the ζ value (support of f*f inside [−1, 1]).

**Simple zeros on the line, Johns's family, T ≤ Q^κ** (β̄ = (2+κ)/(1+κ); certificate = dilated paper window, λ = β̄):
| κ | β̄ | Johns p(β̄;F₁) | paper window, plain | with certificate | gain |
|---|---|---|---|---|---|
| 0 (polylog) | 2 | 0.932282 | 0.927152 | 0.926455 | −0.005828 |
| 1/2 | 5/3 | 0.894985 | 0.891233 | 0.891696 | −0.003290 |
| 1 | 3/2 | 0.865673 | 0.863006 | 0.863958 | −0.001717 |
| 2 | 4/3 | 0.824355 | 0.822935 | 0.824273 | −0.000083 |
| 3 | 5/4 | 0.797213 | 0.796367 | **0.797840** | +0.000626 |
| 5 | 7/6 | 0.764149 | 0.763731 | **0.765285** | +0.001135 |
| 10 | 12/11 | 0.727484 | 0.727231 | **0.728799** | +0.001314 |
| 20 | 22/21 | 0.703011 (QP) | 0.702722 | **0.704266** | +0.001255 |
| 100 | 102/101 | 0.679185 (QP) | 0.678781 | **0.680282** | +0.001096 |
| ∞ | 1 | 0.672501 | 0.672053 | **0.673538** | +0.001038 |
The κ = 3, 5, 10 values were recomputed with mpmath (30 digits): 0.7978400007, 0.7652852362, 0.7287992526.
The crossover lies between κ = 2 and 3 (β̄ ≈ 1.30); between grid points the claim is numerical only.

**Distinct zeros** (first window, c = 67/20, (x₁, x₂, b) = (0.054350, 0.069433, 6·0.006219), dilated):
N_d/N ≥ (3 − C(f_λ) + x₂ − x₁ − λb)/(2 − 2x₁ + x₂); μ, ν, r_h, r_k, Λ < c/2 and 2 + 2√β < c do not depend on λ.
| κ | Johns (1+p)/2 | with certificate | gain |
|---|---|---|---|
| 0 | 0.966141 | 0.953683 | −0.012459 |
| 1 | 0.932836 | 0.929950 | −0.002887 |
| 2 | 0.912177 | 0.912050 | −0.000128 |
| 3 | 0.898606 | **0.899634** | +0.001027 |
| 5 | 0.882074 | **0.884022** | +0.001947 |
| 10 | 0.863742 | **0.866246** | +0.002503 |
| ∞ | 0.836250 | **0.839005** | +0.002755 |

**Single prime modulus** (Hua–Yang 2609.27808: bandwidth 1/(1+θ) at T = q^θ; form factor α):
| θ | λ | plain optimum (= Hua–Yang) | with certificate | gain |
|---|---|---|---|---|
| → 0 | 1 | 0.672501 | **0.673538** | +0.001038 |
| 0.05 | 0.9524 | 0.637556 | **0.638208** | +0.000652 |
| 0.1 | 0.9091 | 0.601315 | **0.601465** | +0.000150 |
| 0.25 | 0.8 | 0.486267 | 0.484450 | −0.001817 |

**Weaker family inputs** (paper window, best λ; for comparison): Alpöge–Furman-type large sieve, form factor α on [0,1]
and 1.263α on (1,2): 0.822037 plain, 0.822366 with certificate (their announced 0.811 must use a weaker functional);
GUE value on (1, 1/0.87338) from the 0.87338 zero-free half-plane, 1.263α beyond: 0.8410. Johns's F₁ supersedes both.

## 3. Dilation lemma (exact)
Let f be an even probability density on I = [−½, ½], K = f̂, and let (x, b) be a critical-line certificate built as in
2610.08965 (local data for eight points with charge α, b = 7α; separation θ; majorant g ≥ f with ĝ supported in (−θ, θ),
Λ = ∫g < 2; K² ≥ β ≥ x on [0, θ]). For λ > 0 put f_λ(u) = f(u/λ)/λ, so K_λ(t) = K(λt). Local data for K_λ with
separation θ/λ: the same W_ij and κ_i, charge λα, storage functions B_r(λ·); substitute y' = λy in (eq:local8):
λα(y₇ − y₀) = α(y'₇ − y'₀), K_λ(y_j − y_i) = K(y'_j − y'_i), gaps λg ≥ θ. The majorant g_λ(u) = g(u/λ)/λ ≥ f_λ has
ĝ_λ(t) = ĝ(λt), supported in (−θ/λ, θ/λ), and ∫g_λ = Λ; K_λ² ≥ β on [0, θ/λ]. So (x, λb) is a certificate for f_λ,
and Thm `thm:line-general` gives s = (2 − C(f_λ) − λb)/(1 − x). Labelled data (first window) scale the same way.

## 4. What remains to be proved
**Transfer lemma (routine; to write).** Johns's family F(Q, w), w ∈ {w♯, w_smooth}, η → 0; fix κ and
λ₀ < β̄(κ). Let f be an even probability density on I, √f smooth inside I, with a critical-line certificate (x, b) as in §3.
Then uniformly for ℓ^{a₀} ≤ T ≤ Q^κ,
liminf_{Q→∞} N₀ˢ/N ≥ (2 − Q_{F₁}(f_{λ₀}) − λ₀b)/(1 − x) − ε, and the values of §2 follow as λ₀ ↑ β̄(κ).
A single λ₀ serves every T ≤ Q^κ because β̄_T ≥ β̄(κ) > λ₀; no monotonicity in λ is needed. Steps:
(i) in Johns Prop. `prop:perchi`, replace the rank–trace inequality by the threshold lemma at c = 2:
‖Â_χ‖² ≥ 4 tr Â_χ − 4 tr P_χ + tr φ₂(P_χ) − 4(s₂ + p) ≥ 2(N_χ − s₁) + tr φ₂(P_χ) − (trace defect), using
N_χ − s₁ ≥ 2(s₂ + p); (ii) bound the trace defect N_χ − tr Â_χ on weighted average by o(N) (Johns Prop. `prop:trace` and
the exterior-zero tails); (iii) the normalised Gabor Gram matrix of the simple on-line zeros equals (K̃(y_i − y_j)) for
the actual window ṽ = f_ε (the paper's own smoothing f χ_ε²/m_ε is admissible for Johns: its square root is smooth
with compact support inside I) by Poisson summation, up to an explicit edge error; zeros within the edge layer go to
R_χ (o(N) of them); finite truncation lowers the matrix in PSD order, so U_SS ⪯ (Λ/m_ε)I < 2I; (iv) sum over χ with
ω_χ ≥ 0 and use Σ ω_χ spn(Y_χ) ≤ (1 + o(1))N; (v) Johns Prop. `H:prop:second` with Prop. `H:prop:TIsharp`.
All per-character constants (7x + Ω, edge terms) are O(1) or O(√T ℓ) per character, so o(N) after weighting.
The same lemma with Hua–Yang's (prime modulus) or Alpöge–Furman's Theorem B (fixed χ) energy gives the 0.673538 rows,
and with the labelled Gram bound (Thm `thm:gram`) the distinct-zero rows.
**For a record at polylogarithmic height (open):** a certificate (x, b) for a window f near Johns's optimiser with
(2 − Q_{F₁}(f₂) − 2b)/(1 − x) > 0.932282: new kernel table, majorant, LP and box search (the paper's took 535·10⁶ boxes).

## 5. Ceilings from trial states (ceil_fam.py; upper bounds for any certificate that sees pairs at index distance ≤ s)
| Window | plain | s = 7 (8-point frames) | s = 40 (all pairs) |
|---|---|---|---|
| paper window, λ = 1 (check; paper: 0.6735660 / 0.6738216) | 0.672053 | ≤ 0.673616 | ≤ 0.673854 |
| paper window, λ = 1.999 | 0.927072 | ≤ 0.927653 | ≤ 0.927710 |
| Johns's optimum, β = 2 | 0.932282 | ≤ 0.933069 | ≤ 0.933354 |
| Johns's optimum, β = 3/2 | 0.865674 | ≤ 0.867395 | ≤ 0.867935 |
| Johns's optimum, β = 4/3 | 0.824356 | ≤ 0.826060 | ≤ 0.826583 |
Our state library is slightly weaker than the paper's (check row: 5·10⁻⁵ above its proved ceilings). So with the paper's
window no certificate can reach Johns's 0.932282 at polylogarithmic height, and with Johns's optimal windows a new
8-point certificate could add at most 0.0008 there and 0.0017 at T ≤ Q or T ≤ Q².

## 6. Novelty pass
- arXiv math.NT listings 2026-11 (140 entries to 5 Nov) and 2026-10 (854), scanned by title; abstracts read for every
  Dirichlet/zero/family hit (2610.28527, 2610.29018, 2610.31744, 2610.26841, 2610.41394, 2610.23467, 2610.04045,
  2610.24612, 2610.05227): none transfers 2610.08965's certificates to a family or to L(s, χ).
- arXiv search (arXiv mode): CLLR 1211.6725 (GRH, 11/12; GRH used in its Lemma `lem:sumnonprincChar` for the
  off-diagonal), Wu 1206.1679, Sono 2105.07422, Dickinson 2211.06264, Kandhil–Languasco–Moree 2411.19762, BGST 2306.04799
  and 2501.14545, Alpöge–Furman 2608.13637v2, Lamzouri 2609.02882, Hua–Yang 2608.16034 and 2609.27808, Ji 2605.09282.
- Zenodo (title queries Dirichlet ∧ zeros / critical line / simple / family): Johns 23070760 (only version; GitHub
  nic410/dirichlet-critical-zeros has README-only commits after 6 Oct), Prüzelius 22071915 (prime modulus, 2/3),
  Fabbian 23025218 (weight-2 newforms); none uses 2610.08965.
- 2610.08965 mentions no L-function other than ζ (rg over its TeX).
- Ours: NOGO.md, TREE.md, COROLLARIES.md, NOVELTY files: nothing on Johns, Hua–Yang or a family transfer; TREE N8 and
  R318 still quote CIS 56% (out of date).
- Review: reviewer verdict on the dilation lemma: sound. On the transfer: correct in outline, but the threshold step for
  the Gabor compression, the edge and trace errors, and the uniformity in T must be written out (now in §4); the
  improvement is claimed for the listed κ only.
- Rounding: 2610.08965's constant 1669159/2478195 = 0.6735382 is "at least 67.353%" when rounded down; TREE.md and
  the brief write 67.354%.
- Verdict: new as far as found, but a routine combination with gains of order 10⁻³, resting on an unrefereed input.
  Classification: COROLLARIES.md, not DISCOVERIES.

## 7. Files
- `cfam.py` — C_fam(f_λ) for the paper's, flat and Montgomery–Taylor windows under form-factor models (ζ, F₁,
  Alpöge–Furman-type, zero-free-half-plane), Gauss–Legendre, accurate to 10⁻¹²; ~2 min; `logs/cfam.log`.
- `pbeta.py` — Johns's p(β̄; F₁) by a 400-cell QP and the dilated certificate per κ; ~1 min; `logs/pbeta.log`.
- `pbeta_small.py` — bandwidth ≤ 1 (single prime modulus, fixed χ); ~40 s; `logs/pbeta_small.log`.
- `distinct.py` — distinct zeros with the first window; ~3.5 min; `logs/distinct.log`.
- `ceil_fam.py` — trial-state ceilings (args: zeta | fam | mid); < 1 min each; `logs/ceil_*.log`.
- `check_mp.py` — mpmath recomputation of the κ = 3, 5, 10 values; ~7 min; `logs/check_mp.log`.
Run: `OMP_NUM_THREADS=1 uv run --no-project --with numpy --with scipy python3 <script>` (check_mp.py: `--with mpmath`).

RH relevance: N8
