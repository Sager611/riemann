NO RECORD: the paper's 8-point certificate is LP-saturated for its window (re-optimised data ≤ 0.673540 at the paper's margin; proved ceiling 0.673566); more needs a new window or longer frames, each a new 0.5–10-billion-box proof

# N8Opt — re-optimising arXiv 2610.08965 v1 (≥ 67.354% of zeta zeros simple and on the line)

Node N8 (critical-line proportion). Run 2026-11-05, 09:50–11:55 UTC, one single-threaded process at a time.

## Short answer
Re-tuning cannot visibly beat 67.354%.
- For the paper's window, a cutting-plane LP over all 121 local-data parameters gives at most **0.6735404** at the
  paper's verification margin (paper: 0.6735382, so +2·10⁻⁶ at best), and at most 0.6735490 even with zero margin,
  which no finite box search can certify.
- The paper itself proves 0.6735660 as the limit of every 8-point certificate for that window (Prop. 7.3, in Lean).
- A step of the window coefficients along the LP gradient gave no real gain (below); θ is not binding.
- More room exists only with longer frames (this window: ≤ 0.673599 with 9 points, ≤ 0.673667 with 12, ≤ 0.6738216 with all pairs) or another window (≤ 0.68185 for every
  window with V(f) ≤ 10B(f); ≈ 0.674–0.675 for the window families studied so far). Each needs a new kernel table, majorant check, LP and a box search of
  ≥ 0.5 billion boxes (≈ 65–90 min per 8-point run on this machine; 9-point frames need an 8-dimensional search).

## 1. The optimisation problem (task 1)
Fixed by the theory: A = Σ ν_z F_z⟨F_z̄,·⟩, tr A = N, tr A² = (C(f)+o(1))N (BGST, unconditional); threshold lemma at
level c = 2 (forced: an off-line pair or double zero is credited 4c − c² ≥ 4 only for c = 2); only simple zeros on the
line interact. A certificate (x, b) means tr φ₂(U) − n ≥ x n − b·spread(Y) − o(n) − O(1). Result (Thm 4.1):

  liminf N₀ˢ/N ≥ s = (2 − C(f) − b)/(1 − x).

Free parameters (paper's values) and constraints:
1. Window f(u) = 1 + Σ_{j≤J} a_j cos 2πju on [−½, ½], J = 10, a in units of 10⁻⁶:
   (122290, −45012, 16241, 1887, 2711, 3731, −5307, 2982, −5153, 5358). Need f ≥ 0 (1 − Σ|a_j| = 49333/62500).
   K(t) = sinc t + Σ (a_j/2)(sinc(t−j) + sinc(t+j)); C(f) = 4/3 + ½Σa_j² + π⁻²Σ((−1)^j a_j/j² − a_j²/4j²)
   = 1.32794727…, rational cap C₀' = 6639737/5000000.
   Rounding forced by the exact checks (C++/Lean): a_j, α, W, κ, knot values in 10⁻⁶; table entries in 2⁻²⁸ on cells
   of 1/16384 with π enclosed by 40-digit rationals; θ rounded down to the grid; γ₀ and C₀' rational.
2. Separation θ = 17/20 (grid 13926/16384) and majorant g = γ₀(sinc θ(u−½) + sinc θ(u+½))², γ₀ = 13/20:
   f ≤ g on I (margin 10⁻⁴), Λ = 2γ₀(1 + sinc θ)/θ = 1.78942 < 2; close pairs need x ≤ β = min_{[0,θ]} K² = 13707941/2²⁸.
3. Local data for frames of 8 points (gaps g_r ≥ θ): α ≥ 0; W_ij ≥ 0 (28) with Σ_i W_{i,i+s} ≤ 2 (s = 1…7);
   rewards κ_i ≥ 0, Σκ_i = x (paper: κ₃ only); storage functions B_r (r = 0…6), piecewise linear on the knots
   ¾ + k/8 (k = 0…13), zero at the end knots, Σ_r B_r ≡ 0; all integers in units of 10⁻⁶. For every frame:
   α(y₇ − y₀) + Σ W_ij K(y_j − y_i)² + Σ B_r(g_r) ≥ x. Then b = 7α. Proof: exact-integer branch and bound on a
   table v_k ≤ 2²⁸K² (cells 1/16384, t ≤ 16): 535,332,163 boxes.
4. Objective: maximise s. Paper: α = 627·10⁻⁶, x = 4361/500000, s = 1669159/2478195 = 0.6735382 (exact, rechecked).

## 2. What we ran (task 2)
| Variant | s | Status |
|---|---|---|
| Paper data | 0.6735382 | proved; our oracle finds its smallest slack 1.265·10⁻⁵ at gaps (1.033, 1.038, 1.966, 1.045, 1.966, 1.038, 1.033) |
| Same window, LP over all local data, margin 1.265·10⁻⁵ | ≤ 0.6735404 | numerical upper bound (relaxation after 17 cutting-plane rounds, still falling) |
| Same, zero margin | ≤ 0.6735490 | = previous + 0.679·margin; unprovable by box search |
| Window moved 0.005 along the LP gradient (a₁…a₁₄; C(f) = 1.3278761) | ≤ 0.6735482 after 11 rounds | original window: 0.6735438 at the same stage; still falling, gain ≤ 4·10⁻⁶ |
| Any 8-point certificate, paper's window | < 0.6735660 | proved in the paper (Prop. 7.3) |
| Paper's window, frames of 9 / 10 / 11 / 12 points | ≤ 0.6735989 / 0.6736298 / 0.6736472 / 0.6736665 | ours (ceil_range.py): ceilings from two periodic trial states, floating point; reproduces the paper's 0.6735660 (8 points) and 0.6738216 (all pairs) |
| Any pair-energy certificate, paper's window | < 0.6738216 | proved in the paper |
| Any window with V ≤ 10B, θ ≤ 0.998, this scheme | < 0.68185 | proved in the paper (Cor. 7.9) |

- LP (lp8.py): variables α, κ₀…κ₇, 28 weights, 84 knot values (|B| ≤ 0.0015; paper's |B| ≤ 0.00085), budget and
  Σ_r B_r = 0 rows, x ≤ β; objective by Dinkelbach; separation oracle = multistart pattern search over the 7 gaps.
  The LP value of every round is an upper bound for certificates of this form within these bounds (fewer constraints).
- Window gradient (grad.py, envelope theorem on the LP duals): |∇s| = 0.012, mostly a₂ (+0.0105). The stepped window
  did not deliver the predicted first-order gain (+6·10⁻⁵): at equal rounds it is ahead by 4·10⁻⁶ and the lead
  shrinks as cuts are added. Duals of an unconverged cut set overstate the gradient.
- Best admissible θ: for this window θ must lie in (0.5205, 0.9685) (necessary: K² ≥ x on [0, θ]; Λ < 2 at u = ½). All near-binding
  frames of the LP have gaps in [1.021, 2.936] (stepped window [0.984, 2.958]), and the paper's ceilings hold for every
  θ ≤ 1. So θ changes nothing.
- Verifier: the paper's C++ verifier built here (Apple clang 17, arm64, -O3) runs 138k boxes/s on these cores.
  We added checkpoint/resume (patch_cp.py; resumed 2M + 2M boxes = straight 4M, identical saved stack) and re-ran
  the paper's full 8-point search in 150M-box chunks. Chunk 1: 150,000,000 boxes (74,999,979 accepted), no failure,
  1458 s (103k boxes/s under load 11). The remaining ~385M boxes (~62 min) did not fit the time limit; the saved state
  is in `repro/` (resume: `verify8cp repro/second_window_params.txt kernel_lower.bin 150000000 cp_2.txt repro/checkpoint_150M.txt`,
  then chain; expected final counts 535,332,163 / 267,666,082).

## 3. What limits it
| Limit | Effect | Way around | Cost |
|---|---|---|---|
| Window ceiling (8 points) | s < 0.6735660 | none with 8 points | — |
| Certificate form (separable storage, knots 1/8) | LP ≤ 0.673549 even at zero margin | storage function of several gaps | new verifier tests |
| Verification margin | 10⁻⁵ of slack costs 0.68·10⁻⁵ in s | — | — |
| Frame length | 9 points ≤ 0.673599, 12 points ≤ 0.673667, all pairs < 0.6738216 | 8- to 11-dimensional box search | est. 10–30× more boxes per extra point |
| Window | known families < 0.674034 (Gebendorfer), all < 0.68185 | new f, kernel table, majorant, LP, search | ≥ 2–3 h per candidate here; our gradient step: no gain |
| Theory | c = 2 forced; other multiplicity credited 2 per unit without interaction | new idea (e.g. let off-line pairs very close to the line interact) | research |
| Energy input | BGST pair correlation only for test functions supported in [−1, 1], so C(f) ≥ 1.3275 (plain bound ≤ 0.6725) | wider support (Montgomery barrier) | open problem |

## 4. Novelty check (task 4)
arXiv math.NT listings Nov 2026 (140 entries to 5 Nov) and Oct 2026 (854): no newer simple-and-on-line record;
2610.08965 is still v1. Related: 2610.26841 (≥ 13/20 of zeros isolated, hence simple and on the line: weaker),
2610.24612 (exposition of 0.6725), 2610.23467 (simple or on the line > 8/9), 2610.22927 (consecutive gaps).

## 5. Files and runtimes
All scripts expect the e-print unpacked in ./src: `curl -s -L -o e.tgz https://arxiv.org/e-print/2610.08965; mkdir src;
tar xzf e.tgz -C src`. Python via `uv run --no-project --with numpy --with scipy`, single-threaded.
- `lp8.py` — cutting-plane LP (args: margin rounds seed free|k3 Bbound points theta; env WINDOW=a1,…). 11–17 rounds: 11–16 min.
- `grad.py` — window gradient from a saved cut set (`lp_cuts_*.npz` written by lp8.py); 3 s.
- `patch_cp.py` — writes verify_cp.cpp (checkpoint/resume) from the paper's anc/verifier/verify.cpp; build
  `clang++ -std=c++23 -O3 -DNDEBUG -DVERIFY_POINT_COUNT=8 -o verify8cp verify_cp.cpp`; usage
  `verify8cp params kernel_lower.bin [budget checkpoint_out [checkpoint_in]]`, exit 7 = checkpoint written.
- `repro.sh` — chained 150M-box chunks (each under 1800 s); `repro/` — parameter file and the checkpoint after 150M boxes.
- `ceil_range.py` — trial-state ceilings for frames of s+1 points (s = 7…11); 1 min; log in `logs/ceil_range.log`.
- `logs/` — LP runs (lp_run3.log: same window 17 rounds; lp_run4.log: same window 11 rounds; lp_step1.log: stepped window).
