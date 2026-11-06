RECORD: lim sup_{x→∞} M(x)/√x > 2.0739001417 and lim inf_{x→∞} M(x)/√x < −2.0764178466, conditional on the correctness of Weisz's published balls for zeros 1–14400 (his own record has the same dependency; §4) (before 2.0617536515 / −2.0647155435, B. Weisz, Zenodo 23651459, 25 Oct 2026, unrefereed; refereed record 1.826054 / −1.837625, Hurst 2018).

How: Weisz's own two points y± and his zero data, with a better kernel. (i) The Jurkat–Peyerimhoff kernel is replaced by the best positive-definite kernel for the given phases (largest eigenvalue of a Toeplitz matrix). (ii) The kernel's support is widened from T = γ₁₄₄₀₀ = 13550.64 to T₁ = 18292.5. This brings in the 5,911 zeros 14401–20311, whose phases at y± are unknown (their ordinates are known only to ~25 digits, the phases need ~9,650), so each enters with its worst case −2|k(γ)||c|. Net gain: +0.0121 (lim sup) and +0.0117 (lim inf). A new point y (a stronger lattice search) is out of reach in 2 hours on one core (§3).

Node N10. Agent mlimsup (Orca), 2026-11-05/06, 23:38–00:40 CET. No Lean, no commits.

## 1. Statement and certificate
**Theorem used** (Odlyzko–te Riele 1985, p. 144, read from the scan; Ingham 1942). Let K ∈ C²(ℝ), K ≥ 0, K(−y) = K(y), K(y) = O((1+y²)⁻¹), and k(t) = ∫K(y)e^{−ity}dy with k(0) = 1 and k(t) = 0 for |t| ≥ T. If the zeros with |γ| < T are on the critical line and simple, then for every y₀:
lim inf M(x)/√x ≤ h_K(y₀) ≤ lim sup M(x)/√x, with h_K(y) = Σ_ρ k(γ)e^{iγy}/(ρζ′(ρ)) = 2 Σ_{0<γ<T} k(γ) Re(c_ρ e^{iγy}), c_ρ = 1/(ρζ′(ρ)).
- If RH fails or some zero is multiple, both limits are infinite (OtR p. 141). So the bounds are unconditional.
- The zeros with γ < 18292.5 are on the line and simple. Zeros 1–14400 are certified by Weisz's Arb data; zeros 14401–20311 are isolated here by Arb (acb.zeta_zeros, Turing's method), and N(18292.5) = 20311 (acb zeta_nzeros). The first 1.5·10⁹ zeros are known to be simple (van de Lune–te Riele–Winter 1986).

**Kernel** (`kernel_plus.txt`, `kernel_minus.txt`): a step function g on [−T₁/2, T₁/2] with M = 1350 cells of width h = 13.55, so T₁ = Mh = 18292.5. Its integer cell heights b_j are stored in the files. Then k = (g∗g̃)/‖g‖²: k is the linear interpolation, at the nodes t = dh, of R(d)/R(0) with R(d) = Σ_j b_j b_{j+d}, and k = 0 for |t| ≥ T₁.
- K = |ĝ|²/(2π‖g‖²) is ≥ 0, even (g is real), entire, and O(y⁻²) (g is a step function).
- k(0) = 1 and k vanishes for |t| ≥ T₁. So every hypothesis holds.
- Best–Trudgian (arXiv 1209.3843, §3.2) state the same admissibility rule: a nonnegative Fourier transform and support in [−T, T].

**Certified values** (`certify.py`, FLINT/Arb, 282 s, 0.9 GB; `certify.log`):

| Point | Zeros 1–14400, exact phases | Zeros 14401–20311, worst case | Certified bound |
|---|---|---|---|
| y₊ (Weisz) | +2.0766736799154 ± 4·10⁻⁴⁰ | −0.0027735381556 | h_K(y₊) ≥ 2.07390014175980135016775 ± 3.2·10⁻²⁴ |
| y₋ (Weisz) | −2.0789843657638 ± 5·10⁻⁴⁰ | +0.0025665191528 | h_K(y₋) ≤ −2.07641784661091652204704 ± 4.5·10⁻²⁴ |

- **Cross-check.** The same run re-evaluates Weisz's Jurkat–Peyerimhoff sums at y±. It reproduces his certificates (2.0617536514585282157830…, −2.0647155435291726514920…) to 10⁻⁴⁰.
- **Intermediate result** (optimal kernel only, T = 13550, no extra zeros): 2.0645146474 and −2.0674922729, both certified in the first run.

## 2. Method
1. **Phases.** Each phase γᵢy± mod 2π is computed in Arb at 32,819 bits, as in Weisz's verifier.
2. **Optimal kernel.** Since k(dh) = R(d)/R(0) and k is linear between nodes, Σᵢ aᵢk(γᵢ) = Σ_{d≥0} q(d)R(d)/R(0), where aᵢ = ±2|cᵢ|cos(γᵢy + arg cᵢ) and q(d) = Σᵢ aᵢ max(0, 1 − |γᵢ/h − d|).
   - For d > 0, R(d) = Σ_j b_j b_{j+d} counts each pair once. So the symmetric matrix carries q(d)/2 in each of its two off-diagonal positions: Q_jl = q(|j − l|)/2 for j ≠ l, and Q_jj = q(0).
   - Then Σᵢ aᵢk(γᵢ) = bᵀQb/bᵀb, a Rayleigh quotient of a symmetric Toeplitz matrix. The best kernel is the top eigenvector of Q (`kernel.py`).
   - By the Krein / Boas–Kac factorisation, every admissible kernel is some g∗g̃ with supp g ⊂ [−T/2, T/2]. So this is the optimum over all admissible kernels, up to the cell width.
   - The optimum has converged by M ≈ 1000. Complex (non-even) kernels add only 10⁻⁶.
3. **Wider support.** Zeros 14401–N(T₁) enter Q with aᵢ = −2|cᵢ|, using a rigorous upper bound for |cᵢ| (`ext_zeros.py`, Arb at 96 bits). k ≥ 0 at those zeros; the certificate uses |k| anyway.
   - Scan over T₁ (bound for y₊): 13550 → 2.06451, 14905 → 2.06910, 16260 → 2.07206, 17615 → 2.07356, 18292.5 → 2.07390.
   - The gain is still rising, but its steps shrink by a factor of about 0.7 per 340 units, so at most about +0.0004 remains. T₁ = 18292.5 keeps the rerun under 10 minutes.
4. **Kernel shape.** g is a cosine bump with a 4% pedestal at the edges. Compared with Weisz's kernel f(γ/13550.64), the final kernel (y₊) is higher at every aligned zero, and the gap grows with γ:

   | γ | 1000 | 2000 | 3000 | 5000 | 8000 |
   |---|---|---|---|---|---|
   | k | 0.986 | 0.947 | 0.887 | 0.721 | 0.431 |
   | f | 0.975 | 0.905 | 0.802 | 0.544 | 0.191 |
   - Why it helps: Weisz's aligned zeros reach γ ≈ 5000 (median 830, 90% below 1800). A wider support raises k at all of them, and the 5,911 worst-case zeros cost only 0.0026–0.0028.

## 3. Why no new point y
1. **What Weisz's points do.** In `.scratch` analysis, the 1,000 heaviest zeros (weight 2f(γ/T)|c|) are aligned, with rms phase errors 0.004 (top 10) to 0.08 (ranks 501–1000). Their weights sum to 2.0238 of the 2.0617.
   - Ranks 1001–14400 are nearly random and add +0.039 (y₊) and +0.042 (y₋).
   - |y₊| ≈ 1.3·10⁹⁶²⁵ and |y₋| ≈ 5·10⁹⁶²⁶.
2. **The ceiling is lattice reduction, not the size of y.**
   - All 14,400 weights sum to 2.5834.
   - A random-phase large-deviation model puts the best h over |y| ≤ 10⁹⁶²⁶ at about 2.5765.
   - So beating Weisz with a new y means aligning more than about 1,050 zeros. That needs reducing a lattice of dimension 1,050 or more with ~33,000-bit entries. Hurst's 800-dimensional fplll runs took 35 days each.
   - Weisz's γ data (9,860 digits) also caps |y| at about 10⁹⁸⁴⁰. This cannot be done in 2 h on one core in fpylll.
3. **Cheap moves at the given points.**
   - Local polishing y → y + δ: +10⁻⁶.
   - Complex kernels: +10⁻⁶.
   - Combinations such as 2y₋ − y₊ keep the alignment only with doubled errors, and they re-randomise the tail. Worse.

## 4. Trust chain
The record is conditional on Weisz's published balls for zeros 1–14400: γ to 9,860 digits, and c = 1/(ρζ′(ρ)) to 60 digits. His own record rests on the same balls. Everything else is certified here.

- **Certified here** (`certify.py`, Arb). From the input balls to h_K(y₊) ≥ 2.0739001417… and h_K(y₋) ≤ −2.0764178466…, every step is ball arithmetic.
  - This covers the phases γᵢy± mod 2π at 32,819 bits, the kernel (exact integers) and the sums.
  - Zeros 14401–20311 and upper bounds on their |c| are our own Arb computation, with N(18292.5) = 20311 by Turing's method.
  - The points y± need no trust: the theorem holds for every real y.
- **Imported from Weisz:** the γ and c balls of zeros 1–14400.
  - His documentation (data/README.md, THEOREM.md) states that each zero was isolated by Arb's `acb.zeta_zeros` at 16,384 bits, with the zero count checked by Turing's method.
  - It also states that each zero was refined to 32,768 bits by one certified interval Newton step, with N(B) ⊂ B checked for every zero, and that c was computed in ball arithmetic.
  - The package ships the data and a regenerator, not the logs of that run.
- **Rechecked here:**
  - SHA-256 of the npz and of both point files against his `SHA256SUMS.txt`. This checks integrity against the published package, not correctness.
  - 80 zeros at 256 bits: their index n, γ to 10⁻⁶⁰, and c to 10⁻⁵⁵.
  - **50 zeros at 33,010 bits** (`spotcheck2.log`, 1,603 s). These are the 40 heaviest by 2k(γ)|c| for y₊ and y₋ (the same 40 for both: zeros 1–72), plus 10 random zeros with γ from 345 to 13,192. The 40 carry 0.911 of the total kernel weight Σ2k|c| = 2.74, and about 0.91 of the 2.077 sum.
    - For each zero, one interval Newton step on a box B of twice the stored radius gives N(B) ⊂ B, so B holds exactly one zero, on the line by the symmetry of B.
    - Im N(B) lies inside the stored γ ball, and 1/(Bζ′(B)) lies inside the stored Re c and Im c balls.
    - Result: 50/50 on all three checks.
  - Earlier, zeros 7213 and 14400 passed N(B) ⊂ B (`spotcheck.log`).
  - Our pipeline reproduces Weisz's certified Jurkat–Peyerimhoff sums to 10⁻⁴⁰. This checks the arithmetic, not the data.
- **Not done:** a full-precision recheck of the other 14,350 zeros and their c balls. It costs 24–50 s per zero here (heights 14–13,550), about 6 core-days in total, and was out of scope.

## 5. Rerun (one command, about 6 min including the 89 MB download; one core; 0.9 GB)
`OMP_NUM_THREADS=1 uv run --no-project --with numpy --with python-flint==0.9.0 python3 Riemann/Deep/MLimsup/certify.py /tmp/mlimsup-data`
- It ends with "ALL CLAIMS VERIFIED". Then delete `/tmp/mlimsup-data`.
- **Rebuild the kernels** (the integer files come out identical; checked):
  1. `certify.py <dir> --phases p.npz`
  2. `ext_zeros.py 18300 e.npz` (about 4.5 min)
  3. `kernel.py p.npz e.npz` (needs `--with scipy`)
- **Data spot check:** `spotcheck.py <dir> --heavy 40 --random 10` (about 30 min at nice 19). `--dry` prints the selected zeros without computing.

## 6. Novelty pass (2026-11-06)
- `lst.py 'Mertens' math.NT --zenodo '"Mertens function"'`: no arXiv titles in 2026-10 or 2026-11. The Zenodo hits contain no oscillation bound beyond Weisz.
- Zenodo concept 23298194 has only v1.0.0 (11 Oct: 1.92036 / −1.92053) and v2.0.0 (25 Oct: 2.06175 / −2.06471). The GitHub repo wbeni95/mertens-function-bounds was last changed on 25 Oct (v2.0.0).
- HAL (query "Mertens", 20 newest): only Bi–Hou, hal-05838582, which concerns the counterexample size, not lim sup.
- arXiv API (`all:"Mertens function"`, `abs:Mertens AND abs:limsup`, `abs:"Mertens conjecture"`, `abs:Ingham AND abs:kernel AND abs:zeros`, `abs:"Jurkat-Peyerimhoff"`): the newest lim sup / lim inf bounds are Hurst's (1610.08551).
- **Novelty of the method.** Kernels used so far:
  - Fejér (Ingham);
  - Jurkat–Peyerimhoff (OtR 1985 §4.1, Best–Trudgian 2015, Hurst 2018, Weisz 2026);
  - Gaussian with an analytic tail bound, for the size of a counterexample (Pintz, Saouter–te Riele, Kim–Nguyen, Bi–Hou).

  Optimising the kernel for fixed phases, and charging the worst case for the zeros past the high-precision range, are routine in spirit (compare Lehman-type tail bounds). We found no use of either for these oscillation constants.
- **Caveat.** The new numbers rest on Weisz's points and data. His record is unrefereed and AI-produced. We spot-checked his data independently (§4).

## 7. Files
- `certify.py`: the rigorous certificate.
- `kernel.py`: the kernel optimiser.
- `ext_zeros.py`: zeros past 14400 with bounds on |c|.
- `spotcheck.py`: the data check.
- `kernel_plus.txt`, `kernel_minus.txt`: the kernels (1,350 integers each).
- `certify.log`, `spotcheck.log`, `spotcheck2.log`.
- Scratch analysis (Weisz's structure, random model, polishing): `.scratch/mlimsup/{analyze,kernopt,polish,cplx,extend_est}.py`.

RH relevance: N10 (explicit oscillations of the Mertens function; quantitative zero records).

DONE (2026-11-06 01:16 CET): review fixes applied (Toeplitz text, trust chain §4, spotcheck2: 50/50 at full precision).
