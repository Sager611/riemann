# DERIVATION: complex amplifiers, the midpoint twist, and the corrected ceiling of Bui–Milinovich's criterion

Sources: BM = Bui–Milinovich, arXiv 1410.3635v3 TeX (`.scratch/gapbar/bm/main.tex`, "l." = line); GapCeil =
`Riemann/Deep/GapCeil/` (DERIVATION.md §1-7: X, X_i, τ, c = κπ, φ = I^{2r}f, Γ arrays, certificates). Notation as there.

## 0. Mollifier shifts with real P are admissible as printed
- BM's Lemma main (l.370-384) gives I_{α,β}(M) = ∫ζ(½+α1+it)ζ(½+α2+it)ζ(½+β1−it)ζ(½+β2−it)M(½+α3+it)M(½+β3−it)w,
  uniformly for all shifts ≪ 1/L. In c(α, β) the mollifier shifts enter only through
  y^{−(α3+β3)x−α3(x1+x2)−β3(x3+x4)}.
- With α3 = −β3 = iκπ/(2L) (the mollifier at the midpoint, M(s + iκπ/(2L))) this factor is
  e^{−i(κπ/2)ϑ(x1+x2−x3−x4)} = e^{ic(r2−r1)/2} (X = ϑx, r1 = X1+X2, r2 = X3+X4, c = κπ). This is the twist of §3.
- F(t) = e^{ivtL}ζ(s)ζ(s+iκπ/L)M(s+iκπ/(2L)) vanishes at t_n and t_n − κπ/L, which is all the reduction (l.164-202) uses.
  F' is handled by Q(L⁻¹(∂_{α1}+∂_{α2}+∂_{α3})) (l.402). So h1 > 1 for this F with real polynomial P gives Λ > κ.

## 1. Complex amplifiers are admissible
- BBLR's twisted fourth moment, as quoted by BM (Thm, l.500-517), is stated for Σ_{h,k≤y} a_h conj(a_k)/√(hk) with
  arbitrary complex a_h (a_h ≪ h^ε, y = T^ϑ, ϑ < 1/4).
- BM's evaluation is linear in the coefficient products: I_1 = Σ d_r(h)d_r(k)P[h]P[k]/(h^{1/2+α3}k^{1/2+β3})… (l.531),
  then Σ_{i,j} b_ib_j i!j!/(log y)^{i+j} J_{i,j} (l.537, l.557), P(x) = Σ b_jx^j (l.250). With complex b_j the same steps
  give b_i conj(b_j): the main term is the sesquilinear form with BM's kernel, P(w1)·conj P(w2) (or the conjugate on the
  other side: "orientation"; both give the same sup, and P ↔ conj P swaps them).
- Wirtinger's inequality holds for complex-valued functions (BM l.140-155). The end-piece step (l.194) needs only
  |M| ≪ T^{ϑ+ε}.
- Hence every complex polynomial P is admissible and h1 = c0/(κ²c1) > 1 still gives Λ > κ.
- No real-P assumption appears in BM l.243-255 ("P(x) = Σ b_j x^j is a certain polynomial").

## 2. The kernel is Hermitian, not real
- The swap (X1,X2) ↔ (X3,X4) maps ρ_j to conj ρ_j, so the Taylor arrays satisfy Γ = Γ^†
  (float check: |Γ − Γ^†| ≤ 1.5·10⁻¹⁰ on i, j < 6; `herm.py`).
- GapCeil kept Re Γ (`gcore.py` `sym`, `certr.py`, `gcert*.py`): that is the real-φ problem.
- Im Γ is large: r = 1, κ = 3.3, raw Taylor coefficient a[0,1] = −2.156 − 10.061i.
- For complex φ the forms are Hermitian, Q_L(φ) = Σ Γ_ij conj(I^iφ)^T W_{b+p} (I^jφ). Their top is larger:
  r = 1, κ = 3.2721: complex 1.0226561, real 0.9999978 (`ctest.py`, N = 12).

## 3. The midpoint twist makes the kernel real (exact identity)
With r1 = X1 + X2, r2 = X3 + X4 (U2 − U1 = r2 − r1), σ = τ1 − τ2:

  e^{ic(r2−r1)/2}·U_L conj(U_R) = 2cos(c(X1−X2−X3+X4)/2) − 2cos(c(X1−X2+X3−X4)/2 + cσ).

Proof: the factor sends A = e^{ic(X1−X3)} to e^{ic(X1−X2−X3+X4)/2}, the conjugate of the image of B = e^{ic(X2−X4)},
and C = −e^{ic(X2−X3)}e^{−icσ} to the conjugate of the image of D = −e^{ic(X1−X4)}e^{icσ}. The τ-box and S do not change.
With a = c(X1−X2)/2, b = c(X3−X4)/2, θ = cσ/2 the right side is 2cos(a−b) − 2cos(a+b+2θ), which is the **product form**

  e^{ic(r2−r1)/2}·U_L conj(U_R) = 4 sin(c(X1−X2+σ)/2)·sin(c(X3−X4+σ)/2)

(a left factor times the same function of the right variables). At X_i = 0 it is the multiplication part 4sin²(cσ/2).
By §0 this is BM's own kernel for the midpoint mollifier with real P.
1. **Reduction.** For f = e^{icU/2}g, conj f(U1)·f(U2)·ρ_j = conj g(U1)·g(U2)·ρ̃_j with ρ̃_j real symmetric
   (orientation L; e^{−icU/2} for the other). A Hermitian form with a real symmetric kernel takes the same top on complex
   and on real vectors (Q̃(g) = Q̃(Re g) + Q̃(Im g)). So: **complex BM pencil = real pencil of ρ̃.**
2. **Complex sup = midpoint real sup** (L² closure). U = log h/L, so e^{icU/2} = h^{iκπ/(2L)}: real g ↔ the midpoint mollifier
   M(s + iκπ/(2L)) of §0. After t → t − κπ/(2L), F(t) = e^{ivtL}ζ(s − iκπ/(2L))ζ(s + iκπ/(2L))M(s).
   - The twist does not keep polynomials (reviewer, §9). At a fixed degree the two classes differ slightly
     (r = 2, degree 8, κ = 3.25: midpoint real P 1.01560266, complex P 1.01558240).
   - Both are admissible (§0, §1), and they have the same closure sup.
3. **Smaller Taylor data.** r = 2, κ = 3.335: max |Γ̃_ij| (i, j < 6) = 1.31 against 643 for Γ. The twist removes the
   e^{±icr/2} oscillation; twisted a[0,1] = −0.0602 (real) against −0.0602 − 0.2831i.
4. **Checks.** Float: r = 2, κ = 3.335, 14 modes: real pencil of ρ̃ 0.9932086 = complex pencil of ρ 0.9932086
   (`twcheck.py`). Arb: the imaginary parts of the twisted arrays contain 0 in every run (IMCHECK).

## 4. Consistency check: an exact symmetry of the true mean values
- For real P let G_β(t) = ζ(s)ζ(s+iα)M(s+iβ), α = κπ/L. Then conj G_β(−t−α) = G_{α−β}(t), and |F'| maps to |F'| at the
  same v. Mean values over negative t equal those over positive t (conjugation). So c_j(v; β) = c_j(v; α − β).
- In f-language M(s+iβ) is f·e^{iμU} (μ ∝ β), and μ ↦ c − μ is the symmetry (g ↦ conj g in the twisted gauge).
- `mutest.py`, r = 2, κ = 3.31, real g of degree 10:

| μ/c | 0 | 1 | 0.2 | 0.8 | 0.5 | −0.5 |
|---|---|---|---|---|---|---|
| kernel K (Hermitian) | 0.9610710 | 0.9610661 | 0.9698746 | 0.9698755 | **0.9778261** | 0.9463031 |
| Re K (GapCeil) | 0.9610710 | 0.9380853 | 0.9589095 | 0.9425451 | 0.9508246 | 0.9508246 |

  K respects the identity (to the float accuracy, 5·10⁻⁶) and peaks at the midpoint μ = c/2.
  Re K, used for complex amplifiers, is symmetric about μ = 0 instead and violates the identity (0.9611 against 0.9381).
  So Re K is not the complex extension of the true main term. For real P the two kernels agree.

## 5. Ceilings and degrees (float)
- **Closure ceilings** (complex = twisted; N = 14 φ-modes, D = 36; `cscan.py`, `logs/cscan_r123_N14.txt`):
  r = 1: 3.297505 (real 3.272097); **r = 2: 3.326575** (real 3.302788); r = 3: 3.238189 (real 3.226093).
- **κ_max by degree** (complex P, ϑ = 0.2499; `copt.py`, `logs/cdeg_scan.txt`):

| deg P | 4 | 8 | 12 | 16 | 20 | 24 | 32 | 48 | 64 |
|---|---|---|---|---|---|---|---|---|---|
| r = 2 | 3.216396 | 3.268823 | 3.291084 | 3.302541 | 3.309199 | 3.313404 | 3.318217 | 3.322267 | 3.323885 |
| r = 1 | 3.254406 | 3.280327 | 3.288268 | 3.291705 | | | | | |

  Complex degree 16 already beats the real degree-60 record (3.300010).

## 6. Record certificates (`certmid.py`: midpoint, real P; `certc.py`: complex P)
- `certmid.py` is GapCeil's certr.py with the twisted Γ arrays of gtw.py; real P, no orientation question.
  Results (r = 2, ϑ = 2499/10000; IMCHECK true):
  - degree 80, v = 1.367258, κ = 831/250: **h1 = 1.00056020750 ± 4.5·10⁻¹⁵**, so **Λ > 3.324** (59 s, 216 MB,
    `logs/certmid_r2_d80_k3.324.txt`; float 1.00056049);
  - degree 64, v = 1.366889, κ = 3323/1000: h1 = 1.00071833023 ± 1.1·10⁻¹² (59 s, 196 MB).
  Midpoint real-P κ_max by degree (r = 2): 3.318217 (32), 3.322267 (48), 3.323885 (64), 3.324691 (80).
- certr.py with the Hermitian form: Q_L = Σ ReH_ij S_ij − ImH_ij T_ij, S_ij = y_i·W·y_j + z_i·W·z_j,
  T_ij = y_i·W·z_j − z_i·W·y_j (y = I^i Re φ, z = I^i Im φ), (ReH, ImH) the Hermitian part of the Arb Γ arrays.
  The tail bound ε|φ|² uses |φ|² = Σ(|Re φ_k|² + |Im φ_k|²)ϑ/(2k+1) (the majorant bounds |Γ| itself).
- **Orientation.** If the true form conjugates the other side it equals Q_L at conj φ, and conj P (also admissible)
  then gives the certified value. Λ > κ holds either way.

| r | deg P | ϑ | κ | v | h1 (Arb) | time, peak |
|---|---|---|---|---|---|---|
| 2 | 32 | 2499/10000 | 663/200 = 3.315 | 1.364289 | 1.00261659841 ± 3.7·10⁻¹² | 24 s, 172 MB |
| 2 | 64 | 2499/10000 | 1661/500 = 3.322 | 1.366893 | 1.00153083832 ± 4.8·10⁻¹² | 31 s, 168 MB |

- Exact coefficients (A_n + iB_n of P in the shifted Legendre basis) are printed in `logs/certc_r2_d*_k*.txt`
  and stored in `logs/copt_r2_d*_k*.json`.
- Checks (degree 32, κ = 3.315): Re P alone gives h1 = 0.959514245386, conj P gives 0.927680374808.
  The imaginary part enters at first order.

## 7. Barrier certificates for every complex amplifier (`gtw.py`, `runcerttw.py`)
- GapCeil's real certificate (gcert_r.py, GapCeil DERIVATION §7: weighted Jacobi basis, high-mode bound, Schur step,
  v-covering) runs unchanged on ρ̃. Only the four phase pairs of the Γ arrays change, to
  (c/2, −c/2 | −c/2, c/2), (−c/2, c/2 | c/2, −c/2), (c/2, c/2 | −c/2, −c/2), (−c/2, −c/2 | c/2, c/2).
- All |phases| ≤ c, so the Taylor-tail majorant |Γ̃^(r)_ij| ≤ 4(i+j+8)²((i+2r)(j+2r))^{2r}γ^{i+j}, γ ≥ 2c+2, stands.
- κ²M̃1(v) − M̃0 ⪰ 0 on real L² gives it on complex L², hence for every complex f = e^{icU/2}g: BM's criterion fails at
  that κ for every complex (so every real) amplifier, every v and every ϑ ≤ 1/4.
- A direct Hermitian certificate (`gcertc.py`, real 2N embedding) failed at r = 2, κ = 3.335
  (`logs/certC_r2_3.335.txt`): the Im Γ_01 terms are linear in t_1 (integration by parts handles only the real part),
  and C rose to 1.5-3.0, above δ_D, in 10 of 16 cell operators. The twist removes them: C ≤ 0.98
  (`logs/certTW_r2_3.335.txt`).

- **Results** (REPORT.md (iii) has the full table): certified for every complex amplifier at
  - r = 1: κ ∈ {3.30, 3.335, 3.35, 3.5, 4, 5};
  - r = 2: κ ∈ {3.335, 3.34, 3.35, 3.36, 3.38, 3.42, 3.45, 3.5, 4, 5};
  - r = 3: κ ∈ {3.335, 3.35, 3.5, 4, 5}.
- N = 72, D = 50-56 (r = 1, and r = 2 for κ ≤ 3.5); N = 48 or 40, D = 76 (r = 2, κ = 4, 5); N = 48 or 40,
  D = 66-84 (r = 3). Max C/δ_D was 0.15 (r = 1), 0.40 (r = 2) and 0.74 (r = 3). Peaks were 114-245 MiB.
- r = 4 and 5 at κ = 3.335 failed (N = 40, C/δ_D up to 1.28 and 1.77). The r = 5 run peaked at 316 MB, above the cap.

## 8. Float landscape (`cland.py`, `ceil45.py`)
- Closure ceilings, complex = midpoint (real, GapCeil):
  - r = 1: 3.297505 (3.272097);
  - r = 2: 3.326575 (3.302788);
  - r = 3: 3.238189 (3.226093);
  - r = 4: 3.147354 (3.141892);
  - r = 5: 3.073662 (3.071102).
- sup_v λ/κ², r = 2, complex (real): 1.02194 (1.00229) at 3.30, 0.99321 (0.97419) at 3.335, 0.98131 (0.96257) at 3.35,
  0.87475 (0.85857) at 3.5, 0.63834 (0.62848) at 4, 0.40199 (0.39567) at 5, 0.13906 (0.13769) at 8.
  Monotone in κ; the complex gain shrinks with κ.

## 9. Review
- **Reviewer 1** (run f06df76c, inline argument):
  - (i) CORRECT;
  - (ii) GAP in claim (c) (twist ≠ polynomial), repaired by §0;
  - (iii) CORRECT given the certificates;
  - overall SIGNIFICANT.
- **Reviewer 2** (run bcf576f8, §0 with BM l.360-420): CORRECT, ROUTINE (a specialisation of the printed lemma; the
  numerical gain is a record). Under Σ∂_α and Σ∂_β the shift differences cancel, so c_1 keeps (v − S)².

## 10. What is not covered
- κ between grid points, κ beyond the grid, and large r (no uniform-in-r bound), now for the twisted pencil.
- Ball-κ was tried (`runcerttwball.py`, κ ∈ [3.4, 3.45], r = 2) and fails: entry radii reach 10⁵⁷-10⁶⁴ through the
  cancelling Taylor convolutions. The twist does not cure this; the high-mode bound alone survives (C ≤ 1 ± 0.5).
  Next step: a Taylor model in κ (exact ∂_κ arrays at κ0 plus a majorant remainder). With M1 ⪰ 0,
  κ²M1(κ) − M0(κ) ⪰ κ0²M1(κ) − M0(κ) for κ ≥ κ0, so only the κ-dependence of the operators must be controlled.
- Inputs as printed: BM's Lemma main and Theorem smoothed (with complex P as in §1), BBLR's twisted fourth moment.
