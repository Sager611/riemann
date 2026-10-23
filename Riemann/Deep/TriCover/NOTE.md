# Non-separable families of homothetic triangles: covering factor 1.0633023 (explicit, 6003 copies) and λ ≥ 1.0633083 (limit)

Orca record hunt, 2026-10-22. Previous record: arXiv 2609.06558 (Y. Lian, F. Xue, Q. Yuan, "Non-Separable Homothetic Triangles, Part I", v1, 6 Sep 2026), Theorem 1.1: a 303-member family with covering factor 250000000/235141779 = 1.0631883498…

## 1. The problem
- A finite family of planar convex bodies is **non-separable (NS)** if no line that misses the union has members on both sides.
- For positive homothets F = (ξᵢ + τᵢK) the **covering factor** is λ(F) = min{s : ∪F ⊂ η + sK for some η} / Στᵢ. Goodman–Goodman (1945): λ ≤ 1 for NS families of disks. For a triangle Δ, λ_simplex := sup λ(F) over all finite NS families of positive homothets of Δ; by Bezdek–Lángi this is also the sup over all planar convex bodies.
- Known before: Λ₃ = μ = 2/3 + 2/(3√3) ≈ 1.0515 (Bezdek–Lángi 2016), λ_simplex ≤ 3/2 (Akopyan–Balitskiy–Grigorev 2018), and 2609.06558: families with 4, 5, 6, 9 and 303 members, the best being 250000000/235141779 ≈ 1.0631883; a continuous model with numerical value 1.0633083 that the authors state is **not** a proved lower bound ("its attainability and optimality remain open").
- The record is the best lower bound for λ_simplex, with any number of copies.

## 2. Results
**Theorem A (explicit family).** config.txt lists 2001 rational vectors; with their cyclic images they give an NS family of **6003** positive homothetic triangles with exact covering factor

  **λ = 500000000000000/470233160167623 = 1.063302298421…**

This is 1.139·10⁻⁴ above the previous record. verify.py (standard library, exact integers, 0.2 s, 25 MB) checks it.

**Theorem B (limit).** **λ_simplex ≥ 1.0633083237.** The paper's continuous candidate is a valid lower bound: chains with *equal* steps at the paper's displayed support parameters (c, d) are admissible for all large m and their covering factors converge to 1/S_cont = 1.06330832375…. This answers the paper's question "Do suitably optimized finite chains approach this continuous model and its candidate value?" (yes, even without optimizing the steps). Proof in §4; limit.py does the interval arithmetic.

| bridges per chain m | copies 3(m+1) | λ, exact (verify.py) |
|---|---|---|
| 100 | 303 | 1.063188861057 (paper, same structure, rounded to 10⁻⁹: 1.063188349868) |
| 500 | 1503 | 1.063284254826 |
| 1000 | 3003 | 1.063296278315 (equal steps at the paper's c, d: 1.063296176741) |
| 2000 | 6003 | **1.063302298421** (config.txt) |
| ∞ | – | ≥ 1.0633083237 (Theorem B) |

Optimized chains follow S(m) ≈ 0.9404609911 + 0.0107/m (S = total scale = 1/λ). Our optimizer reproduces the paper's table of optimized chains (§6: m = 5, 8, 12, 20, 50, 100) to all 10 printed digits.

## 3. Construction and check of Theorem A
**Coordinates (paper §2).** Σ = {z ∈ R³ : z₁+z₂+z₃ = 1, z ≥ 0}. For x ≥ 0 with |x|₁ < 1, T(x) = {z ∈ Σ : z ≥ x} is a positive homothet of Σ of scale 1 − |x|₁. If every coordinate has minimum 0 over the family, the smallest containing homothet is Σ itself (Lemma 2.1) and λ = 1/(sum of scales). T(z) meets θT(x) + (1−θ)T(y) ⊂ conv(T(x) ∪ T(y)) iff Σᵣ max(θxᵣ + (1−θ)yᵣ, zᵣ) ≤ 1 (Lemma 2.3).

**The chain (paper §5).** ρ(x₁,x₂,x₃) = (x₃,x₁,x₂). Supports Cᵢ = T(ρⁱ(c,d,0)); bridges P_{k,i} = T(ρⁱp⁽ᵏ⁾), k = 1..m, i ∈ Z/3. Lemma 5.1: the family is NS if
- {P_{k−1,i}, C_{i−1}} → P_{k,i} and {P_{k,i}, P_{m,i−1}} → P_{k−1,i} for 1 ≤ k ≤ m,
- {P_{m,i}, C_{i−1}} → P_{m,i+1},

where {A,B} → K means that K meets conv(A ∪ B). In words: if two supports lie on one side of a line, the red side grows along one chain (a fan around a corner of a support), jumps to the end of the next chain, and runs back down that chain (a fan around the end of the first chain) until it reaches the third support.

**Optimization.** The recurrence of §5 is parametrized by c, d, τ (position of the last vertex on a conic) and the steps h₁ ≤ … ≤ h_m, with two product equations. We minimize S = 3(1 − c − d + Σ h_k E_k): SLSQP at m = 100 from the paper's continuum candidate, then an augmented Lagrangian (L-BFGS-B) along m = 200, 300, 500, 1000, 2000, interpolating the step profile. gen.py does this.

**Exact certificate.** The two product equations are re-solved to 60 digits for (τ, common step factor), the recurrence is run in 60-digit arithmetic, and every nonzero coordinate is rounded down to the grid 10⁻¹⁵ and lowered by one more unit. Lowering coordinates enlarges the triangles, so every vertex contact of the recurrence becomes a strict relation (smallest slack 3.1·10⁻¹⁵). The total scale grows by less than 4·10⁻¹¹.

**verify.py** checks, with integers only: coordinates ≥ 0, scales > 0, all 6003 members distinct, every coordinate minimum 0; all 12003 relations of Lemma 5.1, evaluating M(θ) exactly at θ = 0, 1 and the breakpoints; and non-separability by forced colouring. It does not reuse the proof of Lemma 5.1: it runs the forcing from every colouring of the three supports. "All equal" must force every member to that colour; each 2-1 split must end in a contradiction. So every colouring compatible with the relations is constant, and Lemma 2.2 gives NS.

**Cross-check.** An independent exact test (union of the projections connected in every direction, at every critical direction and between consecutive ones, integer arithmetic) confirms NS for the 18- and 63-member certificates built the same way.

## 4. Proof of Theorem B
Fix c, d as the binary64 values of the paper's 0.3104341279739202 and 0.41586782506519804. For τ ∈ (0,1) let u, v, κ, D, R₀, Y₀, R_f, Y_f be as in §5 of the paper, and F(τ, s) = (log(R_f/R₀) − Ds, log(Y_f/Y₀) + κs) (the endpoint equations of the continuous model, §7).

**Certified by limit.py** (mpmath interval arithmetic, 60 digits):
- (a) a box B of radius 10⁻³⁰ contains exactly one zero p = (τ*, s*) of F (Krawczyk test, ‖I − C·DF(B)‖ ≤ 10⁻²⁶; so DF(p) is invertible). τ* = 0.67376304443337…, s* = 1.12130143018097….
- (b) on B all strict inequalities used by Theorem 5.3 hold: 0 < u < v < d < 1 − u, u + v > d, 0 < c < v, c + d < 1 − u (margin 4.18·10⁻⁸), d(κ+D) > κc; also κ < v.
- (c) g(t) = κR₀e^{Dt} + DY₀e^{−κt} − 1 < 0 on [0, s*], so E = Dκ/g > 0 is smooth there.
- (d) S_cont(p) = 3(1 − c − d + ∫₀^{s*} E dt) ≤ 0.9404609911001809 (trapezoid sum plus the error term s·h²/12·max|E''|, all enclosed).

**Lemma (equal steps).** Put G_m(τ, s) = (log(R_f/R₀) − m·log(1 + Ds/m), log(Y_f/Y₀) − m·log(1 − κs/m)). For all large m, G_m has a zero p_m = (τ_m, s_m) with p_m → p. The chain with c, d, τ_m and m equal steps h = s_m/m is an NS family of 3(m+1) positive homothets of Σ, with container Σ and total scale S_m → S_cont(p).

*Proof.*
1. ψ_m(x) = m·log(1 + x/m) satisfies |ψ_m(x) − x| ≤ x²/m and |ψ_m′(x) − 1| ≤ 2|x|/m for |x| ≤ m/2. Hence G_m → F in C¹ near p. Since DF(p) is invertible, x ↦ x − DF(p)⁻¹G_m(x) maps a small closed ball around p into itself and contracts, for large m. Its fixed point p_m is a zero of G_m, and p_m → p.
2. The inequalities in (b) are strict at p, so they hold at p_m for large m; also 0 < h < 1/κ. With equal steps the two product equations of Theorem 5.3 read (1 + Dh)^m = R_f/R₀ and (1 − κh)^m = Y_f/Y₀, i.e. G_m(p_m) = 0. Theorem 5.3 of the paper gives an NS family with all coordinate minima 0 and total scale S_m = 3(1 − c − d + Σ_{k=1}^m h·E_k), where E_k = Dκ/(κR₀(1+Dh)^k + DY₀(1−κh)^k − 1).
3. For t = kh ≤ s_m: |k·log(1+Dh) − Dt| ≤ D²ht and |k·log(1−κh) + κt| ≤ κ²ht. By (c) and continuity the denominators stay below some −δ < 0 for parameters near p. So E_k = E(kh) + o(1) uniformly in k, and Σ h·E_k is a Riemann sum that tends to ∫₀^{s*}E. Hence S_m → S_cont(p). ∎

So λ_simplex ≥ lim 1/S_m = 1/S_cont(p) ≥ 1/0.9404609911001809 > 1.0633083237.

Numerically, with 60-digit recurrences, the equal-step chains have S_m − S_cont ≈ 0.01074/m (m = 100, 1000, 10000), and A₀ = 1 − u − c − d ≈ 4.2·10⁻⁸ + 0.006/m, so the constraint that is nearly tight at p is not an obstacle. Rounding the m = 1000 equal-step chain as in §3 passes verify.py with λ = 1.063296176741.

Not claimed: that 1.0633083… is optimal, or that λ_simplex equals it.

## 5. Looking for a better pattern (all negative)
- **General contacts.** SLP over all coordinates and all segment weights with the paper's relations (m = 1…50), started at the optimum. No descent beyond 10⁻¹¹ (tolerance noise). The paper's vertex contacts are locally optimal.
- **Broken symmetry.** All 3(m+1) members free, same relations (m = 1, 2, 5): no descent.
- **Relation exchange.** Re-certify NS with any segment relations that hold (best-slack-first forcing), then re-optimize (m = 2): no descent.
- **The boundary c + d = 1 − u** (where the paper's continuum optimum sits) is forced by this pattern. In the case where a support is alone on its side, the lowest member known to be on the other side is the end of the previous chain. So the support's top vertex must reach that height. A connector or an extra pivot that relaxes this costs at least as much as enlarging the support.
- **Slack against the projection bound.** In every direction, the bridges' total projected length is at least 1.88 times the gaps between the supports' projections (ratio at most 0.53, at 108.75° in the planar picture; m = 100 chain). So the bridges are far from the trivial per-direction bound; a different contact pattern might gain, but we found none.
- **Random global search** with the exact projection constraint (n = 4, 6) is weaker than the known values (1.0521 < 1.0533 and 1.0501 < 1.0572), so it says nothing.

## 6. Files and reproduction
- config.txt: the 6003-member family (integers over 10¹⁵).
- verify.py: exact check of Theorem A (python3 verify.py; standard library only).
- gen.py: optimizer and certificate writer (numpy, scipy, mpmath). python3 gen.py 2000 rebuilds config.txt byte for byte (about 7 s CPU, 82 MB).
- limit.py: interval enclosures (a)–(d) for Theorem B (mpmath; about 4 s CPU, 24 MB).

## 7. Novelty check (partial)
arXiv API on 2026-10-22: "non-separable" ∧ homothetic, "covering factor" ∧ triangles, Goodman ∧ non-separable. Only 2609.06558 v1 and the Bezdek–Lángi papers (1602.01020, 2512.23651) come up; no Part II yet. A wider online check is still to do.

