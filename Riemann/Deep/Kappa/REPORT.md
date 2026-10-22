NO REFUTATION: OpenAI's fourth-moment lemma (7/8 paper, Lemma 18.1) does extend from 3/4 ≤ κ ≤ 1 to 13/18 ≤ κ ≤ 1. Xu–Sang and both Liu papers prove this correctly by re-running OpenAI's induction, and they only need κ ≥ 0.74977.

# Fourth moment below κ = 3/4: audit (2026-10-21)
Sources: OpenAI 7/8 TeX (build/paper.tex; §18 is l. 12477–14984), Xu–Sang 2610.15291v7, Liu 2610.12234 and 2610.26835 (arXiv e-prints), and our Lean port Riemann/Deep/QRH/OAI. Details are in [NOTE.md](NOTE.md); the exact margin checks are in kappa_margins.py (0.5 s, 16 MB).

## Where OpenAI's proof uses κ (TeX line numbers)
- **Prime-sum bound** (l. 12833–12898): the contour sits at (1+κ)/2, under β_* ≤ (1+κ)/2. The logarithmic-control lemma (l. 1531) holds for any line ≥ 1/2. Survives.
- **Width floor** (l. 12929–12936): z ≤ M/(6κ) ≤ 2M/9. The second inequality fails below 3/4 (it becomes 3M/13 at 13/18), but it is never cited again. Only κz ≤ M/6 is used. Survives.
- **Centered comparison** (l. 12947–12984): "6κ−1 ≥ 7/2, so z < M/21 < M/20". The 7/2 and the M/21 both fail below 3/4. Only z < M/20 is used, and that holds exactly when κ ≥ 13/18. It gives the printed bounds A_comp ≤ 23M/30 + ξ and A_comp + (6κ−1)z ≤ 14M/15 + ξ, with margin M/15. Survives, and is tight at 13/18.
- **Other steps that need only κ ≥ 1/6:** greedy slot removal for children (l. 14214–14280), A ≤ M (l. 14768), support errors 0 ≤ 6κ−1 ≤ 5 (l. 14827), and the parameter choices (l. 14885). Survive.
- **No κ at all:** Lemmas 18.2 (l. 13953) and 18.3 (l. 14545), both Poisson transforms, the Gauss-sum amplifier and the exceptional-row count.

## Lean port
- **The lemma is there.** It is the width induction in Energy/: terminal_certificate (CertifiedExistence.lean:42) has the hypothesis 3/4 ≤ κ. It is used only at κ = 3/4 + 2Δ > 3/4 (DetectorPlainMomentParameters.lean:24). Nothing below 3/4 is machine-checked anywhere.
- **Where 3/4 ≤ κ is consumed.** The hypothesis is passed through about 50 declarations, but only four arithmetic lemmas use it, all in the centered comparison: slot_budgets, balanced_reference_geometry, positive_slot_width_drop (z ≤ M/21) and robust_slot_budgets (16M/21 and 13M/14, reserve M/14). Everywhere else it only supplies κ ≥ 0, κ ≥ 1/6 or 6κ−1 ≥ 0.
- **Exact LP check at κ = 13/18.** With the port's own slack M₀/56, the robust gate gives 659/840 ≤ 5/6 and 799/840 ≤ 1. Porting the extension needs new constants (M/20, 23/30, 14/15, reserve M/15) but no new mathematics.
- **One use outside the lemma.** plain_capacity_bound (≤ 2/27, Hecke/DetectorRowCountCrossing.lean:106) belongs to OpenAI's row count. The refinements replace that count with their own (Xu–Sang ≤ 1/13; Liu ≤ c/12 < 1/5).

## The refinements: need and justification
- **Xu–Sang** need κ only at κ₁ = 29991/40000 = 0.749775. Their first stage licenses it, since β_* ≤ 69991/80000 = (1+κ₁)/2. Their Lemma "kappa-extension" and App. A (TeX l. 1867–1971) recheck exactly the steps listed above, then rerun the induction. Complete.
- **Liu 2610.12234** §28, Prop. "iv-plain" (l. 4333–4621): the same checklist, with explicit inequalities. Needs κ ≈ 0.74991. Complete.
- **Liu 2610.26835** Part II, Thm "M2-assembled" (l. 6971), with the comparison at l. 5814–5832: the same checklist. Needs κ♭ = 0.749773. Complete.
- **Room to spare.** At the κ actually needed, 6κ−1 ≥ 3.4986, so z < M/20.99, and OpenAI's own printed margins already suffice.
- **Side remark (routine).** The κ-steps alone allow any κ > 1/2, with comparison margin (6κ−3)M/(6(6κ−1)). This is useless, because the method's 13/15 floor forces κ ≥ 11/15 (Liu 2610.26835 l. 10211).

**First unproved or false step:** none found. Suggested NOGO line: "Kappa: L18.1 extends to 13/18 ≤ κ < 3/4; refinements' checklists complete".
