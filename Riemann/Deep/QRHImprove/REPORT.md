NO DISCOVERY: the best published exponent, B♭ = 0.874886647040 (Liu, arXiv 2610.26835, 18 Oct 2026), is not beaten. Only parameter-level and hypothetical-input checks were done; no new rigorous argument was found.

# Quasi-RH exponent: how the proof works, where it is stuck, what was tried (2026-10-20)

## How the proof works
- **The sum.** A smoothed sum over F = Q(√−3) of cubic-theta coefficients (cubic Gauss sums) at the indices c·n³. It is twisted by the target character η and by sextic residue symbols, and multiplied by K selected primes of total length ℓ. Averaging lengths are X = Z^lx and Y = Z^ly; b = ly − lx.
- **Two estimates.** Theta reflection, an additive Gram bound and Cauchy–Schwarz give |J| ≪ Z^(lx/2+b/12). Poisson summation turns the sum into rows u, each carrying L(w,χ(u))/L(s,ηχ(u)). Row u = 1 is a Mellin integral of 1/L(s,η) with exponent C(s) = s − (4+b)/6.
- **Contradiction.** A power saving that does not depend on η continues 1/L(s,η) to the left of β* (the largest real part of a zero in the family), which is impossible. Matching the upper bound to C(B) gives B = 11/12 − ℓ/4 when lx + ly = 1 − ℓ.
- **Bad rows.** Rows of norm U are grouped by zero class δ = 2a − 1 and prime-sum size x. A zero forces large inverse and plain character polynomials. Three moment bounds then give #rows ≪ U^Rκ(δ,x): the inverse moment, sixth-power amplification, and the plain fourth moment (which needs β* ≤ (1+κ)/2). The numerator is averaged by Hölder against the amplified fourth moment (Xu–Sang).
- **Bottleneck inequality** (Liu (14.8), at row scale h = (1+b+3ℓ)/2). For every δ in [1/8, 5/6] and x in [0, 1/2]:
  E = −1/4 + b/6 + 5ℓ/4 + δ/2 + δx(ℓ + (h/4 − ly/2)/(3κ)) − (3h/4)(1 − Rκ(δ,x)) ≤ 0.
  It binds at δ = 0.4691, x = 1/2, where E does not depend on b. Setting κ = 2B − 1 closes the loop and gives B♭. There Rκ = 0.591, against 1 − δ = 0.531 under the density hypothesis, so the bad-row count holds the exponent back.

## What was tried
Scripts: `refined_system.py` (~35 s, 50 MB) and `.scratch/qrhimp/sens3.py`.
1. **Reproduction.** Rebuilt the system: B = 0.874886647040 = B♭ to 12 digits, same binding point.
2. **Freedoms Liu's optimality proof leaves out** (his Remark 14.13). None helps:
   - Unequal lengths (lx + ly ≠ 1 − ℓ): at fixed B every shift pushes max E above 0 (about +3·10⁻⁴ per 5·10⁻⁴ in ℓ). B♭ sits at a kink of the low-side bound.
   - Pointwise or averaged numerator chosen class by class: no change.
   - The shift θ = q/(6κ) is already optimal. Extra detector cutoffs give independent witnesses, so they gain nothing.
3. **Floor for this sum.** For any lx, ly, ℓ, b the upper-bound side alone forces B ≥ 13/15 ≈ 0.8667. The cause is the reflected-energy excess (5ℓ−1)₊/4 once ℓ > 1/5.
4. **Value of hypothetical new estimates** (same pipeline and fixed point; published system gives 0.874887):
   - row count at the density hypothesis, 1 − δ: 0.866667 (hits the 13/15 floor)
   - numerator averaging like Lindelöf over the rows: 0.866667 (floor)
   - numerator via a sixth moment: 0.869351
   - inverse-moment prime cost r+2z → r+z: 0.872756
   - amplification exponent 5/6 → 2/3: 0.873269; 5/6 → 0.8: 0.874627
   - plain-moment prime cost 6κ → 3κ: 0.874379

   Each needs a new mean-value theorem for the sextic family; none exists in the literature.
5. **Outside inputs.**
   - Guth–Maynard bounds large values in height (n^(it)). The count here is a zero-density bound in conductor for a thin family (U sextic characters of conductor ≍ U, height ≤ Z^τ), and no transfer is known.
   - Halász–Montgomery needs rδ ≥ 1 (r ≥ 2.1 here), and full-family Hecke zero-density gives an exponent ≥ 2(1−δ) > 1. Both are useless.
   - The Pa^(1/6) Gram loss is a genuine main term from sixth-power frequencies (Xu's matching lower bound), so it is not slack.
6. **Other fields and kernels.** The method needs explicit metaplectic theta coefficients, known only for cubic covers; Liu already covers every number field through K(ζ12). The argument has no positivity step, so there is no kernel to swap.

## Sources
OpenAI 7/8 (github.com/openai/math/tree/main/preprints, 30 Sep 2026); Liu arXiv 2610.12234 (0.874957019); Xu–Sang arXiv 2610.15291 (0.8748867); Liu arXiv 2610.26835 (B♭, §14, Theorem 14.12, which also cites Newton's 0.87488664962). A web search on 2026-10-20 found nothing below B♭.
