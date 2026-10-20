NO DISCOVERY: B♭ = 0.874886647040 (Liu, arXiv 2610.26835) stands. None of de Faveri–Pandey (arXiv 2610.01035), Guth–Maynard, Chen–Gupta–Li or Heath-Brown's twelfth moment fills a slot of the refined system, and none would lower B even if transferred at its exact published exponent.

# New moment and large-value inputs vs the quasi-RH record (2026-10-20)

## Why nothing plugs in
- Every moment, large-value and count input of the refined system is an average over one family: the sextic residue characters ψ_u (u sixth-power-free, N(u) ≍ U) over F, varying in conductor. Heights are at most Z^τ with τ → 0 (OpenAI 7/8, Lemma 8.1), and the target L(s,η) stays fixed while Z → ∞ (Prop. 2.1; constants may depend on η). Height (t-aspect) bounds for ζ are never used.
- de Faveri–Pandey, Guth–Maynard and Heath-Brown bound ζ in height only; no transfer to this family is known. Chen–Gupta–Li (arXiv 2507.08296) carry Guth–Maynard to all characters mod one q (density exponent 7/3). The ψ_u have different moduli; a family containing them has ≍ U² members, and there even the density hypothesis gives a count U^(2(1−δ)) > U at the critical δ.
- de Faveri–Pandey's main tool (Cor. 1.4, large values of sums with phases √x + ηx^(3/2), from Atkinson's formula) has no counterpart: the pipeline's direct bound is a mean-square Gram estimate.

## Even a perfect transfer at the exact exponents changes nothing (`new_inputs_check.py`)
Binding class: δ = 0.46934, x = 1/2, κ = 0.749773, count R = 0.59033, numerator K = (3R+1)/4 = 0.69275. There E does not depend on b (Liu Lemma 14.10), so B moves only if R or K drops at this point.
- Numerator, by level sets: Σ_B |L| ≪ max_v U^(v + min(R, n(v))), where U^n(v) bounds #{u : |L| ≥ U^v}. The fourth moment gives n = 1 − 4v; the max sits at v* = (1−R)/4 = 0.1024.
  - de Faveri–Pandey's large values (Thm 1.3: T^(24/17) V^(−128/17) for V in [T^(7/60), T^(5/38)], which contains V = T^(1/8); Thm 9.1) start at v = 7/60 = 0.1167 > v*, so K is exactly unchanged.
  - Their moment bounds follow from these and also lose. Hölder with M_k ≪ T^(1+η_k) helps only if η_k < (k−4)(1−R)/4:
    k = 6 needs 0.2048, has 1/4 − 1/60 = 0.2333; k = 8 needs 0.4097, has 0.4885; k = 12 needs 0.8193, has 1 − 520/445121 (Heath-Brown: 1).
- Row count: a class has at most U^(A(a)(1−a)) rows, a = (1+δ)/2 = 0.7347, if A is a density exponent for the family. This beats R only if A(a) < 2R/(1−δ) = 2.2249.
  - Best known for ζ at a: Guth–Maynard 15/(3+5σ) = 2.2477 (Tao–Trudgian–Yang Table 2); also 30/13 = 2.3077, Ingham 2.3709, Chen–Gupta–Li 7/3. Count 0.5964 > 0.5903.
- Reruns over the full rectangle (b free, ℓ bisected, κ = 2B−1): the baseline, each input alone, and all of them together give B = 0.874886647 (identical to 10⁻¹²).
- Independent check: the tangent system of (14.8) alone gives B = 0.874886647040, the quartic root to 4·10⁻¹⁴; max E on the rectangle at (b*, ℓ*) = (0.086535, 0.167120) is −8·10⁻⁹.

## Where a Guth–Maynard-type estimate would matter (hypothetical)
Guth–Maynard's large-value theorem itself (N²V⁻² + N^(18/5)V⁻⁴ + TN^(12/5)V⁻⁴) would enter the inverse and plain counts e_I, e_S (Liu arXiv 2610.12234 (34.6)). Suppose it held verbatim for the sextic family and were applied to products of the two detecting polynomials and the selected prime sums. Then the critical count would drop to 0.5647 and B to about 0.8744 (`gm_transfer_count.py`, coarse grid, which gives 0.87489 without it). No such theorem exists. It would also presuppose a sharp mean-value theorem for the family, which is open too (REPORT.md item 4).

## Sources
de Faveri–Pandey arXiv 2610.01035 (Thms 1.1–1.3, 9.1, 10.1, Cor. 1.4); Guth–Maynard arXiv 2405.20552 (Thms 1.1, 1.2); Chen–Gupta–Li arXiv 2507.08296 (Thms 1.1, 1.2); Tao–Trudgian–Yang arXiv 2501.16779 (Table 2); Heath-Brown, Quart. J. Math. 29 (1978) 443–462; Liu arXiv 2610.26835 (§14) and 2610.12234 (§34); OpenAI 7/8 (Prop. 2.1, Lemma 8.1).
arXiv API searches on 2026-10-20 ("quasi-Riemann", "zero-free half-plane", Guth AND Maynard) found nothing below B♭.
Scripts: `new_inputs_check.py` (100 s, 55 MB) and `gm_transfer_count.py` (70 s, 30 MB); outputs in .scratch/qrhimp2/.
