NO DISCOVERY: no two-piece or mollified detector beats B♭ = 0.874886647040, and no detector or reflected-energy change can break 13/15: that floor is also forced by the count-free generic-row end of the error sum.

# Detector restructuring vs the quasi-RH record (2026-10-20)
Sources: OpenAI 7/8 (Lemma 7.1 table (7.16), (14.5), Lemmas 14.2-14.3, 17.1-17.2, 18.1); Liu arXiv 2610.12234 (Lemmas 16.1, 22.1, Prop. 22.2); Liu arXiv 2610.26835 ((14.7)-(14.10), Thm 14.12).

## How a detector prime acts (OpenAI Lemma 7.1)
- Standard slot: restrict the theta index to p | c, stretch Z by q_p, subtract the rescaled probe. Row u gets q_p^(z-1)·χ_p(u); the principal term gets q_p^(z-1). The row range h grows by ℓ; the principal gains Z^(ℓ/6) (count of sixth powers). The subtraction is forced: it removes a leftover u-independent term.
- Every other local restriction, derived from the same table:
  - Gauss-modulus slot (p | s, Y → Y·q_p): row factor χ_p(u)/q_p, no principal gain, no row growth. Direct bound = standard probe with ly → ly+ℓ₂, but E is larger by ℓ₂(g+θ) ≥ 0 (pointwise version: ℓ₂(g+δ/2)). Numerically: ℓ₂ = 0.01 → B = 0.875545; ℓ₂ = 0.03 → 0.876862.
  - Cube slot (p | n): the leading term R = a_p²Q^(4-6x-6z) has no χ_p(u), so nothing is detected. The first detecting term (A = p^7) costs Q^(-2) per prime.
  - Marked-only slot (no subtraction): the leftover term q^(z-1-w) raises E by ℓ₂(δ/2-g) ≥ 0.
  - Stretch Z by q_p^γ (rescaled term at Z·q_p^(γ-1)): this is the same probe at scale Z^(1+(γ-1)ℓ), only a reparametrisation.

## (a) Two-piece detector
- Both pieces standard: this is Liu's system with unequal slots. The term with every slot marked (J = ∅) dominates, so the reflected excess depends only on the total ℓ. Long slots also break the slot-size limit of the plain fourth moment (OpenAI Lemma 18.1).
- Long piece in the Gauss modulus: dominated by simply lengthening Y, which Liu's free-length optimum already covers. Cube, marked-only and stretched long pieces: dead or equivalent (see above). Result: B ≥ B♭ in every case; the standard detector alone gives B = 0.874886646 (reproduced).

## (b) Mollified detector
- Weights λ_n on squarefree marked n: local multiplicativity plus Mellin separation turn the detector into a finite sum of standard products ("shapes").
- Cross terms: Cauchy-Schwarz over the finitely many shapes in the theta factor; each shape's tuples stay inside Lemma 14.3. The principal coefficient λ_nμ(n) ≥ 0, so there is no cancellation.
- A shape of length ℓ_s < ℓ costs E(ℓ_s) - (ℓ-ℓ_s)/6, which is at most -0.0095 at the optimum. So B equals B♭ and never goes below it.
- Non-product weights break the product-form hypotheses (Lemmas 14.1, 17.2, 18.1) and lose the prime boost of the row count. Mollifying the numerator through weights 1*λ on s leaves every exponent unchanged.

## (c) The 13/15 floor
- Direct side: an exact LP over the local exponents of Lemma 14.3 reproduces Liu's k(m,l) to 6·10⁻¹⁶. So the excess (1+3ℓ-2m)/4 is genuinely attained (extremal case: trivial dual, every mark active, H = (1-ℓ)/2); it is not a bookkeeping loss.
- New step (sketch, not in the papers): the active slots carry phases (q/p)₃² between marks ((14.5)). Split P = P1·P2; by reciprocity the P2 phase is a cubic symbol in P1·n. Cauchy-Schwarz over P1, then Heath-Brown's cubic large sieve, removes the L term of Lemma 14.2. The excess becomes (3+5ℓ-6m)/12 and the direct-side floor drops from 13/15 to 28/33 = 0.848485 (ℓ = 3/11, M = 8/11).
- The floor has a second side. Liu (14.7) at δ → 0 (a = 1/2, R = 1, K = 1, g = 0) gives B ≥ 2/3 + ℓ + b/6, which is Liu's floor endpoint C0 = -1/4 + b/6 + 5ℓ/4 < 0. This bound uses no count and no moment, so it holds even under GRH for every row L-function. Together with B ≥ 11/12 - ℓ/4 it forces B ≥ 13/15 for every length choice. The optimum over free lengths is 0.866667 with either direct-side excess, and also under a density-hypothesis count.
- So only cancellation across rows can go below 13/15; the method sums rows in absolute value. B♭ is unaffected either way: at its optimum ℓ = 0.1671 < 1/5, no floor term is active.

## Is the argument complete?
- Complete at the exponent level: the reductions of the Gauss-modulus, cube, marked-only, stretched and mollified detectors to points of Liu's refined system, whose infimum is B♭ (Thm 14.12).
- A sketch: the trilinear step. It still needs the exact tuple phases after sector fixing and the D, r reductions in the proof of Lemma 14.2. It moves only the direct-side floor.
- The inverse-moment cap r+2z ≤ m comes from the u = za column case of (17.13): the dual grows by 2za there. It is genuine, and the trilinear step does not touch it.

Scripts: reflected_lp.py (7 s, 30 MB) and detector_variants.py (170 s, 53 MB); outputs in .scratch/qrhimp3/.
