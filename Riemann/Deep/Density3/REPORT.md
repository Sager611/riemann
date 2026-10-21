NO DISCOVERY: the ring class density exponent e(σ) = 2(1−σ)/(3−2σ) cannot be improved with second moments (proved optimal for that toolkit); known subconvexity improves it only for σ > 0.9995 (refereed input) or σ > 0.9915 (Yang's claimed ¼ − 1/120), which is routine and leaves the almost-all threshold 3/2 unchanged.

# Density3 report (2026-10-21)
Proofs: [NOTE.md](NOTE.md). Brief: `.scratch/orca/density3-brief.md`. Paper proofs, not refereed, no Lean. Notation: D < 0, N_D(σ;t) = number of zeros with β ≥ σ, |γ − t| ≤ 1, summed over the ring class characters of D; e(σ) is the current per-D exponent (RCApps2), already below Combos' 2σ(1−σ).

## What was proved
- **Prop 1 (barrier).** With the per-D large sieve on k-th powers, the first moment, the twisted second moment and any mollifier length, the best exponent is exactly e(σ). Combos' own Carlson argument, re-optimised (cutoff |D|^{1/(2(3−2σ))}, best k per block), already gives e(σ): its worst block has length |D|^{1/(2(3−2σ))}, where the large sieve and its square give the same bound. Checked numerically (`qa1.py`).
- **Thm 2 (fourth moments).** If Σ_{χ≠1}|L(½+iu,χ)|⁴ ≪ |D|^{a+ε}, then N_D(σ;t) ≪ |D|^{min(e(σ), (1+a)(1−σ)/(2−σ))+ε}; this helps iff a < 1/(3−2σ). Only a = 1 is known (large sieve; also the full GL₂ family by positivity), and a = 1 never helps.
- **Thm 3 (subconvexity, unconditional instances).** If L(½+iu,χ) ≪ (|D|(|u|+3)²)^{¼−δ+ε} for all χ ≠ 1, then N_D(σ;t) ≪ |D|^{2(1−2δ)(1−σ)+ε}T^{O(1)}, better than e(σ) exactly for σ > (1−3δ)/(1−2δ). Inputs: δ > 0 (Michel–Venkatesh 2010), δ = 1/1889 (Blomer–Harcos–Michel 2007, as quoted by Yang; its uniformity not checked), δ = 1/120 (Yang, arXiv 2607.04476, July 2026, unrefereed). Ranges: σ > 0.99947, resp. σ > 0.99153; gain about 1% of the exponent there. Classical method ("sup bound × mean value"); new only as a statement.
- **Prop 4 (what an improvement costs).** Improving e(σ₀) through fourth moments forces |L(½,χ)| ≪ |D|^{1/(4(3−2σ₀))} for every χ; at σ₀ = 5/6 that is |D|^{3/16}, Burgess strength, far from the known ¼ − 1/120. Halász–Montgomery at the critical block would need |D|^{1/16} at σ₀ = 5/6 (beyond Weyl).

## Consequences
- RCApps2's almost-all exponent A₀ = 3/2: unchanged. Its binding range is σ ≈ 0.81–0.84. A per-D fourth moment lowers it iff a < 3/4 (A₀ = 1.4919 at a = 0.74, 1.2944 at a = ½), and by Prop 4 any such a implies Burgess-strength subconvexity for all class group L-functions.
- Conditional exponents (1.5997 under Liu's half-plane, ring class PNT range x^{0.6251}): unchanged (they use σ ≤ 0.8749).
- Unconditional least primes (record |D|^455): untouched; nothing here is log-free.

## Why the brief's routes fail (NOGO)
- Halász–Montgomery with the ideal lattice: below |D|^{1/2} the large sieve is an exact Parseval identity on primitive ideals, and the off-diagonal sums are ≈ N^{1/2}|L(½,ψ)|. No Pólya–Vinogradov analogue exists, and there is nothing to subdivide (all ψ ≠ 1 have conductor |D| for fundamental D). An exceptional-set version with the known fourth moment is exactly neutral.
- Huxley-type large values: they need pointwise off-diagonal bounds that equal subconvexity of the needed strength (Prop 4(b)).
- Fourth moment / subconvexity on average: the known bound is a = 1. A Lindelöf fourth moment over this thin family is an averaged form of the Michel–Venkatesh mixing conjecture (known only under GRH: Blomer–Brumley–Khayutin).
- Better or longer mollifiers, moments off the critical line: no gain (NOTE §8).

## Lead (not verified)
One amplifier fits the whole family: λ(p²) = 1 at every inert p for every dihedral form θ_χ, and λ_π(p²) ≤ 1 for every form of level |D| with nebentypus χ_D. An amplified fourth moment of length |D|^θ for that GL₂ family would give a = 1 − θ, twice the subconvexity saving. Reaching the middle range needs θ > ¼, Burgess strength again.

## Novelty check
22 arXiv API queries (class group, ring class or dihedral × fourth moment, moments, zero density, subconvexity, amplification, mixing). Read Yang 2607.04476 in TeX: uniform subconvexity ¼ − 1/120 for class group characters of orders, and a full-family fourth moment q^{1+ε}, but no thin-family fourth moment or density. Nothing found with a per-D fourth moment below |D|^{1+ε}, or a per-D density exponent below e(σ).

## Significance (for the reviewer)
No new theorem worth DISCOVERIES.md. Thm 3 is a routine, microscopic COROLLARIES.md item. Prop 1 and Prop 4 are useful NOGO facts: e(σ) is the end of second-moment technology, and the next step is Burgess-strength input.
