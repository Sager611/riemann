NO RECORD: on [0.75, 0.856] every A(σ) record is fixed by mean-value terms (the diagonal of the k-th-power large-value inequality, i.e. the Montgomery cap 3σ − 1 behind 3/(2σ); Heath-Brown's double-zeta term; Guth–Maynard's energy bounds); the half-plane only gives sup-norm bounds for zeta sums, and every channel that uses them returns exactly the old range max(τ_B, 3σ − 1), so σ* = 0.856092 stays and no table entry moves.

# DensityB2 report (2026-11-05)

In words: DensityB lowered A(σ) on (0.8561, B) by feeding the Lean-checked zero-free half-plane Re s > B into the
Halász–Montgomery inequality. This round asked whether the same input, used inside averaged large-value arguments,
lowers the threshold 0.856 or beats any other record on [0.75, 0.856]. It does not. Every record there is limited
by a term that does not depend on how large zeta sums can be: the diagonal of a fourth- (or higher-) power mean
value, or Heath-Brown's double-zeta bound, which is already as small as its diagonal at all binding points. The
half-plane says only that zeta sums of length N never exceed about N^B (and N^{1−B}T^{B−1/2} for long sums), far
above their typical size N^{1/2}; every way of inserting that (pointwise, Hölder, moment interpolation, exponent
pairs, reflection principle, double zeta sums) gives back exactly the range DensityB already found.
Proofs: [PROOF.md](PROOF.md). Certificate: [verify.py](verify.py) (exact rationals, ALL PASS, 0.01 s). Numerical
cross-check: [explore2.py](explore2.py). Paper arguments, not refereed, no Lean.

## Table
| σ range | old record (source) | new bound | derivation (binding input; why the half-plane does not help) |
|---|---|---|---|
| [0.75, 0.76) | 15/(3+5σ) (Guth–Maynard 2024) | unchanged | GM large-value estimate; binding at τ = 6/5 are GM's energy bounds (their Props. 10.1, 11.1; not re-derived here, stated in ANTEDB without any zeta-sum exponent); its double-zeta inputs are diagonal-optimal (Prop. 5); τ_B channel gives 3/(τ_B+1−σ) ≥ 2.31 (verify.py A5) |
| [0.76, 0.7727) | 9/(8σ−2), 15/(13σ−3), 6/(5σ−1) (Ivić 1984, m = 6, 5, 4) | unchanged | cap τ_m = 4σ−2+(2−2σ)/m from the m-th-power diagonal (Prop. 1); start c_m from Heath-Brown's double-zeta term; zeta pair leaves c_m (A4) |
| [0.7727, 0.7735) | 2/(9σ−6) (TTY, improved Bourgain) | unchanged | Jutila k = 3 + Bourgain's theorem (double zeta sums at length N, diagonal-optimal) |
| [0.7735, 0.7777) | 9/(7σ−1) (Ivić, m = 3) | unchanged | as the m-rows |
| [0.7777, 0.7954) | 9/(8(2σ−1)) (TTY, improved Bourgain) | unchanged | as the Bourgain row |
| [0.7954, 0.8) | 3/(2σ) (Bourgain 2000) | unchanged | Jutila k = 2 terms inside Bourgain's theorem; zeta pair in Ivić's lemma only reaches σ ≥ 0.798639 (A4) |
| [0.8, 0.856092] | 3/(2σ) (Ivić 1984, gzd m = 2) | unchanged | Montgomery cap 3σ − 1 = fourth-power diagonal (Prop. 1); half-plane gives τ_B ≤ 3σ − 1 there (Prop. 2) |
| (0.856092, B) | DensityB: 3(B−1/2)/((5/2−B)σ − (5/2−2B)) | unchanged | no channel extends the Montgomery range beyond τ_B (Prop. 2(d)) |

## Answers to the brief
- Records and inputs (feasibility step, done by 22:10 UTC): ANTEDB main (= TTY Table 3) as in the table; Kerr's
  unpublished 3/(2σ) for σ ≥ 23/29 is excluded there. Zeta enters through (i) the diagonal of the k-th power of
  the dual zeta sums (caps τ_k), (ii) Heath-Brown's double zeta sum (starts c_m, Bourgain, Guth–Maynard), (iii)
  exponent pairs in Ivić's lemma (11.74) and zeta-from-exp, (iv) the twelfth moment for Type II zeros.
- Lemma 1 / the pair (B − 1/2, 1/2): (i) is independent of any zeta bound (S(N', W) ≥ |W|N'² holds for every
  set); (ii) is diagonal-optimal at every binding point; (iii) with the zeta pair, Ivić's lemma gives 3/(2σ) only
  from 0.798639, inside Bourgain's range; (iv) has slack on every row. Halász–Montgomery, Corollary mu-lv and the
  reflection principle with the pointwise bounds all give MC exactly for τ < τ_B.
- LV(σ, τ) slightly above τ_B: Lemma hl-improv would need zeta large values on a thin strip below g_B(τ'); by
  reflection these are short zeta sums at levels just below B, and the best available count (Montgomery at
  τ''/2, squared) gives exactly the cap 3σ − 1 again (Prop. 2(d)). Above σ* the same computation shows the
  DensityB range τ_B cannot be extended this way.
- Heath-Brown-type bounds with the twelfth moment: moving his Perron line to Re w = 1 − s and interpolating the
  twelfth moment with the sup bound at B gives a Montgomery condition affine in the interpolation parameter;
  the end points are Heath-Brown's 11σ − 8 and τ_B (Prop. 3(b)). Same for the fourth moment and for Hölder
  between the k = 1 and k = 2 forms of the reflection principle.

## What a new argument would have to supply (PROOF.md §8)
1. A bound for Σ_{t≠t'} |D(t − t')| over a large-value set W (D = dual zeta sum of length T/N) that avoids
   Hölder over all pairs: the fourth moment over W − W is at least |W|(T/N)² from pairs at distance ≤ 2.
2. Or, for the k = 1 reflection principle (its diagonal allows the Montgomery range up to 2σ; Type II would then
   bind): an off-diagonal double-zeta bound for short sums, X^{−1} Σ_{t≠t'} |Σ_{m∼X} m^{−i(t−t')}|² ≪ N^{1−δ} for |W| = N^{2−2σ},
   X = T/N, τ just above 3σ − 1 ("Lindelöf on average over difference sets"). Heath-Brown gives it for
   τ ≤ 5σ − 3, the half-plane for τ < τ_B.
3. Or zeta large values of short sums (length T^{0.28} to T^{0.40}) at levels just below B that beat the squared
   Montgomery bound 2(2 − 2σ'').

## Side findings
- Bug in our DensityB/repro/explore.py: Ivić's zeta-from-exp terms were combined with min instead of max.
  Harmless for DensityB (identical values at σ = 0.80, 0.84, 0.86, 0.87); at σ = 0.78 it produces a spurious
  2.00000 below the table value 2.00893. explore2.py uses the correct max. (File outside this round's write
  scope; the orchestrator may patch line "min(t-6*(s-0.5), 5*t-32*(s-0.5))" to max.)
- Typo in the ANTEDB blueprint, Theorem "Reflection principle" (additive_energy.tex): the printed conclusion
  "ρ ≤ min(2 − 2σ, k(3−4σ)/2 + s' − 1)" is vacuous in the diagonal case; the correct form is
  "ρ ≤ max(2 − 2σ, k(3−4σ)/2 + (s' − α)/2)" (PROOF.md Lemma R; verify.py B7). The ANTEDB Python code does not
  use this theorem. Not reported anywhere (no contact, no posting).
- The half-plane does improve Heath-Brown's double-zeta bound for short sums (τ > 7/4 at ρ ≤ 1/2), only at
  points that bind nothing.

## Numerical cross-check (explore2.py, 22 values of σ)
new − old = 0 for every σ ≤ 0.856 (0.750, 0.755, …, 0.855, 0.856); −0.00407 at 0.858, −0.01008 at 0.860,
−0.02178 at 0.865, −0.02712 at 0.870 (DensityB reproduced). Grid maxima (step 0.004) can undershoot the table by
up to about 0.002; the old/new difference is the test. Values in PROOF.md §9.

## Rerun (about 2 minutes)
`cd Riemann/Deep/DensityB2 && python3 verify.py && python3 explore2.py`
(standard library only; explore2.py reads ../DensityB/repro/explore.py and best_beta.json.)

## Novelty pass (for the record table; this report claims no new bound)
- `python3 .scratch/lst.py 'zero.density|large values|density estimate|Dirichlet polynomial|zero-density|large value' math.NT math.CA --zenodo '"zero density"'`:
  6 title hits in 2026-10/11 (2611.05000, 2610.15989, 2610.18923, 2610.24958, 2610.36885, 2610.37371); abstracts
  read: none improves A(σ) for ζ on [0.75, 0.86] (2610.18923: trilinear sums via Guth–Maynard and Chen–Gupta–Li;
  2610.15989: short strips). Zenodo (15 newest of 156): nothing relevant.
- arXiv API, newest 40 math.NT zero-density papers and 25 "large values" + "Dirichlet polynomial" papers
  (2024–2026): no new A(σ) record on [0.75, 0.86] (2508.02041 near σ = 1; 2507.08296 q-aspect; 2609.22624
  explicit short intervals; 2607.04632 survey).
- ANTEDB main, read 2026-11-05: table as above.
- Our files: DensityB is the only related entry; its claims are unchanged.

## Suggested NOGO line
Half-plane Re s > B inside averaged large-value arguments (Jutila, Heath-Brown, Ivić, Bourgain, Guth–Maynard,
reflection principle, moment interpolation): binding terms on [0.75, 0.856] are mean-value diagonals; no gain
below σ* = 0.856 (DensityB2).

RH relevance: N7 (zero-density exponents: how many zeros can lie off the critical line).

Runtime of this agent: 21:51–22:30 UTC. Downloads deleted.
