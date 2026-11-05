RECORD: C41b (upper) 2.2195 -> 1.9 (Riemann/Deep/OptConst/sofa_certificate.json)

# Optimization-constants hunt (Tao et al. repository): report

## Result: the ambidextrous moving sofa constant is at most 19/10
- **Constant.** [41b](https://teorth.github.io/optimizationproblems/constants/41b.html): largest area of a
  sofa that can pass a unit-width Z-shaped hallway (a right turn and, far away, a left turn).
  Before: upper 2.2195 on the page (from Baek's announced, unverified proof for 41a; the page for 41a stars it),
  2.37 published (Kallus–Romik 2018, via 41b ≤ 41a), 2√2 (Hammersley-type); lower 1.64495521 (Romik 2016).
- **New.** C_41b ≤ 19/10, unconditional (uses neither Baek's nor Gerver's theorems), for the Z-hallway formulation
  and for the "turn right and turn left from one placement" formulation.
- **Method.** Kallus–Romik's finite-angle containment bound, extended to two hallways (their 2018 paper, end of
  §5, proposes exactly this and leaves it open). New ingredients: a case analysis by orientation lift +
  intermediate value theorem + Gerver's strip bounds, so that arbitrary rigid motions are covered and only "two
  full quadrant windows of hallway angles" are needed; a compact translation domain from connectivity (x-extent
  ≤ W) and clamping; one hallway angle per window (φ = arctan(20/21), rational cos/sin); eight cases, each closed
  by a branch-and-bound whose every leaf bound is recomputed exactly. Proof: [PROOF.md](PROOF.md).
- **Certificate.** [sofa_certificate.json](sofa_certificate.json) (178 KB, sha256
  0f6890c014847f777d4f6413b2873a0952f72a76adcd28dfc035aa3bcdeca6e0): T = 19/10, W = 155/32, 8 cases, 528,486 leaves.
  Check: `python3 verify.py sofa_certificate.json` (standard library only, exact integer/Fraction arithmetic,
  about 4 minutes, exit 0). Output of the final run:

```
self-test ok (exact clipping vs grid sampling, max deviation 0.0001)
pair 01          quadrants [0, 1]:    2800 leaves, max leaf bound 1.899969272386 <= 19/10
pair 02          quadrants [0, 2]:  513503 leaves, max leaf bound 1.899999822499 <= 19/10
pair 03          quadrants [0, 3]:    3006 leaves, max leaf bound 1.897382659004 <= 19/10
pair 12          quadrants [1, 2]:    2626 leaves, max leaf bound 1.899684792473 <= 19/10
pair 13          quadrants [1, 3]:    2100 leaves, max leaf bound 1.899969272386 <= 19/10
pair 23          quadrants [2, 3]:    2787 leaves, max leaf bound 1.899982419468 <= 19/10
same quadrant 1  quadrants [1, 1]:     831 leaves, max leaf bound 1.897720593221 <= 19/10
same quadrant 3  quadrants [3, 3]:     833 leaves, max leaf bound 1.897898899969 <= 19/10
CERTIFIED: every ambidextrous sofa (Z-shaped hallway, or right and left turn from one placement) has area <= 19/10 = 1.9
```
- **Checks done.** Exact areas agree with independent grid sampling (12 random boxes, deviation ≤ 1.2e-4 at grid
  resolution); the same tree with T = 1.85 is rejected; a tree with one flipped bit is rejected; the float search
  (C) and the exact checker agree on every reported maximum to 12 digits. Tightness: for the natural case {0,2}
  there are translations of area 1.8298 (2√2 − 1 ≈ 1.8284 at 45°), so one angle per turn cannot go much below 1.83.
- **Novelty check (5 Nov 2026).** 41b page and README re-downloaded (still 2.2195); repository PRs and issues
  (none on 41b/sofa after #164); Semantic Scholar; GXL evidence search; Zenodo. Closest work: O. Mohammed,
  Zenodo 23667275 (26 Oct 2026, files embargoed; code at github.com/OmarMohammedSTEM/Ambidextrous-Sofa-Problem)
  states "rigorous" values U_6 = 1.6906 … U_30 = 1.6505 for a two-family relaxation, but they come from SciPy
  maximisation (his src/upper_bound.py: "a *lower* bound on U … to claim a proven ceiling we would need
  branch-and-bound"; src/interval_bound.py only checks balls around his optimum) and assume monotone 90° rotation
  paths. They are not upper bounds and are not on the page. Our 1.9 is, as far as we can find, the first certified
  bound below 2.2195.
- **Caveats.** Computer-assisted and unrefereed; proof and code written by an AI agent (Claude); not yet reviewed
  by a human.
- **Next.** Two angles per turn should give about 1.73–1.80 (float search closes 1.80 with 28M leaves; the
  relaxation value is ≥ 1.729): needs a faster exact checker (C with exact integers, or interval arithmetic).

## What was tried
| constant | bound | best known (incl. open PRs) | what we did | outcome |
|---|---|---|---|---|
| 41b ambidextrous sofa | upper | 2.2195 (2.37 published) | two-hallway Kallus–Romik relaxation, exact branch and bound | **1.9 certified** |
| 1d signed 3rd autocorrelation | upper | 1.454727664528751 (PR #312, 81920 cells) | upsample to 163840 cells; p-norm continuation from p0 = 4 (L-BFGS, FFT gradients); local refinement from p = 1.3e5 | no gain: p0 = 4 path heads to ≈ 1.45487; refinement converges to ≈ 1.4547277 |
| 3b Kakeya sum-difference | lower | 1.77898884 (13-point law) | multi-law entropy mixtures: minimax over multipliers ν, global law search at the KKT ν* (s* = 0.2468) | no gain: the 13-point law is the minimax optimum |
| 3c 4-slope sum-difference | lower | 1.674894211307116 (two laws) | KKT multipliers of the two-law certificate, global search for a third law | best other law 1.67349; no gain |
| 3e single-set sum-difference | lower | 1.27155 (simplex limit) | simplex density λ; single type-class laws near the simplex law | λ = 1/2 is optimal; no gain |
| 1c Barnard–Steinerberger | lower | 0.4103779774 (Chen–Gu, arXiv:2610.18888) | analysis of the record structure | remaining gain inside that family ≲ 1e-8; skipped |

Table of all constants and current bounds: .scratch/optconst/table.md; open-PR claims: .scratch/optconst/pr_claims.md.

## Files
- [verify.py](verify.py): exact checker (stdlib only). [sofa_certificate.json](sofa_certificate.json): data.
- [PROOF.md](PROOF.md): the argument. [sofabb.c](sofabb.c) + [make_certificate.py](make_certificate.py): the
  floating-point search that writes the trees (`cc -O2 -o sofabb sofabb.c -lm && python3 make_certificate.py 19 10`,
  about 2 seconds).
