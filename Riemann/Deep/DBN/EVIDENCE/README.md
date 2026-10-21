COMPLETE: Λ ≤ 0.148408913066829 is backed by Planat's unchanged producer and auditor (their files match his arXiv source byte for byte) plus separate exact certificates for the new kernels; Λ ≤ 0.148368783506957 is backed by a superset of his auditor (all his checks kept, mixtures and kernel geometry added). Both are conditional on Planat's premises A1–A3, exactly as his 0.148612 is.

# Evidence package for the dbn4 bounds (2026-10-21)
"Planat's auditor" is the exact-arithmetic checker he ships with arXiv:2610.18214: it replays a step-by-step history of the comparison argument and rejects any step that does not follow. "Kernel" means one inequality 2K_y ≥ λ(y/h)J that turns a source floor into a contraction rate; "mixture" is the same with two shifts at once.

## The four gaps the reviewer named, and what closes each
| Gap | Evidence | Result |
|---|---|---|
| 1. Floors byte-identical? | `02_floor_identity.txt` | Each evidence file is Planat's refined source with one insertion of 99661 bytes at byte 7,465, right before the closing bracket of his kernel list; the insertion parses as exactly the 786 new kernels; no original byte is removed or changed. All 6,080 floors and 608 heights per provider are identical. Every floor, height and own-pair cap consumed by both histories matches his original files (0 mismatches). Their sha256 values equal the hashes in both history headers. |
| 2. Kernel certificates complete? | `04_*`, `04a_identities.txt`, `certs/single_leaves.txt` | All 843 single kernels (his 57 + our 786): two independent polynomial constructions agree with his `check_geometry.poly`; the identity 2K_y/y − λJ = 4P/(Q*·D4) holds exactly at 2,529 rational points; the 9,944 leaf boxes tile the unit cube exactly (volume 1, no overlaps); every Bernstein coefficient, recomputed on each box straight from the polynomial, is ≥ 0. Denominators: Q* > 0 away from the zero's own pair, D4 > 0 because ρ < 1, all scalings are positive. His unchanged `prove()` passes 830 kernels; the other 13 need depth 34 against his cap of 30 and pass with only that number raised. 22 symbolic identities (scaling, evenness, the reduction, the homogenised blocks, own pair, ODE step, terminal time) all reduce to 0. |
| 3. What his unchanged auditor checks on the new kernels | `03_planat_auditor_148409.txt` | A fresh run of his unchanged producer reproduces the 0.148409 history byte for byte (sha256 e8b13ddc…). All 18,997 steps use new kernels (751 distinct; none of his 57). The table lists all 28 checks in his auditor with line numbers and call counts. Eligibility, validity, caps, rates, the outward ODE step and the dyadic lattice run on every step; the terminal, strip-cost and start checks run once each. 0 failures. His 9 negative controls and 12 more aimed at new-kernel steps are all rejected. A forged history with an invalid kernel (λ = 21/10) and consistent hashes PASSES his auditor: kernel geometry is outside it by design ("kernel inequalities have separate certificates"). That is why gap 2 matters; the certificates reject that kernel. |
| 4. The 0.148369 variant's own auditor | `05_*`, `06_*`, `certs/mix_leaves.txt`, `scripts/audit_superset.py` | The old custom auditor (`dbn4/audit_mix.py`) skipped several of Planat's checks and let a negative index wrap around (now fixed, and superseded). It is replaced by `audit_superset.py`, generated from Planat's verified file by exact edits. The diff with indentation ignored removes 3 lines, and each is replaced by a widened version (signature, schema, return dict). All other lines are his. It adds a mixture branch (same cell, provider and height; both floors and caps from the source; declared distinct shifts; λ₁, λ₂ ≥ 0; domain; exact rate) and checks every kernel and mixture used against its leaf certificate. With certificates on it passes Planat's own 0.148612 history, the 0.148409 history and the 0.148369 history (9,012 mixture steps; records and mixture table unchanged, only his header fields added). It rejects the forged λ history, his 9 controls, and 34 further damaged histories. One of these raises λ₁ with every rate kept consistent; it passes only when the certificates are switched off. All 345 mixtures: three polynomial constructions agree (including the producer's own); 1,035 exact point identities hold; the 6,160 boxes tile the cube with all coefficients ≥ 0. |

## Files
- `01_planat_inputs_sha256.txt`: the arXiv tarball (fresh download, sha256 4599b663…) and all 54 entries of Planat's SHA256SUMS.json. arXiv stores 8 members unzipped; re-gzipping them as his code does (level 9, mtime 0) restores his hashes exactly.
- `02_floor_identity.txt`, `03_planat_auditor_148409.txt`, `04_single_kernels_summary.txt`, `04_single_kernels.tsv` (one row per kernel), `04a_identities.txt`, `05_mixtures_summary.txt`, `05_mixtures.tsv` (one row per mixture), `06_superset_auditor.txt`, `06_auditor_diff.txt`.
- `certs/single_leaves.txt`, `certs/mix_leaves.txt`: the certificates (leaf boxes). The checker never trusts them: it rebuilds each polynomial from the hashed source or history.
- `scripts/`: `common.py` (paths), `leafcert.py` (polynomials, prover, certificate checker), `s1`–`s6` (the steps above), `derive_identities.py`, and `audit_superset.py` (generated by `s6`).

## Exact commands (repo root; Python ≥ 3.10; standard library only except the last step)
```
mkdir -p .scratch/dbn4v/dl
curl -sSL -o .scratch/dbn4v/dl/e-print_2610.18214v1 https://arxiv.org/e-print/2610.18214v1
cd Riemann/Deep/DBN/EVIDENCE/scripts && export PYTHONDONTWRITEBYTECODE=1
nice -n 10 python3 s1_inputs.py          # verify + extract Planat's source (5 s)
nice -n 10 python3 s2_floor_identity.py  # rebuild the two evidence sources, byte diff (2 min)
nice -n 10 python3 s3_planat_run.py      # his unchanged producer + auditor, coverage, controls (2-9 min)
nice -n 10 python3 s4_kernel_certs.py    # 843 kernel certificates (3 min)
nice -n 10 python3 s5_mix_certs.py       # 345 mixture certificates (2.5 min)
nice -n 10 python3 s6_superset_run.py    # superset auditor on 3 histories, 44 controls (1.5 min)
<python with sympy> derive_identities.py # symbolic identities (15 s)
```
Each script sets its own hard time limit (SIGALRM) and writes its EVIDENCE file. Scratch goes to `.scratch/dbn4v` (about 70 MB, mostly the unpacked arXiv source). Measured here: about 10 CPU-minutes in total, one niced process at a time; peak memory 60 MB for the steps above, except 375 MB for `s6`, which holds three histories at once.

## Not covered here
The premises A1–A3 themselves: Raymond's source floors, the wall and first contact. REPORT4 judged them proved, but they are computer-assisted and unrefereed. Planat's native Arb runs that produced his refined floors are also not re-done; we only show his files are used unchanged. Both bounds are conditional implications in exactly his sense.
