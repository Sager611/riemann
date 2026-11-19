AUDIT PASSED — every item on the page holds as written: t_Z([0,1]) ≥ 0.4222354990 follows from arXiv 2610.00386's Theorem thm:cert (re-proved by hand), is rigorously certified (rerun OK) and is independently confirmed (true min A = V + 1.7e-9; all 206 G_i agree to 1e-16); only wording slips, and two dependents (RECORDS.md:8, TREE.md:14) round the bound UP to 0.4222355.

Object: our record (discoveries/integer-chebyshev.md; Riemann/Deep/IntCheb/), built on Theorem thm:cert of arXiv 2610.00386
(v3, 12 Oct 2026, still the newest version). Full chain with line numbers: FULLPROOF.md. Auditor: intchebaudit, 19 Nov 2026.

## Items of the page, one by one
| # | item (page) | status | evidence |
|---|---|---|---|
| 1 | t_Z([0,1]) ≥ 0.4222354990 | HOLDS | exp(V/M) = 0.42223549904149014 (exact V, M); claim is a truncation (rounded down), margin 4.15e-11 |
| 2 | bracket upper 0.4226846976 (the paper's) | HOLDS | upcert rerun: t_Z ≤ 0.422684697504807; 0.4226846976 rounded UP, the right direction for an upper bound |
| 3 | the paper's theorem and 206-polynomial library | HOLDS | thm:cert re-proved (FULLPROOF §1); library = the paper's ancillary library (same 206, same order) |
| 4 | new measure, 1832 cells; new multipliers | HOLDS | K = 1832, all masses > 0, support [0,1], 17 active multipliers; cell description: see §Other |
| 5 | inf A ≥ −0.862192066 on [0,1/2], 20,328 intervals | HOLDS | lowcert sound (FULLPROOF §3), rerun OK; independent global min A = −0.862192064266 (V + 1.73e-9); "20,328" counts bisection nodes, not leaves |
| 6 | exact symmetry x ↦ 1−x | HOLDS | proved for the objects used (cells, masses, library ↔ reflected library, equal c); exact tests by verify.py and independently by sympy; numerically A(x) = A(1−x) |
| 7 | G_i bounds by exact resultants | HOLDS | gcert correct; root-based mpmath G_i (no resultants) agree for all 206 within [4.6e-20, 1.1e-16], never below gcert; min G_i = −0.862191163250 (V binds) |
| 8 | mass condition B < M | HOLDS | exact: M − 1 = 6.78e-17, B − M = −1.0799e-9; it is what the M-normalisation needs (FULLPROOF §2; reviewer CORRECT) |
| 9 | both published certificates re-verified | HOLDS | rerun: published lower ⇒ 0.42222860001915, upper ⇒ 0.422684697504807; published lower also reproduced by independent code: min A = V + 1.0e-8, all 206 G_i within 1.1e-16 of gcert, min G_i = −0.862208282175 as in the TeX |
| 10 | Novelty: nothing at or above 0.42223549 | HOLDS (19 Nov) | see §Novelty |
| 11 | verify.py "about 40 s" | WORDING | 51 s for the new certificate alone, 105 s with both published ones |

## How the theorem gives the number (direction of every step)
t_Z ≥ exp(min(inf_[0,1] A, min_i G_i^res)/M) for any positive μ on [0,1] of mass M with B < M (scaling μ' = μ/M, c' = c/M
puts thm:cert's hypotheses in force with B' < 1; G_i^res = M·G'_i; the middle term −Σ c_j log a_j = A(∞) is needed only
when B' = 1 and would not bind anyway: −0.6897). inf A ≥ V and min G_i^res ≥ gmin are lower bounds, exp is increasing,
M > 0: so t_Z ≥ exp(min(V, gmin)/M) = exp(V/M) since gmin > V. Computed in exact rationals and mpmath interval arithmetic.
The deviation from the TeX statement (no middle term, M ≠ 1) went to a reviewer: CORRECT.

## Defects found (none load-bearing)
1. Rounded UP in two dependents: RECORDS.md:8 ("at least 0.4222355") and TREE.md:14 ("t_Z([0,1]) ≥ 0.4222355").
   The certificate gives 0.42223549904149 < 0.4222355 (short by 9.6e-10), so as stated these lower bounds are unproved.
   Fix (orchestrator; this audit may not edit them): write 0.4222354 or 0.4222354990.
2. NOTE.md/REPORT.md call minleaf − V = 1.2e-11 the "true margin" of the published certificate. It is the slack of the
   bisection's last leaf bound; the true margin is min A − V = 1.0e-8 (published) and 1.73e-9 (ours), by independent search.
3. "20,328 intervals" (page, REPORT.md) and "27,522" count all visited bisection nodes (4096 initial pieces + internal + leaves).
4. verify.py exp_lower returns mp.mpf(iv.a), which rounds to nearest at 53 bits, so the printed 15-digit value is not
   rigorous in its last ulp (≤ 3e-17). The 10-digit claim has a 4.15e-11 margin: no effect.
5. TeX (proof of thm:cert, lines 313-315): "A(∞) = −Σ c_j log a_j" holds only for B = 1 (for B < 1, A(∞) = +∞). Harmless:
   the theorem stays true, the middle term is a superfluous weakening.
6. Falsification note for future checks: plain grid scans of A are unreliable at the 1e-6 level. A 2·10^5 uniform grid
   and an edge-adapted 4.5·10^5 grid (16 lowest minima refined) both missed the true well (width ~1e-6, curvature
   ~2.6e6 from density jumps 1.9e-5 away) and reported min A ≈ V + 6.4e-7; only refining ALL sampled minima found it.

## Dependents
RECORDS.md:8 and TREE.md:14 (both: fix the rounding, item 1). NOGO.md:618 (explore35) and
Riemann/Deep/ExplorePotential/REPORT.md mention integer Chebyshev methods only in general; neither uses the number.
The N10 barrier itself (t_Z([0,1]) > 1/e, so integer-Chebyshev methods cannot reach the prime number theorem) needs only
Pritsker's 0.4213 (2005) or older bounds; the record moves a constant, not the barrier.

## Other checks
- Cell description: each of the paper's 23 support intervals in [0,1/2] holds exactly twice the paper's number of cells
  in our measure (23 of 23; cells.log).
- Hypotheses: support [0,1], all masses > 0, c ≥ 0, 206 polynomials irreducible, primitive, lc ≥ 1, distinct (python-flint
  in verify.py and sympy independently); B < M exact. None violated, in either certificate.
- Memory/time: every process under bash .scratch/memcap.sh 290; peaks 0-87 MB (logs). verify.py rerun 105 s.
- Side remark: thm:cert also holds with g_i = (1/d_i)∫log|Q_i| dμ ≥ G_i (for a library factor H_i = g_i exactly; the
  resultant subtraction only weakens). No gain here: the binding x, x−1, 2x−1 have |Res| = 1 with all 17 active Q_j (resone.py).

## Novelty (19 Nov 2026)
arXiv 2610.00386 is still at v3 (12 Oct 2026; abs page). Monthly listings 2026-10 and 2026-11 of math.NT, CA, CV, SP, PR,
FA, CO, NA, OC, MG, GM, AG (titles matching integer Chebyshev | transfinite diameter | Chebyshev constant): only
2610.00386. arXiv search (3 queries): nothing newer. Zenodo "integer Chebyshev", "transfinite diameter": 0 hits.
Nothing at or above 0.42223549 found; the record stands.

## Files (Riemann/Deep/IntChebAudit/)
FULLPROOF.md (chain with line numbers); indep.py (independent mpmath/numpy code: modes exact, rand, scan, escan, gscan,
G, leaves); localmin.py; gen_inputs.py (verify.py's stdin files); lowcert_dbg_patch.py (prints tight leaves); cmpG.py;
cells.py; resone.py; logs: verify_rerun.log, exact_{new,pub}.log, rand_new.log, scan_new.log, escan_new.log, gscan_{new,pub}.log,
G_{new,pub}.log, gcert_{new,pub}.out, cmpG_{new,pub}.log, tight_leaves_new.txt, localmin_new.log, cells.log, novelty_lst.log.
Rerun: cd Riemann/Deep/IntChebAudit; bash ../../../.scratch/memcap.sh 290 900 uv run --no-project --with mpmath --with gmpy2
--with numpy python3 indep.py ../IntCheb/CERTIFICATE_LOWER_new.json gscan 48  (80 s, 81 MB). Runs on the published
certificate need its ancillary JSON (anc/ in https://arxiv.org/e-print/2610.00386; deleted after use).

RH relevance: N10
