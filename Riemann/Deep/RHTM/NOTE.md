# RHTM audit notes (arXiv 2609.30306, 120-state RH machine)

## Chain checked

1. Criterion (paper Thm 4.6): every step was re-derived by hand. Lemma beta: M_d = b+...+b^d, and l_b-1 <= beta_b <= l_b. Lemma center: 1 <= a_n <= 1+log n; A_n is not a power of n for n >= 3. Prop large: Delta_n < -(E+a_n)/log n < sqrt(n) log n/(8 pi) < 0.11 (q+1) 2^q. Prop small: the table rows check out. Converse: Landau with g = theta-x+C sqrt(x) log^2 x >= 0 has Mellin transform P(s)/s - 1/(s-1) + 2C/(s-1/2)^3. No error found.
2. Source (rh120_source = rh_source + width regrouping). Semantics were read by hand from Wade's builder:
   - (decnz, X) runs X only on the positive outcome; the zero outcome pops one pair.
   - break.k and continue.k erase PC bits.
   - The macros give exactly the paper's lemmas: transform x(j+v+1)+j; divide (q, r, J-1-r) with split divisor; halve; seed; compare halts iff f > q.
3. Backend: Wade TMBuilder, PC word on the left, unary registers U(r) = 1^(r+1)0 from cell 5. Each successful dec shifts the tail left and appends one zero block, so the tail grows forever. That is harmless.
4. Dispatcher: a four-zeros scanner. The paper's certificate gives a max zero-run of 3. Every run here agreed with it.
5. Reductions: bisimulation (merges 2b.reg.-2.inc with inc.shift_2), then a projection merging N24 and 6.continue.2 (needs (6.continue.2,0) unreachable; window-certified). 119 = 120 + merge of prep_1 with 6.continue.1 (needs (prep_1,1) unreachable [normal form] and (6.continue.1,0) unreachable [window]).

## Tests (evidence/logs)

- natural_120.log, natural_119.log: blank-tape stages 2-5, trace hashes equal to ref.py. The paper's own literal checks stop at n=4.
- surgery_small_120.log (109 cases) and surgery_big_120.log (14 cases): tape edits at the starts of the width, digits, compare and seed phases. They cover threshold boundaries for q = 1..6, the q-boundary values of n, and fake A_n.
- variant_inflate.log: a falsified source compiled with the paper's own builder halts exactly at the predicted stage.
- reductions.log, mathcheck.log, pipeline.log.
- Not done: stage 6 natural run (357,669 ops, about 2e11 steps). It was CPU-starved at load 25-50 and stopped.

## Novelty / record

- bbchallenge Cryptids lists RH = BB(744) (2016) and BB(120) (this paper).
- arXiv API: there is no other RH Turing machine paper after Aaronson-Yedidia 2016 (whose RH machine has 5,372 states). The 744-state machine is Matiyasevich-O'Rear, not Aaronson-Yedidia as the brief said.
- GitHub: jshunia/rhtm (RH_119.tm, note dated 2026-09-28); aimir/tiny-riemann (275 states, Lean-verified, uses the 744 criterion).

## Rerun

1. Get the source: `curl -sL -o src.tgz https://arxiv.org/e-print/2609.30306` and untar it into src/.
2. Extract the code: `python3 extract_appendix.py src/rhtm.tex`.
3. Copy b08 to vendor/builder.py; b07 to dispatch.py; b04, b05, b06 to baseline_source.py, rh_source.py, rh120_source.py; b09 to minimize.py; b10 to tables.py; b11 to window_certificate.py; b12 to compatible_merge.py.
4. Rebuild the tables: `python3 build_only.py`. This writes machines/RH_121_raw_named.tm and the names maps.
5. Simulate: `cc -O2 -o tmsim tmsim.c`, then `./tmsim RH_120.tm roles120.txt --maxstage 5 --check 13`.
6. Tape-edit tests: make snapshots with `--save K snap_<phase>.bin`, using K = 49 for digits, 76 for width and 86 for compare. Then run `python3 tests.py RH_120.tm roles120.txt small`.
7. Other checks: `python3 reductions.py` (paths are inside the script), `python3 mathcheck.py 2000000`, and `python3 variant_build.py` run in the builder directory.
