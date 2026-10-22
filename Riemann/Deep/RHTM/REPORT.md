NO DISCOVERY: the 120-state RH machine of arXiv 2609.30306 (Shunia) survives an independent audit; no false step found. Record already listed (bbchallenge Cryptids: RH = BB(120)) and since cut to 119 by the same author.

## What the machine tests

- A_n = n! * prod_{p<=n} p, built by the self-sieve A_J = J*A_{J-1} (J | A_{J-1}) else J^2*A_{J-1}.
- beta_n(B) = number of steps B -> floor((B-1)/n) until 0; Delta_n = max(n-1-beta_n(A_n-1), 0); q_n = floor(log_4 3n).
- Claim: RH <=> Delta_n < (q_n+1)*2^q_n for every n >= 2; the machine tests n = 2,3,4,... and halts at the first failure.
- All arithmetic is exact integer arithmetic on unary counters: no approximations of e, gamma, log or H_n, so no rounding-direction bug is possible.

## Math audit (hand re-derivation + exact checks)

- Forward: Schoenfeld 1976 Thm 10 is cited correctly (RH => |theta(x)-x| < sqrt(x) log^2 x/(8 pi), x >= 599). Slack at large n is the factor 4 pi/log 4 ~ 9.
- Small range 2 <= n <= 2656: the table is right. The tightest row is q=5: 1365/L = 189.6 < 192. Exact big-integer check: the criterion holds on the whole range, with max Delta/threshold = 0.0375.
- Converse: a one-sided bound theta(x)-x >= -C sqrt(x) log^2 x gives RH by Landau's positivity theorem. Standard, and correct as written.
- Clock lemma checked exhaustively to n = 10^6. Float bound to n = 2*10^6: Delta <= 151 against a threshold of at least 24,576. Machine RH-checks are hopeless at reachable n.

## Machine audit (own C simulator evidence/tmsim.c, own interpreter evidence/ref.py)

- The paper's appendix code rebuilds the printed table byte for byte (sha256 b4b0f076...).
- Bootstrap reproduced exactly: 89,775,610 steps, head at -3818, exact 1912-cell support.
- Blank-tape run of the literal 120-state table, stages n = 2..5 (5.03e9 steps): every inc/dec with its outcome matches the interpreter's trace (hash-identical). Registers were decoded from the tape at every op for n <= 4 and sampled for n = 5. Stage 6 was not finished (CPU starved).
- Falsified instances on the literal table: the tape was edited at phase starts, 123 cases. Width cases sit at threshold-1 and threshold for q = 1..6, up to n = 2656. Digit cases use fake A_n with Delta at the threshold boundary. Compare and seed phases were also covered. The machine halts exactly when the formula says, only through the real halt (77,1), and its op trace matches the interpreter.
- Falsified variant: deficit inflated to 2(n-1)+1 and recompiled with the paper's builder (123 states). It halts in stage n = 5 at step 5,032,616,178, as predicted.
- Reductions 122 -> 121 -> 120: own homomorphism check plus own window fixpoint on the raw table. The only disagreements are (prep_1,1) and (6.continue.2,0). The window excludes (6.continue.1,0) and (6.continue.2,0). Exclusion of (prep_1,1) rests on the paper's normal-form argument, not on a finite certificate; it held on every run.

## Record

- Prior: 744 states (Matiyasevich-O'Rear 2016). aimir/tiny-riemann has a Lean-verified 275-state machine using the 744 criterion.
- Author's repo jshunia/rhtm has RH_119.tm (merges states 6 and 37). I checked that it equals the 120->119 projection, then reran stages 2..5 and 27 tape-edit tests on it: all pass.
- No smaller RH machine found (arXiv API, GitHub search, bbchallenge Cryptids).
- Nothing clearly smaller is in reach. Window invariants of radius 1..6 on the 120 table exclude only (37,0), so no further merges come from them. Going below 119 needs a new source or backend plus a full proof.

## Caveats

- This is testing, not proof. The unbounded parts (register microcode for all values, unreachability of (prep_1,1)) remain the paper's hand proofs.
- Details and rerun steps: NOTE.md. Evidence: evidence/ (code, tables, logs, 128 KB).
