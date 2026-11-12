/-!
# Round 460 — node N11 (GRH families and Siegel zeros), round 16, the unseen Zenodo N11 records: Zeraoulia's normalization defect in Lemma 7.1 of Zhang's 2007 Landau–Siegel preprint is confirmed in the TeX (REFUTATIONS); Daugherty's Lean archive (0.0621/log log q) is clean and states what the node line cites; Jaswal's conditional h ≤ 3000 is already in our class-number chain; Hickey's Lean 'fulcrum' (Siegel zeros give twin primes) is vacuous under the Lean-checked half-plane; six Landau–Siegel 'proofs' are cranks; Weingarten calculus, zonotopes and Gelfand pairs transfer no tool

* Node: N11 (16th round; R426 before). Round type: [check] of the unseen Zenodo records next to the node (Zenodo without a date filter on
  "Siegel zero", "Landau-Siegel", "exceptional zero", lesson 216: 58 records, 46 unseen, 10 on the node), plus three exploration scouts.
* [check] R. Zeraoulia, Zenodo 22865905 (21 Sep 2026), 'A normalization defect and a hidden frequency-density requirement in Lemma 7.1 of
  Zhang's 2007 Landau--Siegel zero preprint': CONFIRMED in the TeX of arXiv 0705.4306v1 (Y. Zhang, 'On the Landau-Siegel Zeros Conjecture',
  v1 only, unpublished). With alpha = (log Q)^(-1) (l. 231, (2.10)), phi_s(x) = Q^(xs), (5.3) phi_s * phi_s' = alpha(phi_s - phi_s')/(s - s'),
  (7.9) tau_1 = K_+(rho,psi)^(-1) phi_0 + h_1 (h_1 in S_1) and (7.10) sum_{T_2} A*(theta)/theta = -K_+(rho,psi), the phi_0 coefficient of
  tau_1 * h_2 is alpha, while the sketch (l. 1704-1765) concludes tau_1 * h_2 = phi_0 - h with h in S; h_1 * h_2 lies in the span of phi_theta,
  theta in T_1 u T_2, and 0 is not in T (approximating phi_0 = 1 by S is the lemma's point, l. 1664-1668). Repair: alpha^(-1) = log Q in (7.10),
  which with A*(theta) << 1 needs an unproved density of frequencies (Zeraoulia). Reviewer: correct (0 not in T), routine: a defective printed
  derivation, not a false lemma (h_1 * h_2 adds no phi_0 term when T_1, T_2 are disjoint; phi_theta * phi_theta = x phi_theta). A broken claimed
  proof: REFUTATIONS. (The conjecture itself
  follows from the Lean-checked half-plane: every zero has beta <= 0.8734, so (1 - beta) log q >= 0.139 for q >= 3.)
* [check] Daugherty, Zenodo 23459046 (qrh_formal v2.0.0, 18 Oct 2026; Lean 4): 29 .lean files, no sorry/admit/native_decide/axiom outside comments;
  littlewood_lower_0621 (Unconditional.lean l. 67-70): q >= 10^100, primitive real chi: L(1, chi) >= 0.0621/log log q, as the node line cites;
  real_zero_gap (SiegelGap.lean l. 36-43): (log 3)/8 <= (1 - beta) log q, from real_zero_le_seven_eighths (the 7/8 theorem, imported): HOLDS.
* [check] Jaswal, Zenodo 23547020 v5 / data 23546913 (21 Oct 2026): h <= 3000 (36,993,853 fields, largest |D| 2,562,211,723) assuming
  Kimura's real-zero bound: KNOWN (discoveries/class-number-5000.md and -10000.md cite it; found by harvest2).
* [check] Hickey, Zenodo 23106948 (2 Oct 2026), 'Twin Primes and Siegel Zeros: a Fulcrum' (Lean 4, standard axioms): Siegel zeros plus an
  'engine hypothesis' give infinitely many twin primes; vacuous now (no Siegel zeros under the half-plane); no node record.
* [check] Zenodo 19503555 (Prop. 1, an ineffective sieve bound for primes in J_n, effective when no Siegel zero exists): effective via the
  half-plane (routine, no explicit constant); Lee 23452649 (local Sato-Tate for a non-CM elliptic curve in growing progressions): off the node,
  not checked. Cranks: Zhu 23558340 and 23598408 (AI-assisted 'proofs', c = 1/2 ... 1/120), Huang 20329369, Nagasawa 21381231, Aksman 22115872,
  Zhai 23496520. Scouts s1-s3 and s5 timed out; s1's question was rerun narrowly (s9) and checked in the TeX by the orchestrator.
* Exploration: Weingarten calculus (Collins-Sniady Gram positivity; N1/N3): needs a probability space with a Haar-observable Gram
  representation of the Weil or Li form, none canonical: no transfer. Zonotopes (Alexandrov-Fenchel on zonotopes built from zero
  coordinates; N13): true convex inequalities with huge slack (+5983 on the first four zeros), insensitive to RH: no transfer. Gelfand pairs
  (Godement-Bochner positivity of spherical functions; N13): no arithmetic kernel maps spherical parameters to zeros without re-encoding the
  explicit formula: no transfer.
-/

namespace Fayed.R460
end Fayed.R460
