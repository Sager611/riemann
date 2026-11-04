NO FLAW FOUND: all 5 sections of OpenAI's "The circulant Hadamard conjecture" (23 Sep 2026, TeX + PDF from github.com/openai/math) checked step by step by hand and by exact computation (4 scripts, OUTPUT.txt); no false step, gap, case gap or quantifier slip. A Lean proof of the full statement exists in openai/math (static scan clean; not compiled here).

# Audit: OpenAI's proof of Ryser's circulant Hadamard conjecture (agent "ryser", 2026-11-04)
Source: openai/math preprints/The-circulant-Hadamard-conjecture-September-23-2026 (paper.pdf, build/sections/*.tex). It is the "corrected proof" linked from Morris's arXiv 2302.08346v2 (withdrawn 30 Oct 2026: "Error found by AI").

## The argument
1. Prop 2.2: order n > 1 forces n = 4u^2, u odd (Turyn's 2-adic descent in group-ring form).
2. Lemma 3.1 (new core): if F*K = p^e*b in A[C_{p^e}] (A unramified local at p, b a unit), then F(rho)/F(1) is a unit = 1 mod (rho-1), rho of order p. In C_{u^2} the p-part of the group order equals v_p(u^2), so for xx* = u^2 the alternating product Delta(x) = prod_S x_S^{(-1)^|S|} over the characters of order dividing rad(u) is a unit at every prime with |sigma(Delta)| = 1 (exponents sum to 0): a root of unity of p-power order for every p | u (Prop 3.4).
3. Section 4: c, d, g = (row at Z = 1, -1, i)/2 have norm u^2; R = Delta(d)/Delta(c), W = Delta(g)/Delta(c) are odd-order roots = 1 mod (1+i), so R = W = 1; but the +-1 entries force lambda_2(R) + lambda_(1+i)(W) = u^2/c_empty != 0 in the residue field at 2.

## Verdicts by hand
- Lemma 2.1 (unramified local ring, basis 1..xi^(D-1), valuation <-> coefficient divisibility): standard, correct.
- Prop 2.2: correct; indices k = 2s-j, t = s-j stay in range; final count 2^(s+1) odd squares = 0 mod 8 vs 4u^2 = 4 mod 8.
- Lemma 3.1: correct. e >= 2: valuations at a primitive p^e-th root sum to e, so one factor has v >= 1; its remainder mod Phi_{p^e} has coefficients in pA and Phi_{p^e} projects to p, so q(F)/p is integral; both ratios survive projection. e = 1: one augmentation is a unit. Remark 3.2 (sharpness) is right.
- Prop 3.4: correct. Every partial evaluation satisfies the lemma's hypothesis (b = u^2/p^e is a unit); pairing S with S+{p} gives residue 1 at all primes over p; units elsewhere; intersection of all B_m is B; Kronecker; residue-1 roots of unity in residue characteristic p have p-power order.
- Section 4: d - c = -2b and g - c = (1+i)(ib - J - (1-i)U) are exact; lambda_t is a homomorphism; J_S = 0 for S nonempty. Correct.
- Section 5 (Barker): even-length reduction classical and correct; odd lengths cited (Turyn-Storer, Schmidt-Willms).

## Computation (exact integer/fraction arithmetic, single thread, < 70 MB)
- check_comparison.py: all 108 F in Z[C_9] with FF* = 9 (18 of type 3g, 18 of type (3-2J_3)g, 72 six-term) pass; 1300 random pairs F*K = p^e*b over Z and Z[i] in C_9, C_27, C_25, C_49 (v(F(1)) from 0 to e) pass; the conclusion fails, as it must, for p-2J in C_p, p^e-2J (scalar exponent too big) and 3-2J_H in C_3 x C_3 (not cyclic).
- check_alternating.py: every integral x with xx* = u^2 and order-only character values has Delta = 1 (u = 3, 5, 9, 15, 21; 1152 cases at u = 15). Control with group C_rad(u): Delta in {-1, +-i} occurs (384 of 576 at C_15), so the hypothesis matters.
- check_binary.py: random sign arrays at n = 36, 100, 196, 324, 900 (full alternating product, one prime over 2, also where 2 splits): identity lambda_2(R) + lambda_(1+i)(W) = u^2/c_empty holds 510/510; R = W = 1 never; order 4 gives W = i (paper's Remark 4.2).
- check_barker.py: table rows are Barker; exhaustive lengths <= 20 give {2,3,4,5,7,11,13}; P(t) = n mod 4 for n <= 16; no circulant Hadamard row of order 16.

## Does the argument prove something false? Not on the near-misses tried.
Menon-Hadamard difference set in C_4 x C_3 x C_3 (a 12 x 3 perfect binary array, built and verified): all norms hold, but the odd part is not cyclic. With K1..K4 the subgroups of order 3, c = J - K1 - K2 - K3 + K4 has c(chi)/c(1) = -1 for three of the four character pairs (d(chi)/d(1) = -1 for the fourth), so Delta has even order and R or W is -1 or +-i for every order-3 character, while the binary identity still holds. The chain breaks exactly at Lemma 3.1, the only place where cyclicity is used. Circulant weighing matrices gave no applicable test: the lemma needs the weight's p-part to equal the order's p-part, lifted examples are invisible to order-p characters, and I know no proper example of that shape.

## Prior review and Lean
- openai/math: issues and discussions disabled, PR conversations locked; no response found on arXiv (search, Nov 2026).
- Lean: lean/OAI/LinearAlgebra/CirculantHadamard (68 files, Lean v4.34.1) proves OAI.CirculantHadamard.exists_iff_order_one_or_four. The Comparator challenge states exactly "a real +-1 circulant H with H*H^T = n*I exists iff n is 1 or 4" (Fin subtraction is cyclic); permitted axioms propext, Quot.sound, Classical.choice. rg finds no sorry, admit, axiom, native_decide or implemented_by; imports are Mathlib, its own folder and omega only. Not compiled here (Mathlib needs several GB; brief caps 300 MB). Last check worth running elsewhere: in openai/math/lean, lake exe cache get && lake env comparator ComparatorChallenges/CirculantHadamard.json.

Verdict: the proof looks correct. If the Comparator run passes, Ryser's circulant Hadamard conjecture (and the even half of the Barker conjecture, also formalized) is settled; external result, not a discovery of this project.

