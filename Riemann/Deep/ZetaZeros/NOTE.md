# NOTE — audit of arXiv 2610.08965 v1 (Knausgård, 6 Oct 2026)

Claim: liminf N_d/N ≥ 1645064/1960733 = 0.839005 (distinct), liminf N_0^s/N ≥ 1669159/2478195 = 0.673538
(simple AND on the line), unconditional; plus 0.8368 / 0.8892 / 0.8893 (Thms 1.4, 1.5).
Sources: e-print TeX `knausgard839.tex` (2584 lines) + ancillary files (anc/, 354 files, all SHA-256 OK;
GitHub kristianmk/839 = arXiv copy + one README link).

## 0. Context the brief missed
"2/3 simple needs RH" is out of date. Unconditional 2/3-type results exist since 2026:
Alpöge–Furman arXiv 2608.13637, Lamzouri arXiv 2609.02882, Lean library zeta23
(anthropics/formal-math @ fbdc36b; README: hypothesis-free headline theorems, only propext/choice/Quot.sound,
Comparator-checked): 0.67250 simple-on-line, 0.83625 distinct. Our own NOGO.md l.162 / TREE.md l.47 already
record this. The paper is an increment (+0.0010 / +0.0020) on that framework.

## 1. Energy input (task 1) — no RH
- Paper Prop. 3.1 (TeX l.814–846) quotes BGST Lemma 5 = Baluyot–Goldston–Suriajaya–Turnage-Butterbaugh,
  Acta Arith. 214 (2024), arXiv 2306.04799. Checked in BGST TeX: F(x,T) = Σ_{ρ,ρ'} x^{ρ−ρ'} w(ρ−ρ'),
  w(u)=4/(4−u²), sum over ALL zeros 0<γ,γ'≤T with multiplicity; Thm 1: F(α) = T^{−2α}(log T+O(1)) + α
  + O(1/√log T) uniformly 0≤α≤1, no hypothesis. x^{ρ−ρ'} = x^{β−β'} x^{i(γ−γ')}: NOT unimodular, NOT
  ordinate-only. Lemma 5 (test fn r even, supp [−1,1], Lipschitz at 0) is exactly what the paper uses.
- Paper's operator (l.770–792): z_ρ = i(ρ−½)L/2π (real iff β=½), F_z(u)=η_ε(u)e^{−2πiuz},
  A = Σ ν_z F_z⟨F_{z̄},·⟩. Checked: A Hermitian (zeros symmetric under ρ↦1−ρ̄ ⇔ z↦z̄), tr A = N,
  tr A² = Σ ν_zν_w K_ε(z−w)² with COMPLEX arguments.
- Weight removal: (R − R''/4L²)^ (x) = (1+π²x²/L²)K_ε(x)²; at x = z_ρ−z_ρ' factor = 1−(ρ−ρ')²/4 = 1/w. Checked.
  Main term R(0)+2∫αR = ∫f²+∬|s−t|ff = C(f). Checked (also closed form for cosine windows, see consts.py).
- In Lean the bandwidth λ<1 is used, so BGST's zero-free-region error x^{1−2η} needs nothing deep.

## 2. Where β = ½ enters (task 2)
Not from ordinates. Off-line pair {ρ,1−ρ̄} = points z, z̄ (Im z = (β−½)L/2π ≠ 0) gives the term
ν(vw*+wv*) = (ν/2)((v+w)(v+w)* − (v−w)(v−w)*): at most ONE positive eigenvalue, trace 2ν.
Lemma 2.1 (threshold, l.524–541; proof re-derived, correct) with p = #(such terms):
tr A² ≥ 2c tr A − 2c tr P + tr φ_c(P) − c²p. At c=2 (line, l.1072–1083): N−n = tr Q ≥ 2p ⇒
tr A² ≥ 2(N−n) + tr φ_2(U): every off-line zero (or multiple zero) is charged like a double zero.
This is the AF/Lam mechanism; correct for ANY β≠½, however close to ½.

## 3. Other written steps re-derived (all correct)
pinching (Jensen), Lemma maj (Selberg majorant; ĝ supported in [−θ,θ], Λ = 2γ₀(1+sinc θ)/θ),
Lemma pair, Lemma transfer (budget, rewards −6x₂, telescoping Φ(σ)−Φ(σ'), checked identity),
Thm gram (matching ⇒ S θ-separated; G_SS ⪯ cI so no clipping), Thm general (μ,ν,r_h,r_k inequalities),
Thm line-general, counts N = n₁+2n₂+Σj+Σ2j, N_d = n₁+n₂+h+2k.


## 4. Constants (scripts/consts.py)
C(f): 1.33270495828 (cap 1.332705, margin 4.2e-8); 2nd 1.32794727284 (cap 1.3279474, 1.3e-7);
3rd 1.33246363317 (cap 1.3324638, 1.7e-7); closed form = quadrature to 5e-12.
min(g−f) on I: 5.0e-4 (1st, at u=½), 5.2e-2 (2nd), 2.0e-3 (3rd). Λ: 1.67105 (<1.675), 1.78943 (<2), 1.67931 (<2).
min K² on [0,θ]: 0.065799 ≥ β=0.065747; 0.051086 ≥ 0.051067; 0.066829 ≥ 0.066807.
μ,ν,r_h,r_k = 1.960733, 3.684917, 1.792984, 0.0688 (all ≥0). q, s₀ exact as fractions.
Further counts (Thm 6.2 a,b,c) re-derived; t=27/8: 4t−t² = 2.109375 ≥ 2+2x₁ = 2.1087.

## 5. Numerical searches
localmin.py (random starts + pattern search, all labellings): first 9.750e-5 at 7 simple points,
gaps ≈1.025 (paper: ~1e-4 there); second 1.265e-5 at gaps (1.03,1.04,1.97,1.04,1.97,1.04,1.03)
(paper: 1.3e-5); third 2.0e-5 (paper: 1.0e-5 and 2.0e-5; kinks at knots make local search incomplete).
groundstate.py (periodic labelled configurations, chunk difference removes boundary terms, clipping by
eigenvalues): min defect/point +2.81e-3 (distinct, e.g. labels 111221 gaps 1.03,1.87,1.07,1.96,1.05,1.01),
+3.29e-4 (line, gaps ≈1.04/1.97 mixtures). Nothing below the certificates.
opcheck.py: operator layer with off-line pairs, 400 random configurations, all identities and Lemma 2.1 hold.

## 6. Lean and code
axioms.txt: each headline theorem uses propext, Classical.choice, Quot.sound + one
*.check_passes._native.native_decide.ax_1_1. rg for implemented_by|extern|unsafe|csimp|axiom|sorry|elab|macro:
only partial defs in RunCheck.lean (diagnostic counters, not in the proof). Defs.lean tables = inputs.json =
paper tables (checked by eye for window 1). Checker compiled via precompileModules (standard native_decide trust).
upstream_relax.patch: 14 zeta23 files; b ≥ ½ → ≥ 1/8, a ≥ ½ → ≥ ¼, constants 4π→16π, 8π→32π in o(·) bounds.
Not done: Lean build (forbidden), re-running the 3 searches (authors report C++/Python/Lean counters agree).

## Scripts
scripts/*.py expect the e-print unpacked in ./src (curl -s -L -o e.tgz https://arxiv.org/e-print/2610.08965;
mkdir src; tar xzf e.tgz -C src). Run with nice, single thread; each < 5 min, < 60 MB.
