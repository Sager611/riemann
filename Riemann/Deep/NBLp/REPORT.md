THEOREM: For 1<p<2, ζ(s)≠0 on Re s>1/p ⟺ 1_(0,1) lies in the L^p(0,1)-closure of span{ρ_{1/k}: k≥2} (⟺ that span is dense in the subspace S_p of functions constant on each (1/(n+1),1/n]; literal density in L^p(0,1) is always false). Hence unconditionally for every 1<p≤8/7; for 1<p<8/7 the explicit approximants 1+Σ_{k≤N}μ(k)(1−k/N)ρ_{1/k} have L^p norm ≪ N^{−η} for every η<1/p−7/8.

Node N5 (Nyman–Beurling / Báez-Duarte). Full proof: [PROOF.md](PROOF.md). Agent `nblp`, 2026-11-05. Paper proof only; no Lean.

## What the brief asked, and the correction

The brief asked for density of `span{ρ_{1/k}}` in `L^p(0,1)`. That is false for every `p`: on `x ∈ (1/(n+1), 1/n]`
one has `ρ_{1/k}(x) = {n/k}` (PROOF.md, Prop. C(a)), so the span lies in the closed, infinite-codimensional
subspace `S_p` of step functions on those intervals. The right statement, and the one Báez-Duarte's criterion
uses, is `1_(0,1) ∈ closure` (`d_{N,p} → 0`), which is equivalent to density in `S_p` (Prop. C(d)).

## Results (PROOF.md `1)

* **Theorem A** (`1 < p < 2`): `ζ ≠ 0` on `Re s > 1/p` ⟺ `d_{N,p} → 0` ⟺ span dense in `S_p`. Under the zero-free
  hypothesis, `‖1 + Σ_{k≤N} μ(k) k^{−ε}(1−k/N) ρ_{1/k}‖_p ≤ A_p ε^{1−1/p} + 2ε` for `N ≥ N_0(p,ε)`.
* **Theorem B** (open range, explicit): if `ζ ≠ 0` on `Re s > θ` and `p < 1/θ`, then
  `‖1 + Σ_{k≤N} μ(k)(1−k/N) ρ_{1/k}‖_p ≤ C N^{−η}` for every `η < 1/p − θ`.
* **Corollaries.** Unconditional for `1 < p ≤ 8/7` (from the Lean-checked 7/8 half-plane). With `Θ = sup Re ρ`:
  `1_(0,1) ∈ closure` ⟺ `p ≤ 1/Θ` (for `p < 2`); never for `p > 2`. RH ⟺ the criterion holds for all
  `p ∈ (1,2)`; under RH the approximants converge at a power rate `N^{−(1/p−1/2)+}`, unlike `d_{N,2} ≍ (log N)^{−1/2}`.

## Proof summary

1. *Necessity* (easy): a zero `w_0` with `Re w_0 > 1/p` gives the bounded functional `F ↦ ∫_0^1 F x^{w_0−1}`, which
   kills every `ρ_{1/k}` (Mellin transform `−ζ(s)(k^{−s} − k^{−1})/s`) but not `1_(0,1)`.
2. *Hölder bridge* (Lemma M): `‖F‖_{L^p(0,1)} ≤ (1 − r(a−½))^{−1/r} ‖x^{a−½}F‖_{L²}`, `1/p = ½ + 1/r`, `a < 1/p`.
   Weighted L² is Mellin–Plancherel on the line `Re s = a`.
3. *Open range* (Theorem B): on a line `θ < a < 1/p`, the Riesz–Perron formula plus the Lindelöf-type bound
   `ζ^{±1} ≪ |t|^{η'}` inside a zero-free half-plane (Lemma A: Borel–Carathéodory + three circles) give
   `|1/ζ − P_N| ≪ N^{−η}|t|^{η'}`, hence a power rate.
4. *Endpoint* (Theorem A, the new step): put `H_ε = Σ μ(k)k^{−ε}ρ_{1/k}` (as a function: `−E_ε([1/x])`,
   `E_ε(m) = Σ_{j≤m} Σ_{d|j} μ(d)d^{−ε} − m/ζ(1+ε)`). Work on the line `Re s = σ' = 1/p − ε/2`, mirror image of
   `1/p + ε/2` across the line `Re s = 1/p` that bounds the zeros. The Hadamard product gives
   `|ζ(s)| ≤ |ζ(s+ε)|·|γ(s+ε)/γ(s)|` there (every factor `|ρ−s−ε|/|ρ−s| ≥ 1`). Expanding
   `|1 − R|² = 1 − 2Re R + |R|²` (`R = ζ(s)/ζ(s+ε)`) and evaluating `∫ R/|s|² dt` exactly by residues
   (`= 2π[R(2σ')/(2σ') − 1/((2σ'−1)ζ(1+ε))]`) gives `‖x^{σ'−½}(1 + H_ε)‖²_{L²} ≤ K_p ε` (up to a `ρ_1` term of
   size `ε`). Hölder costs only `ε^{−1/r}`, so `‖1 + H_ε‖_p ≤ A_p ε^{1−1/p} + ε → 0` because `p > 1`. Finally
   `H_ε` is the L^p limit of the explicit combinations `Σ_{k≤N} μ(k)k^{−ε}(1−k/N)ρ_{1/k}` (step 3 on the line
   `1/p − ε/4`).

## Numerical check (one process, M = 10⁹ terms, exact integer sieve; code `nb.c`, log `run1e9.txt`)

`‖1 + Σ_{k≤N} μ(k)(1−k/N) ρ_{1/k}‖_{L^p(0,1)}` (i.e. `c_k = −μ(k)(1−k/N)`), computed as
`(Σ_m |v(m)|^p/(m(m+1)))^{1/p}`, `v(m) = 1 + Σ_k μ(k)(1−k/N){m/k}`; interval = partial sum to `M` plus the
rigorous tail bound `max|v|^p/(M+1)`:

| N | p = 1.1 | p = 2 (reference) |
|---|---|---|
| 10 | [0.3459480954, 0.3459480986] (exact via period 2520: 0.34594809800) | 0.36463013 |
| 100 | [0.0592261752, 0.0592262243] | [0.1326155, 0.1326187] |
| 1000 | [0.0150505588, 0.0150512962] | [0.1206823, 0.1210579] |

At `p = 1.1` the norm falls 5.8× and 3.9× per decade (local exponents 0.77, 0.59), consistent with the
RH-conditional rate `N^{−0.409+}`; at `p = 2` it stalls, as expected (`d_{N,2} ≍ (log N)^{−1/2}`). Cross-checks
(`checks.py`): the `N = 10` value agrees with the exact infinite sum over the period `lcm(1..10) = 2520`
(digamma weights); the identity `ρ_{1/k} = {m/k}` on `(1/(m+1), 1/m]` was checked in exact rationals (0 mismatches).
Endpoint family (illustration, `M = 10⁷`, tail not bounded rigorously): `‖1 + H_ε‖_{1.1} = 0.360, 0.205, 0.110,
0.0573` for `ε = 0.4, 0.2, 0.1, 0.05`. The residue identity of step 4 was checked numerically at `p = 8/7`
(`J = 3.2800455` numeric vs `3.2800452` formula at `ε = 0.1`; `3.43308` vs `3.4330798` at `ε = 0.05`; `jcheck.py`).

## Novelty

Searched: arXiv abs/e-print pages; Semantic Scholar citation lists of Báez-Duarte 2003 (math/0202141, 101 citing
works), Bagchi 2006 (math/0607733), Báez-Duarte 2005 (math/0505453) and Dimitrov–Oliveira 2016; OpenAlex; arXiv
math.NT/CA/FA listings for 2026-10 and 2026-11. (Evidence API out of quota; arXiv search API throttled.)
Closest prior work:
* Dimitrov–Oliveira, arXiv 1608.07887 (2016): for Dirichlet L-functions and `p ∈ (1,2]`, `L ≠ 0` on `Re s > 1/p` ⟺
  a *Hilbert-space* (weighted L² on the line `Re s = 1/p`, twisted functions `κ`) natural-number criterion. With
  Lemma M this yields our open range `p < 1/θ` in a few lines, so that range is **routine**. It cannot give the
  endpoint (weighted L² on the line `1/p` controls only the Lorentz norm `L^{p,2} ⊋ L^p`).
* Gallardo-Gutiérrez–Seco, arXiv 2206.02654 (Canad. Math. J.): in the sequence spaces `ℓ^p_α` (for `α = −2` this
  is our `S_p`) they prove only *necessity* (approximation ⟹ zero-free half-plane), their Theorem A.
* Ghosh–Kremnizer–Noor–Santos, Adv. Math. 455 (2024), arXiv 2206.00434: necessity in Noor's Hardy-space model
  `H^p(𝔻)`, `1 < p ≤ 2`, with the remark "It is unclear whether the converse … holds" (a different space; not
  resolved here).
* Báez-Duarte (math/0505453) left the L^p strengthening open; Bagchi, Burnol, Nikolski, Balazard–Saias, de Roton,
  Vasyunin, Darses–Hillion, Alouges–Darses–Hillion and Manzur–Quintero (2611.03924) work in Hilbert spaces.

To our search, the sufficiency direction in `L^p(0,1)` at the endpoint, hence the exact equivalence of Theorem A
and the unconditional `L^{8/7}` statement, is new; the explicit Riesz–Möbius approximants with a power rate
(Theorem B) are new as stated but routine. Two parallel `reviewer` passes (Lemma H, Γ, M, A1–A5; Lemma PW, P,
A6–A7) found the core argument sound and two minor, fixed slips (Lemma Γ(b) needs `ε ≤ 1`; the `σ − η > 1` case
of Lemma P now has the explicit constant `1/(1−η) + 1/η`). A further independent review is still advised before
any DISCOVERIES entry.

## Limits

The theorem moves no RH node: for `p < 2` the criterion is equivalent to the zero-free half-plane `Re s > 1/p`,
so it is no stronger than the known 7/8 region at `p = 8/7` and no weaker than RH as `p → 2`. It gives a
Banach-space (non-Hilbert) form of Báez-Duarte's criterion in which RH reads "the criterion holds for every
`p < 2`", and its approximants converge at a power rate rather than a logarithmic one. `d_{N,2} → 0` is untouched.

RH relevance: N5
