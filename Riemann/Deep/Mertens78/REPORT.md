RECORD: |M(x)| ≤ 1.84·10^6 · x^0.9382 for all x ≥ 1 — the first explicit power-saving bound for the Mertens function M(x) = Σ_{n≤x} μ(n). It beats Ramaré's |M(x)| ≤ 0.013 x/log x from x ≈ 10^174 on. Certified family down to |M(x)| ≤ 1.01·10^382 · x^0.8992 (all x ≥ 1); also |M(x)| ≤ x^0.95 for x ≥ 10^395.

# Explicit power saving for the Mertens function from the zero-free half-plane Re s > 7/8 (node N10)

Agent mertens78, 2026-11-05 (16:56–18:52 UTC). Files in this folder:
- PROOF.md: full derivation (Lemmas A1–A4, B1–B6), certificate checks, two self-review passes, results.
- cert.py, run.py, batch.py, chk.py: the rigorous checker (Arb ball arithmetic via python-flint 0.9.0). all.json: the nine certified configurations. results.json, verify.log: their output (fresh re-run from this folder).
- c_*.json: the eight LP certificates (finite lists of decimal masses). lp.py, mkcert.py, opt.py: floating-point search that produced them (untrusted; only the checker matters). sanity.py: numerical sanity checks.

## 1. Short answer

1. **Theorem.** For every x ≥ 1, |M(x)| ≤ 1.84·10^6 x^0.9382. The simplest version, with a one-line certificate (Harnack's inequality, no linear programming): |M(x)| ≤ 7.91·10^5 x^0.9409.
2. **Only input about zeros:** ζ(s) ≠ 0 for Re s > 7/8 (OpenAI 2026, Lean-checked in Riemann/Deep/QRH). Other inputs: Rademacher's explicit convexity bound for |ζ| (Math. Z. 1959), the Herglotz representation of positive harmonic functions on a half-plane, |1/ζ(s)| ≤ ζ(σ)/ζ(2σ) for σ > 1, and Arb.
3. **Exponent versus constant.** Exponents from 0.9409 down to 0.8992 are certified; the constant grows from 10^6 to 10^382 (table §2). Harnack certificates alone cannot reach θ = 0.92 (the infimum, PROOF.md Remark after A4); the signed certificates found by linear programming go further.
4. **Practical range.** Our bound is smaller than Ramaré's 0.013 x/log x only for x ≥ 10^174 (headline row). For all smaller x, Ramaré's bound is stronger. The true size of M(x) is about √x, so the result is a first explicit power saving, not a sharp estimate.
5. **Novelty.** We found no explicit bound |M(x)| ≤ C x^θ with θ < 1 and explicit C, x₀ (arXiv math.NT listings for Oct and Nov 2026, arXiv keyword search, Zenodo, our own files; §5). The nearest result is Zenodo "Kernel-checked explicit Chebyshev bounds from a zero-free half-plane" (2026-10-19): |ψ(x) − x| ≤ 700 x^{61/64} for x ≥ 2^256. It covers ψ only. Bounding M needs bounds for 1/ζ, which are harder than bounds for ζ′/ζ.

## 2. Certified results

| certificate | θ | C, valid for x ≥ x_a | C′, valid for all x ≥ 1 | beats 0.013x/log x from |
|---|---|---|---|---|
| H093: Harnack, σ₁ = 0.93, σ₀ = 1.02, η = 0.1 | 0.9409 | 6.10·10^5 (x ≥ 10^163) | 7.91·10^5 | 10^175.8 |
| ML093g5, σ₁ = 0.93, η = 0.1 | **0.9382** | 1.50·10^6 (x ≥ 10^160) | **1.84·10^6** | **10^174.0** |
| ML092g8, σ₁ = 0.92, η = 0.1 | 0.9296 | 1.53·10^8 (x ≥ 10^169) | 1.85·10^8 | 10^181.4 |
| ML092g8, σ₁ = 0.92, η = 0.05 | 0.9277 | 6.22·10^8 (x ≥ 10^172) | 7.54·10^8 | 10^185.2 |
| ML091g12, σ₁ = 0.91, η = 0.03 | 0.9192 | 2.94·10^12 (x ≥ 10^200) | 3.54·10^12 | 10^211.9 |
| S090g20, σ₁ = 0.90, η = 0.01 | 0.9121 | 2.75·10^20 (x ≥ 10^275) | 3.35·10^20 | 10^287.0 |
| S090g40, σ₁ = 0.90, η = 0.005 | 0.9092 | 1.03·10^26 (x ≥ 10^326) | 1.23·10^26 | 10^340.0 |
| S089g400, σ₁ = 0.89, η = 0.005 | 0.9024 | 5.60·10^110 (x ≥ 10^1173) | 6.70·10^110 | 10^1190.1 |
| S089g3000, σ₁ = 0.89, η = 0.005 | 0.8992 | 8.20·10^381 (x ≥ 10^3826) | 1.01·10^382 | 10^3847.6 |

Pure powers |M(x)| ≤ x^θ′ (best row for each θ′):

| θ′ | 0.95 | 0.94 | 0.93 | 0.92 | 0.91 | 0.90 |
|---|---|---|---|---|---|---|
| holds for x ≥ | 10^395 | 10^600 | 10^1142 | 10^2409 | 10^14573 | 10^477393 |

Before: no explicit θ < 1. The explicit records were Ramaré's 0.013 x/log x (x ≥ 1,078,853), Cohen–Dress–El Marraki's x/4345 (x ≥ 2,160,535), and the sub-power savings of arXiv 2208.06141. The non-explicit M(x) = O(x^{7/8+ε}) is Lean-checked in Riemann/Deep/PNTError.

## 3. Method (details and proofs in PROOF.md)

1. **Harmonic function.** Z(s) = (s − 1)ζ(s)/s has no zeros or poles on Re s > 7/8, so u = log|Z| is harmonic there. Rademacher's bound gives an explicit harmonic majorant h(s) = A + log|(1+s)/s| + κ log|1+s| with κ = (1/8 + η)/2 (Lemma A1). So v = h − u ≥ 0 is harmonic.
2. **Lower bound for |ζ| (upper bound for 1/ζ).** For a positive harmonic function, the value at σ₁ + it is controlled by values on lines σ_j > 1, where |ζ| is trivially bounded above and below. In the simplest case this is Harnack's inequality v(σ₁+it) ≤ ((σ₀ − 7/8)/(σ₁ − 7/8)) v(σ₀+it). In general it is a finite signed "certificate" λ with Σλ_j P_{x_j}(w − y_j) ≥ P_{x₁}(w) for every real w, where P is the Poisson kernel (Lemma A2, via Herglotz). This gives |1/ζ(σ₁ + it)| ≤ C_∞ |t|^a for |t| ≥ 1000, with a = (Σλ_j − 1)κ (Theorem A3, Corollary A4).
3. **Smoothed Perron formula.** ∫₀^y M = (1/2πi)∫ y^{s+1}/(s(s+1)ζ(s)) ds is absolutely convergent and moves to Re s = σ₁. No residue appears: 1/ζ vanishes at s = 1. Averaging M over a window of length h costs h/4 + 1 (Lemma B4). For 0 ≤ t ≤ 1000 the true |1/ζ(σ₁+it)| is bounded with Arb; it is at most 5.7 there.
4. **Assembly.** With h = c_h x^θ and θ = (σ₁ + a)/(1 + a), every term of the bound divided by x^θ is nonincreasing in x. So one ball-arithmetic evaluation at x_a proves the bound for all x ≥ x_a, and |M(x)| ≤ x covers 1 ≤ x ≤ x_a (Theorem B6).

Where the losses are:

| step | loss | could be improved by |
|---|---|---|
| Rademacher's majorant | log(3ζ(1+η)) ≈ 3–6 in A; κ ≥ 1/16 | sharper explicit convexity on Re s = 7/8 (e.g. explicit approximate functional equation) |
| reference line σ₀ > 1 | Σλ⁺ log(ζ(σ₀)/ζ(2σ₀)) ≈ 10–12 | an explicit lower bound for |ζ(1+it)| (reference at σ = 1) |
| Harnack vs LP certificate | Harnack: Λ = x₀/x₁; LP: Λ ≈ 2 with larger Σ|λ| | larger LPs (more lines and points) push θ toward 7/8 at rising cost |
| unsmoothing | h/4 + 1 | squarefree density 6/π² (effect < 10% on C) |

| approach | θ reached | constant | proof weight |
|---|---|---|---|
| Harnack only | > 0.92 (best crossover row: 0.9409) | 10^6 | one inequality, no computation for the certificate |
| LP certificates (this work) | 0.9382 … 0.8992 | 10^6 … 10^382 | finite certificate, checked by Arb bisection |
| Borel–Carathéodory + three circles (Titchmarsh §14.2) | any θ > 7/8 in principle | far larger in practice | classical, sub-polynomial 1/ζ bound |

## 4. Certificates and how to reproduce

Run from this folder (single-threaded, about 12 min):
uv run --no-project --with python-flint python3 batch.py all.json   (writes results.json, prints one line per configuration)

What it checks in ball arithmetic (120 bits) for each configuration:
- (C1) and (C2) for the certificate. Harnack: analytic. LP: the tail bracket at |w| = Y and adaptive bisection on [−Y, Y]. S089g3000 needs about 10^6 balls.
- the constants of Corollary A4;
- J₀ (the Arb integral of |1/ζ| for 0 ≤ t ≤ 1000);
- C at x_a (Theorem B6).
LP certificates were regenerated with constraints out to |w| = 10^5 after an earlier version failed the far-tail check. Only certificates that pass the checker are reported.
Numerical sanity checks (sanity.py, not part of the proof): Rademacher's bound is at least 8.1·|ζ| on a test grid; the bound Φ(t) is at least e^{13.8}·|1/ζ(σ₁+it)| on sample t from 0 to 10^5.

## 5. Novelty pass (2026-11-05)

- arXiv math.NT monthly listings 2026-10 (854 entries) and 2026-11 (140) were scanned by title for Mertens, Möbius, reciprocal, explicit, zero-free, half-plane and quasi-Riemann. We read the abstracts of 2610.12234 and 2610.15291 (improved half-planes 0.874957 and 0.8748867), 2610.18491 (class-group applications of quasi-RH), 2610.23660, 2610.17720 and 2611.03356. None bounds M(x) explicitly.
- arXiv keyword search (four queries) returned only 2208.06141: explicit 1/ζ bounds in the classical, Littlewood and Korobov–Vinogradov regions, and M(x) bounds with sub-power savings.
- 2604.24505 (an effective Newman Tauberian theorem with an "explicit error term" for Mertens) was already checked in round R336: it gives no power saving.
- Zenodo (API, best match for "Mertens function", "reciprocal of the Riemann zeta" and "zero-free half-plane"): only the ψ-bound record above; no M(x) power saving.
- Our files (DISCOVERIES, COROLLARIES, discoveries/*): no Mertens power-saving entry.
Verdict: new as an explicit statement. Mathematically it is the explicit version of the classical deduction "zero-free half-plane implies M(x) = O(x^{θ+ε})". The certificate method is a linear-programming refinement of Harnack's inequality; we did not find it used for explicit 1/ζ bounds.

## 6. Caveats

- We cite two classical inputs and do not re-prove them. Rademacher's Theorem 4: we checked the transcription in the TeX of arXiv 2101.12263 and with a numerical grid. The Herglotz representation is textbook; it is needed only for signed certificates, and Harnack's case follows from the Poisson formula on discs.
- We used B = 7/8. Bebel's Lean-checked 0.873382… would lower θ by less than 0.001; not done.
- No Lean formalization (excluded by the brief).
- One reviewer pass, run by this agent on the draft of PROOF.md §0–2, returned FLAWED for missing hypotheses only. They were 0 < a < 1 in Lemma B5 and Theorem B6, and σ₁ ≥ 0.89 for the contour shift in B2. All three are now stated; every certificate has 0.09 ≤ a ≤ 0.19 and σ₁ ≥ 0.89. The reviewer found the signed-weight inequalities and the monotonicity argument correct. A formal reviewer verdict for DISCOVERIES is still the orchestrator's step.

RH relevance: N10.
