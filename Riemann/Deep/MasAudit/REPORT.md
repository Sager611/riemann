GAP — the decisive step, identifying the trace of the strip-restricted global geometric kernel with Weil's local terms (§12.3 "global to S-local deconstruction", p. 514: empty, footnote "Include section."; proof of Thm 12.5.1.1, p. 520: one line citing a definition), is unproved, implies GRH, and as written (independent of the Sobolev degree s and the strip width α) is false in an admissible case (K₁₂, α = 0, s ∈ {4, 6}); the GRH deduction itself (Ch. 13) is also unwritten.

# REPORT — masaudit: first audit of A. M. Masullo, *The Music of the Primes* (940 pp., 20 Nov 2026), claimed proof of GRH for Hecke L-functions of every number field

Inputs: `.scratch/masullo/music.pdf` (sha256 e039d3f6cf0a4b14…, 940 pp., pdfTeX, created 20 Nov 2026 18:13 UTC), `music.txt` (pdftotext). Every formula-level point below was read on the rendered PDF page (pp. 440, 514, 519–520, 522; others via text with page/line). Time used: 19:19–20:30 UTC.

## Short answer
1. The manuscript does not contain a proof of GRH. Its main theorem (Thm 12.5.1.1, p. 519) is Weil's explicit formula with the zero sum restricted to a strip |Re ρ − ½| < α; for α > ½ it is exactly Weil's formula in Connes' form (constants verified numerically, P1–P5), for small α it is GRH. All of the content sits in one step — equating the diagonal mass of the geometric kernel with the sum of its orbital pieces and identifying the remnant with Connes' S-local formula — and that step is an empty section (§12.3) plus one line on p. 520.
2. That step cannot be routine: as written it does not depend on s or α, yet with the manuscript's own spectral theorem it gives two incompatible identities for the A₄ field K₁₂ at α = 0, s = 4 and s = 6 (certified triple critical zero of ζ_{K₁₂}). An independent reviewer confirmed this (D1 CORRECT).
3. Chapter 13, which was to state GRH and deduce it, is two pages: no statement ("Insert statement."), footnote "State and prove the Main Theorems". The restated trace formula (Thm 13.3.1.1, footnote "Correct statement.") and §14.1.2 claim more than GRH: every non-central zero has multiplicity ≤ [F:Q] (for Q: all zeros of ζ and of Dirichlet L-functions simple, except at s = ½), derived "with the Positivity Criterion (Theorem ??)", contradicting §1.4.6 / §1.7.2 ("Weil's Positivity Criterion is not used").
4. Nothing checked contradicts known results; the known-results test shows where the proof must pay: the global-to-S-local passage, which Connes (1999) showed equivalent to RH.

## 1. Claims and scope (OUTLINE.md)
- Thm 12.5.1.1 (p. 519, l. 31490–31535): F any number field, all Hecke (Grössen-)characters (Ω(C_F)), φ ∈ D(C_F), φ(C_F^1) = 0; zeros in the interior of ½ + I with full multiplicity. Proof p. 520 (l. 31536–31563).
- Thm 13.3.1.1 (p. 522, l. 31595–31635): only zeros with σ(c) = ½, weights N_spec = min{N(c), [F:Q]} — the α = 0 case of the spectral theorem (cap (10.2.1.6), p. 440) at s = 2[F:Q] (p. 533: "is forced"); test functions in A: Mellin transform vanishing to order [F:Q] at the central point s = ½ − ix₀(χ) of every component ("[F:Q]-codimensional": in fact one block of [F:Q] conditions per component).
- Thm 13.3.1.2 "Divisor of L(s,χ)": empty. Per §14.1.2 (p. 528, l. 31915–31925): GRH + N(ρ) ≤ [F:Q] off the central point; no information at the central point (p. 529).
- Deduction sketch §1.4.2–1.4.5 (p. 28–30): approximate φ_z = t^{−z}1_{t>1}; L′/L(s+½) − δ(1/(s+½)+1/(s−½)) holomorphic on Re s > a for every a > 0.

## 2. Unfinished places (OUTLINE.md, logs/unfinished.log: 76 editorial footnotes, 14 '??')
On the main path: "Finish!" ×2 (§1.3.3–1.3.4, p. 27); 25 × "Return here." in Ch. 4–11 (p. 199–483); Thm 10.4.2.8 "Return here and complete proof." (p. 464); **§12.3 "Include section." (p. 514)**; §12.4: "Extend Definition 9.3.3.1 … Then study the convergence of the integral in (12.4.1.23)", "Return here after the Proposition is updated." (p. 514–515); §12.5 promises "the necessity of the choice a = ½ is explained" (p. 519): not there; Ch. 13: "Correct statement.", "Insert statement.", "State and prove the Main Theorems …" (p. 522); §14.1.2 "Correct this subsection.", "Theorem ??" (p. 528, 533). Appendix B (explicit formulas): "Brian should double-check this.", "Brian needs to recheck the proof, both for factors of 2 and 2π …" (p. 616, 619).

## 3. Crux and falsification (CRUX.md, PROBES.md)
- **Where the strip comes from** (p. 31, §1.4.8; Thm 10.2.1.1, 10.2.2.1, (10.3.2.3) p. 454): trace of convolution on W_I^{s,2}(C_F)/closure(|·|^{1/2}Θ); evaluations at off-strip zeros are discontinuous on W_I, so those zeros leave the spectrum. KNOWN in substance (Connes 1999 Thm 1, L²_δ, same multiplicity cap for even s; Meyer 2005). True for any θ-realizable Dirichlet series: proves nothing.
- **Geometric kernel = complement by definition** (Def. 10.4.3.7, p. 473): k^geom := φ(xy^{−1}) − k^spec, so ∫k^geom(x,x) = −(strip zero sum).
- **Written geometric pieces are local**: Thm 11.3.2.1 (p. 495): stratum trace = (log q_v/μ_v(O_v^×))Σ_k φ(π_v^k) (Poisson summation; 0 for q_v beyond supp φ); Thm 12.4.1.1 (p. 514) = Connes' S-local formula. No dependence on I, s or zeros.
- **Decisive step** (★) ∫k = ∫k_S + Σ_{p∉S}∫k_p and (★★) ∫_{C_F}k_S = local terms − pole terms (the S-local theorem has no pole terms): together they imply GRH; supported only by the definition (11.3.1.4) and an empty §12.3.
- **(a) Davenport–Heilbronn**: the spectral side and §1.4.8's topological argument apply verbatim to Davenport–Heilbronn-type series (true but empty); Thm 11.3.2.1 and 12.4.1.1 use only local Euler factors and hold for every Euler product. The only step that could separate L from a Davenport–Heilbronn function is the missing one. No written step is refuted by DH; the test locates the gap.
- **(b) Multiplicities**: 13.3.1.1 − Weil on A ⇒ GRH + N ≤ [F:Q] off the central point (the manuscript says so, p. 528). Tests: Artin-forced multiplicities d ≤ √[K:Q] (S₃ double zero, A₄ triple zeros: P6) obey the bound; CM central zeros of order ≥ 3 > 2 (E₁₂₅₄, rank 3, w = −1, GZK: P8, unconditional) are excluded by (13.3.1.2) (M̂ vanishes there to order [F:Q]) — no contradiction; F = Q consequence (simple zeros of ζ and Dirichlet L) consistent with data (P7) but open — a claim beyond GRH.
- **(c) Approximation φ_z = t^{−z}1_{t>1}**: spectral side not absolutely convergent (Σ 2/|z−iγ| = 1.18, 2.34, 4.06 for T = 100, 300, 1000; symmetric sums 0.0511, 0.0613, 0.0662: P11). §1.4.3's "a similar reasoning" with Montel lacks the needed local uniform bounds; fixable by classical zero-counting and density estimates; not the real issue (subtracting Weil from (12.5.1.1) gives GRH without approximation).
- **(d) Finite probes** (P1–P15): constants of (12.5.1.1) correct for 5 fields (diff ≤ 2.4e−8; wrong-constant control 0.27); (12.5.1.1) literally false at α = 0 (empty interior; Q: 0.009048 vs 1e−8 error) — the capped edge terms (10.2.2.5) are omitted; central-point definition (p. 528–529) literally false (Q(i): (1−5^{−s})^{−2} has double poles at all 2πik/log 5); Thm 1.7.1.1(1),(3) true (elementary); m_p(ψ) = 1 for the CM character of y² = x³ − x for all 4783 split p < 10⁵.
- **Sobolev test (decisive for the written proof)**: α = 0, F = K₁₂ (A₄ closure of x⁴+8x+12), χ = 1; Hardy Z of L(ρ₃) changes sign in (1.80, 1.83), so ζ_{K₁₂} = ζ_{k₃}L(ρ₃)³ has a critical zero of order ≥ 3; the spectral side counts it 2 times at s = 4, 3 times at s = 6; the asserted geometric side is identical; p∉S strata vanish for every cap. So (★)+(★★) as written are false for some admissible s; the manuscript's stated reason for s ≥ 2[F:Q] (Remark 11.2.2.2) is idle here. Hence §12.3 must carry s- and α-dependence that nothing written supplies. The consistency threshold is global (cap ≥ largest critical-zero multiplicity, s ≥ 5 on K₁₂), not the local s ≥ 2[F:Q]; at s = 2[F:Q] the a-line gluing holds iff all non-central critical zeros have N ≤ [F:Q], the bound it is used to prove.

## 4. Known results and barriers
| step | status |
|---|---|
| (12.5.1.1) for α > ½ | KNOWN: Weil's explicit formula in Connes' form (P1–P5) |
| spectral side, strip/line, multiplicity caps (Ch. 10) | KNOWN in substance: Connes 1999 (L²_δ, Thm 1), Meyer 2005 (Duke) |
| stratum traces (Thm 11.3.2.1) | KNOWN-type: per-place Poisson summation = prime-power terms |
| S-local trace (Thm 12.4.1.1) | KNOWN: Connes 1999 §VII Thm 4 (cited) |
| gluing (★)+(★★), "massive cancellation" (§1.4.10) | NEW, unproved; implies GRH; as written CONTRADICTS the manuscript's own spectral theorem at α = 0 (K₁₂, s = 4 vs 6) — this is the S-local → global passage that Connes showed equivalent to RH |
| GRH from (12.5.1.1) (§1.4.5) | KNOWN-type (trace formula ⇒ RH, Connes 1999); Ch. 13 empty |
| multiplicity bound N ≤ [F:Q] (§14.1.2) | NEW claim beyond GRH (for Q: simple zeros); consistent with P6–P8 data |
| "positivity not used" (§1.4.6) vs "with the Positivity Criterion (Theorem ??)" (§14.1.2) | internal CONTRADICTION about the route |
| Connes–Consani positivity barrier | not engaged: positivity-type φ = ψ∗ψ* violate φ(C_F^1) = 0 (p. 35), so the manuscript bypasses it, and the burden moves entirely to (★)+(★★) |
| Deninger | Ch. 8 "periodic orbits and fixed points" relabels the local Euler-factor poles (log q_v-periodic orbits); no new input |

## 5. Review
Two independent `reviewer` rounds (statements and arguments pasted inline). Round 1: C1 (GAP: the gluing and pole terms are unproved) CORRECT; C2/C3 wording WRONG (equivalence overstated; defect sum incomplete; zero not certified) → fixed (implication only; full defect with positivity; certified sign change P6c). Round 2: D1 (s = 4 vs 6 contradiction on K₁₂ refutes the parameter-uniform gluing) CORRECT; D3 (verdict GAP, D1 a defect of the written proof at α = 0, not a refutation of Thm 12.5.1.1 for α > 0 or of 13.3.1.1 at s = 2[F:Q]) CORRECT; D2 "equivalent" → "implies" (fixed in CRUX §4).

## 6. What would settle it
- The author's §12.3. Test any version on the a-line (α = 0) for K₁₂ at s = 4, 6, 24: a valid gluing must reproduce the caps 2, 3, ≥ 3; and on strips it must use zeros of L in |Re s − ½| > α.
- Checked (PROBES P15): §6.1, §7.5 and B.7 use 1/L only on Re s = 1 (L_S(χ, 1 − s_c), s_c a local pole) and shift contours only inside the classical zero-free region: no hidden GRH input there; the GRH content is not smuggled in elsewhere, it is absent (§12.3).
- Not checked: Ch. 11.2 (Thm 11.2.2.1), Ch. 2–9, Appendices (lctvs, stacks), the proofs of Thms 10.2.2.1, 11.3.2.1 and 12.4.1.1 line by line; Ch. 9 definition of W^{s,2} versus the odd-s cap in (10.2.1.6) (P12 side remark).

## Files
- `Riemann/Deep/MasAudit/REPORT.md` (this), `OUTLINE.md` (claims, skeleton, unfinished list), `CRUX.md` (≤ 1.5 pp.), `PROBES.md` (P1–P15)
- `code/p1_explicit.py` (P1–P5, P13), `code/p2_mult.py` (P6–P10), `code/p6b_a4.py`, `code/p6c_certify.py` (P6, P12), `code/p11_approx.py` (P11)
- `logs/p1_explicit.log`, `logs/p2_mult.log`, `logs/p6b_a4.log`, `logs/p6c_certify.log`, `logs/p11_approx.log`, `logs/unfinished.log`; memcap peaks ≤ 27 MB (one run killed at its 600 s limit, rerun to height 1000)

