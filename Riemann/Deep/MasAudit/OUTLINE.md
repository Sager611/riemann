# OUTLINE — Masullo, *The Music of the Primes* (940 pp., PDF 20 Nov 2026 18:13 UTC)

PDF page = printed page (offset 0). `l.` = line in `.scratch/masullo/music.txt`. Formula claims confirmed on page images (see PROBES.md / logs).

## Structure (TOC, p. 3–5, l. 1–120)
| Ch. | pp. | content |
|---|---|---|
| 1 | 11–56 | Introduction (unfinished by author's note, p. 6 l. ~150); §1.4 trace formula and consequences (27–32); §1.5 Connes/Deninger (32); §1.6–1.8 method, letters to Sarnak (33–56) |
| 2–3 | 57–154 | adelic θ-functions; idelic zeta functions, Weil's identity |
| 4 | 155–220 | periodic Schauder bases (§4.5 without, §4.6 with multiplicity) |
| 5–6 | 221–310 | spaces of θ/Z-functions, Structure Theorems I–V, converse theorems |
| 7 | 311–380 | Schauder bases of idelic zeta functions: polar multiplicities, Paley–Wiener levels, nuclearity, strong co-primality |
| 8 | 381–394 | periodic orbits, fixed points, orbital traces |
| 9 | 395–432 | Hardy–Sobolev topologies, approximate units, (Λ,ε)-cutoffs |
| 10 | 433–474 | trace formula: spectral side (§10.4, 455–474) |
| 11 | 475–496 | shell structure; orbital kernels and traces (§11.3, 491–496) |
| 12 | 497–520 | geometric side; S-local trace (§12.4, 514–518); **Thm 12.5.1.1 (Trace Formula), p. 519–520** |
| 13 | 521–522 | **Zeros of Hecke L-functions: 2 pages; §13.1, §13.2 empty; §13.3 = restated Trace Formula + empty Thm 13.3.1.2** |
| 14 | 523–536 | multiplicities (§14.1, 524–533), function fields, higher rank |
| A–E | 537–932 | ANT, analytic NT (explicit formulas B.4), lctvs, Schwartz–Bruhat, topological stacks |

## Main claims as written
- **Thm 12.5.1.1** (p. 519, l. 31490–31535): F number field, φ ∈ D(C_F) with φ(C_F^1) = 0 [check on image]. Σ_{v real} ½∫′ φ|g|^{1/2}/|1−g| + Σ_{v cx} (1/2π)∫′ … + Σ_{v finite} (log q_v/μ_v(O_v^×))∫′ … = M(φ)(|·|^{1/2}) + M(φ)(|·|^{−1/2}) − Σ_{c,I} N(c) M(φ)(c|·|^{−1/2}); c over zeros of L in the interior of Ω(C_F)_I = {σ(c) ∈ ½ + I}; N(c) = full multiplicity. Proof (p. 520, l. 31536–31563): kernel k = k_S + Σ_{p∤S} k_p (11.3.1.4); Thm 11.3.2.1: ∫k_p(x,x) = 0; Thm 12.4.1.1: ∫k_S(x,x) = local terms − poles; (12.5.1.4) "for I = [−α, α]": ∫k(x,x) = −Σ_{c,I} N(c)M(φ)(c|·|^{−1/2}) "by design of the global geometric kernel (Def. 10.4.3.7) and by Thm 10.4.2.8".
- **Thm 13.3.1.1** (p. 522, l. 31595–31635), footnote 1 "Correct statement.": spectral side Σ over zeros with σ(c) = ½ **only**, weight N_spec(c) = min{N(c), [F:Q]} (13.3.1.1); test functions in "the [F:Q]-codimensional subspace A" cut out by (13.3.1.2): for 0 ≤ k ≤ [F:Q]−1 and all χ, ∂_s^k ∫_{R>0}∫_{C_F^1} χ(τ(t)) t^s φ(τ(t)g) χ(g) dμ¹ dt/t = 0 at s = −i x_0(χ) (χ(τ(t)) = t^{i x_0(χ)}).
- **Thm 13.3.1.2** "Divisor of L(s, χ)" (p. 522, l. 31636–31643): body empty; "Upon combining Theorem 13.3.1.1 and Theorem ??"; footnote 2 "Insert statement."; footnote 3 "State and prove the Main Theorems, also explaining how to deduce them by approximation φ → e^{−s|log(t)|} as Re(s) > 1."
- ⇒ **No statement of GRH and no proof of the deduction exists in the body.** Only the §1.4.5 sketch (p. 27–31).

- Thm 12.5.1.1 vs 13.3.1.1: same geometric side; 12.5.1.1 = strip version (I = [−α, α], interior zeros, full N(c); omits the capped boundary-line terms of (10.2.2.5); false as written for α = 0, where the interior is empty: PROBES P1, local − pole terms = 0.00905 ≠ 0 for Q); 13.3.1.1 = a-line version (I = {0}), i.e. Thm 10.2.1.1's cap N_s = min{N, ⌊(s+1)/2⌋} (p. 440, (10.2.1.6), image checked) at s = 2[F:Q] (p. 533 l. ~32190: "the choice of Sobolev degree s = 2[F:Q] … is forced").
- §14.1.2 (p. 528, l. 31915–31925): Thm 13.3.1.2 (not in the text) "a consequence of the Trace Formula combined with the Positivity Criterion (in its form from Theorem ??)" asserts GRH **and N(ρ) ≤ [F:Q] for every non-central zero**; for F = Q this is simplicity of all zeros of ζ and of all Dirichlet L-functions except s = ½. Contradicts §1.4.6 / §1.7.2 (p. 30, 35): "Weil's Positivity Criterion is not used". Fn 1 p. 528: "Correct this subsection."
- Central point (p. 528–529): c₀|·|^{1/2}, s = ½ − i x₀ on each component; multiplicity there "beyond the reach of the methods"; Prop. 14.1.2.5 (CM elliptic curves, conditional on (14.1.2.4)/BSD) shows it can exceed [F:Q].

## Logical skeleton (as written)
1. Spaces: θ-functions H_{0,cont}(F^×, S₀(A_F)) (Ch. 2, 5); every θ-function = Φ(f_S ⊗ ⊗_{v∉S} 1_{O_v}) (Structure Thm V, §6.2; uses Pic(O_{F,S}) = 1); Mellin transform = Z_S(f_S, c)·L^S(c) (§1.4.7, p. 30).
2. Hardy–Sobolev spaces W_I^{s,2}(C_F) with weight w_{I,s}(g) = max{|g|^{2a+2b}…}(1+|log|g||²)^{s/2} (Ch. 9; (9.2.2.41), quoted p. 473).
3. Spectral side (Ch. 10): quotient H^s_{F,a,I} = W_I^s / closure(|·|^a θ) has Markushevich basis indexed by zeros of L with real part in a + I (Thm 10.2.1.1 p. 440 for I = {0}, cap N_s; Thm 10.2.2.1 p. 445–446, full N inside, cap (10.2.2.5) on b = ±α); trace of A(φ)∘P^spec = Σ_strip N M(φ) (Thm 10.3.2.1, Thm 10.4.2.8 p. 463, proof fn 8 "Return here and complete proof.").
4. Kernels (Def. 10.4.3.7, p. 473, s ≥ 3): k^spec = k·w_{I,s}^{-1}; **k^geom := φ(xy^{-1}) − k^spec**. Hence ∫k^geom(x,x) = −Σ_strip N M(φ)(c|·|^{−1/2}) "by design" (12.5.1.4).
5. Orbital pieces (Ch. 11): divisorial strata α_k of GL₁(A_F); k_{p_k} := (α_{k+1} − α_k)k^geom (11.3.1.4, p. 494, a definition); §11.3.1 is an overview "omitting several items" (l. 29897–29907). Thm 11.3.2.1 (p. 495; a = ½, s > 1, α ≥ 0, φ supported on C_F^{>1}): orbital trace of a stratum = (log q_v/μ_v(O_v^×)) Σ_{k∈Z} φ(π_v^k): the local prime-power term of Weil's formula; zero for v with q_v > T by support. Uses Poisson summation; no input on zeros, I or s.
6. S-local trace (§12.4, Thm 12.4.1.1, p. 514): Connes' S-local formula [Co §VII Thm 4] for k^PW_{S,φ,λ,λ′} on J_{F,S}: local terms + φ(1) log(λλ′) + o(1); **no pole terms**. Proof fns: "Extend Definition 9.3.3.1 … Then study the convergence of the integral in (12.4.1.23)" (p. 514), "Return here after the Proposition is updated." (p. 515).
7. **Gluing (§12.3 "global to S-local deconstruction of the geometric kernel", p. 514): empty, fn 1 "Include section."** In the proof of 12.5.1.1 (p. 520) the gluing is the one line "By (11.3.1.4), this reads: ∫k_S + Σ_{p∤S}∫k_p = ∫k", and k_S (global, on C_F) is identified with the S-local k^PW of Thm 12.4.1.1, and the pole terms M(φ)(|·|^{±1/2}) appear, without argument.
8. Deduction of GRH (Ch. 13): absent. Sketch only in §1.4.2–1.4.5 (p. 28–30): φ_z(t) = t^{−z}1_{t>1}, mollify (ε) and truncate (h), Montel, finite-place sum → L′_f/L_f(z+½) on Re z > ½; conclude L′/L(s+½) holomorphic on Re s > a for every a > 0.
- Global analytic inputs used anywhere: non-vanishing and polynomial growth of 1/L_f on 1 + iR (Thm B.6.1.6, Strong Co-primality Thm 7.5.0.1 p. 374), Titchmarsh-type ordinates T_k (B.5), convexity (B.1), Weil's explicit formula (B.4, fns "Brian should double-check this", p. 616, 619).

## Unfinished places (full list: logs/unfinished.log, 76 editorial footnotes + 14 '??')
- On the main path: §1.3.3–1.3.4 "Finish!" ×2 (p. 27); 25 × "Return here." in Ch. 4–11 (p. 199–483, incl. Ch. 9 p. 414–426 and Ch. 10 p. 435–462); Thm 10.4.2.8 "Return here and complete proof." (p. 464); **§12.3 "Include section." (p. 514)**; §12.4 fns 2–4 (p. 514–515); Ch. 13 fns 1–3 (p. 522: "Correct statement.", "Insert statement.", "State and prove the Main Theorems …"); §14.1.2 "Correct this subsection." (p. 528); "Theorem ??" for the Positivity Criterion (p. 528, 533).
- Appendices: B.4 explicit formulas "Brian still needs to check proof …", "Brian should double-check this.", "Brian needs to recheck the proof, both for factors of 2 and 2π …" (p. 616, 619); B.2 "must either rewrite for accuracy or remove" (p. 603); B.5 "Brian needs to check the next two lemmas" (p. 633); A.3–A.5 "Double-check for the S-local case." (p. 556), "Why? Must say more." (p. 580); C, D, E: further "Must expand proof", "Check this definition is correct", "This is not a sentence" (p. 753, 862, 905).


## Further statement-level defects (image-checked)
- p. 519 (l. ~31470): "This is also where the necessity of the choice a = 1/2 is explained." §12.5 contains only the theorem and a 12-line proof: no such explanation.
- p. 522 vs p. 519: (13.3.1.4)–(13.3.1.5) drop the factors ½ (real), 1/2π (complex) and 1/μ_v(O_v^×) (finite) of (12.5.1.1)–(12.5.1.2) but keep log q_v: mixed normalisations (p. 514 calls ½μ_v, μ_v/2π, (log q_v/μ_v(O_v^×))μ_v the canonical measures, so either the ½ is missing or log q_v is doubled). PROBES P1: with dx/|x| the ½ is needed (factor 1 breaks Q by 0.27). Covered by the author's "Correct statement."
- p. 528: "ρ ∉ R, i.e. ρ ≠ 1/2" while the central point is s = ½ − ix₀(χ) (p. 529), not real when x₀ ≠ 0.
- Thm 12.5.1.1 omits the capped edge-line zeros of (10.3.2.3)/(10.2.2.5) and is false as written at α = 0 (PROBES P13).

