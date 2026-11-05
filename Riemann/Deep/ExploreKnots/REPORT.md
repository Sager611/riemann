NO IDEA: quantum knot invariants see L-functions only through special values, the functional equation of theta series, and Stokes constants equal to Dirichlet coefficients. The field also contains its own Davenport–Heilbronn counterexamples: every torus knot except the trefoil carries a sum of two or more Dirichlet L-functions. For T(2,5), the Lee–Yang combination vanishes at 1.19527 + 5.53931i, while every quantum structure stays intact.

# Exploration explore11: quantum knot invariants and quantum modularity → RH tree (2026-11-05)

Field: Kashaev invariants and colored Jones polynomials at roots of unity, strange identities, Zagier's quantum modular forms (QMF), the Habiro ring, state integrals, and resurgence. Brief: `.scratch/orca/explore11-brief.md`. Derivations: [derivations.md](derivations.md).

Notation:
- χ_12 = (12/·).
- ψ is the even primitive character mod 20 with ψ(3) = i, and ε = τ(ψ)/√20.
- φ is the golden ratio.
- T_k are Glaisher's numbers 1, 23, 1681, 257543, ….
- DH means a Davenport–Heilbronn-type function: a functional equation (FE) but no Euler product.

## Survey (one page)
- **Quantum modular forms.**
  - Zagier (Clay Proc. 12, 2010) gives the examples: Dedekind sums, the Kontsevich–Zagier F(q) = Σ(q;q)_n, Eichler integrals of unary theta series, Kashaev invariants, and Maass period functions.
  - Nonzero weight, limit laws and cotangent sums: Bettin–Drappeau 2210.07854. Cuspidal orbits on Fuchsian groups, Gaussian limits: Bettin–Drappeau–Lee 2610.41973 (Oct 2026).
  - All of Zagier's examples recast as resurgence: McSpirit–Rolen 2505.00799.
- **Strange identities.**
  - Zagier (Topology 2001): F(q) "=" −½Σ nχ_12(n) q^{(n²−1)/24}. The asymptotic coefficients of F at q = e^{−t} are L(−2k−1, χ_12) = −2(−1)^k T_k.
  - Hikami: the same holds for T(2,2m+1), via half-derivatives of Andrews–Gordon (math/0303250).
  - Hikami–Kirillov 2003: torus knots ↔ Virasoro minimal models M(s,t). Gang–Park–Sohn 2512.23122 (Dec 2025) re-derive this through the 3D–3D correspondence.
  - Hikami–Lovejoy 1409.6243: cyclotomic expansions for T(2,2t+1), and a family of QMFs dual to the generalized Kontsevich–Zagier series.
  - Resurgence of Habiro elements via strange identities: Crew–Fantini–Goswami–Osburn–Wheeler 2304.07001. Asymptotics and positivity of generalized Fishburn numbers: Goswami–Jha–Kim–Osburn 2204.02628.
  - Lawrence–Zagier 1999: the WRT invariant of the Poincaré sphere is a radial limit of a false theta series. Witten's asymptotic expansion conjecture for Seifert manifolds was proved in 2510.10678, with the framework in 2508.21710. Strong quantum modularity for Brieskorn spheres: 2512.10768.
- **Hyperbolic knots.**
  - Garoufalidis–Zagier 2111.06645 build a matrix invariant Φ_α(h) with entries in Habiro rings and a refined quantum modularity (an SL₂(ℤ) cocycle).
  - Garoufalidis–Yoon 2211.00270: asymptotic multiplicativity under cyclic covers.
  - The 4₁ invariant under modular transformations: 2510.02188.
  - The growth rate is the volume, a Borel regulator, e.g. vol(4₁) = (3√3/2)L(2, χ_{−3}).
  - Andersen–Kashaev state integrals use Faddeev's Φ_b, which has modulus 1 on ℝ. Teichmüller TQFT: 1109.6295, 2307.12848.
- **Habiro ring.**
  - Garoufalidis–Scholze–Wheeler–Zagier 2412.04241: the Habiro ring of a number field and modules graded by K₃, with Borel, p-adic and étale regulators. Elements come from the infinite Pochhammer symbol.
  - Follow-ups: cyclosyntomic (q-deformed p-adic) regulator 2602.21894 (Feb 2026); explicit Habiro-cohomology classes 2505.19885; q-Hodge complexes 2410.23078 and 2510.04782; a Habiro lift of the colored Jones polynomial 2603.01619 (Mar 2026).
- **Resurgence.**
  - Garoufalidis–Gu–Mariño 2007.10190 and 2012.00062: integer Stokes constants (peacock patterns, 3D index).
  - Fantini–Rella 2404.11550: "modular resurgence". In the single-tower case the Stokes constants are the coefficients of an L-function. Follow-ups: q-Pochhammer symbols and the topological-string/spectral-theory (TS/ST) operators of 2506.08265.
  - Cheng–Coman–Fantini–Rella 2608.08902 (Aug 2026): a vector-valued version, with Stokes constants built from a vector of L-functions; the unary theta / false theta case is worked out in full.
  - Broadhurst–Dorigoni 2507.21352: Lambert series with characters.
- **Central values.**
  - The Vasyunin cotangent sum in the Nyman–Beurling Gram matrix is a weight-1 QMF (Bettin–Conrey 2013).
  - Additive twists L(½, f ⊗ e(x)) form QMFs that give reciprocity for twisted moments: Nordentoft 1909.09665; Maass case 2208.14346.
- **DH-type with FE.** Nakamura 2008.02570: q^s L(s,χ) + i^{−κ}G(χ)L(s,χ̄) has Riemann's FE, and for non-real χ it has infinitely many zeros off the line.
- **Searches.**
  - arXiv API search for abstracts containing "Davenport-Heilbronn" together with knot, quantum, theta or "minimal model": no results.
  - No 2025–26 item links quantum invariants to the location of zeros.
  - This repo: "quantum modular forms / Habiro" appeared only in R328 (no idea). Lax–Phillips, Mayer, scattering and resurgence are already NOGO (R5.5, R48, R58, explore3, R339).

## (i) Where L-functions enter, and what could constrain zeros
1. **Special values.**
   - The strange-identity coefficients are L(−2k−1, χ), or L(−2k, χ) for odd χ.
   - Volumes are L(2, ·) and Borel regulators.
   - The Habiro ring encodes K₃ regulators: Borel, and p-adic in the cyclosyntomic version.
   - All of these are values. A zero enters only through the Hadamard product, at relative size O(δ²/k²) (T2).
2. **Functional equation.** Quantum modularity of these invariants is the modularity of a vector of unary theta series. That is Hecke's FE for the vector (L(s, χ_j))_j, with the S-matrix of the minimal model as the FE matrix.
3. **Stokes constants** are the Dirichlet coefficients χ(m) (D1; Fantini–Rella).
4. **Central values** of additive twists are QMFs. This gives exact reciprocity formulas for twisted moments (N9 tools) and the NB Gram matrix (N5); none of it says where zeros are.

The candidate constraints of the brief, and their counterparts:
- **Positivity.** ⟨4₁⟩_N is a sum of squares, and the Fishburn-type coefficients are positive. Both also hold for DH-type knots (T1).
- **Unitarity.** The modular S-matrix is real orthogonal. This is exactly what produces DH functions (D3).
- **Integrality.** The coefficients are integers for DH-type knots too (T1, T2).
- **Self-adjointness.** TS/ST operators are positive and trace class, but their Weyl law is ~E², not (T/2π)log T. State integrals are matrix elements of unitary operators (|Φ_b| = 1 on ℝ), with no spectrum tied to L-zeros. On the modular surface, ζ zeros are only scattering poles (explore3).
- None of these counterparts constrains zeros.

## Transfers, tests, verdicts

### T1 (N1, N12, N13): quantum-topological structure as an RH input
**Hypothesis.** Four properties of the trefoil's Kashaev invariant force GRH for its L-series L(s, χ_12):
- a strange identity with integral L-value coefficients,
- membership in the Habiro ring,
- quantum modularity with a unitary S-matrix,
- positive coefficients at q = 1 − x (the Fishburn numbers).

**Test (t1, t3c, t5, t6).** All four properties hold for T(2,5):
- Hikami's series Σ_{k₂≥k₁}(q;q)_{k₂} q^{k₁²+k₁}[k₂ choose k₁]_q equals −½Σ n Imψ(n) q^{(n²−9)/40} exactly to t¹³.
- Its coefficients L(−2k−1, Imψ) = −2, 142, −29282, 12485422, … are integers.
- Its S-matrix is the Lee–Yang matrix [[c, s], [s, −c]] with c + is = ε and tan arg ε = 1/φ (fit residual 4·10⁻¹⁶).
- X(1−x) = 1, 2, 6, 23, 109, 621, 4149, … is positive up to x³⁰.

**Zeros (t3, t3b, t3d, t3e; all found by Newton at 30 digits and re-checked at 40 digits with mpmath.dirichlet).**

| L-series | zeros right of the line, t < 60 (argument principle) | first zero |
|---|---|---|
| trefoil L(s, χ_12) (control) | 0 (grid search) | — |
| T(2,5) own series L(s, Im ψ) | 17 | 0.838544451285 + 8.292146000549i |
| T(2,5) Lee–Yang eigencombination H₋ (Λ(1−s) = −Λ(s)) | 10 | 1.195270648642 + 5.539314643116i |
| T(2,5) Lee–Yang eigencombination G₊ (Λ(1−s) = Λ(s)) | 0 below t = 60 (grid); 1 below t = 123.8 (argument principle) | 0.652400414407 + 122.846724764454i |
| T(3,4) = ½[L(s,(24/·)) − L(s,(12/·))] | 16 | 0.986306813311 + 6.982295921607i |
| T(2,7) (3 characters mod 28) | 12 | 0.820173484809 + 9.766689580103i |
| Poincaré sphere Σ(2,3,5), H₋ | 4 (grid search; argument principle not run) | 0.604987666972 + 19.746196916245i |

- Among T(s,t) with s < 8, t < 13, only the trefoil carries a single character; D2 proves this for all torus knots and all Brieskorn spheres.
- By Saias–Weingartner, every other torus knot's series has zeros in Re s > 1 (D2).

**DH test.** Fails inside the field. The field realizes DH functions: torus knots and Brieskorn spheres other than the trefoil give Davenport–Heilbronn-type combinations (the same type as Nakamura 2008.02570), with every listed property intact. Saias–Weingartner gives each of them infinitely many zeros in Re s > 1.

**Known-equivalent test.** No criterion is produced. The only distinguishing property is multiplicativity of χ (a Hecke/U_{p²} eigen-property), i.e. the Euler product, which is the existing DH filter.

**Verdict: dead.** This is a barrier for N1, N12 and N13: unitarity, integrality and positivity of quantum invariants cannot imply RH.

### T2 (N9, N3, N5): strange-identity coefficients and Stokes constants as zero detectors
**Test (t2).**
- T_k = 1, 23, 1681, 257543, 67637281, ….
- L(−2k−1, χ_12) ↔ L(2k+2, χ_12) via the FE, exact for k = 0, 1, 5, 20.
- Peeling T_0…T_139 returns the Stokes constants χ_12(1…13) exactly.
- What a zero would change: move the lowest zero ½ + 3.80462763305087i off the line by δ, as an FE-symmetric quadruple. The relative change of T_k is +0.087δ² (k = 0), −0.011δ² (k = 5) and −0.0011δ² (k = 20). T_0 = 1 would become 1.00087 at δ = 0.1.

**DH test.** Fails. T(2,5)'s coefficients are integers too, its Stokes constants are Im ψ(m), and it has zeros off the line. Integrality pins the function (R37), not where its zeros lie.

**Known-equivalent test.** The only zero-sensitive functional of the T_k is the Riesz sum Σ(−x)^{k+1}/(k! L(2k+2, χ_12)), with L(2k+2, χ_12) = 2πT_k/((12/π)^{2k+3/2} k! Γ(k+3/2)). It is O(x^{1/4+ε}) ⇔ GRH(χ_12). This is Riesz 1916 (R13.4), with zeros at height γ suppressed by e^{−πγ/4}. Li/Keiper is the derivative version (N3).

**Verdict: dead** (known equivalent).

### T3 (N13, N7, N9): hyperbolic knots, state integrals, Habiro ring of the trace field, twisted central values
**Test (t4).**
- ⟨4₁⟩_N = Σ|(q;q)_k|² > 0.
- log⟨4₁⟩_N − [NV/2π + 1.5 log N − ¼ log 3] = 0.5544/N at N = 3200. This matches Garoufalidis–Zagier's 11/(72√−3) term (0.5542).
- V = (3√3/2)L(2, χ_{−3}) = 2.029883212819307.

**DH test and known-equivalent test.**
- L-values appear only at s = 2 and as regulators (complex or p-adic), which DH-type functions also have.
- Twisted central values are QMFs for any modular form, including non-eigenform combinations (DH again).
- What is left is known: NB Gram matrix (N5), moment reciprocity (N9), Hilbert–Pólya barriers (R50, R51, explore3).

**Verdict: dead.**

## Ranking (RH relevance × chance × gain)
| transfer | node | chance | gain if it worked | verdict |
|---|---|---|---|---|
| T1 quantum structure (integrality, unitary S, positivity) ⇒ GRH | N1, N12, N13 | none: T(2,5), T(3,4), T(2,7) and Σ(2,3,5) have the structure and off-line zeros | GRH for χ_12 | dead; barrier (this report's main output) |
| T2 strange-identity coefficients / Stokes constants detect zeros | N9, N3, N5 | none: relative effect O(δ²/k²); the only criterion is Riesz | a new equivalent | dead; known (R13.4) |
| T3 volumes, Habiro ring, state integrals, twisted central values | N13, N7, N9 | none: values at s = 2, regulators, moment identities | spectral or positivity input | dead; known |

**New relative to R328** (which said "no idea"):
- the exact location of L-functions in the field;
- the trefoil is the only torus knot with an Euler-product L-series, and no Brieskorn sphere has one;
- explicit off-line zeros of knot L-series and of the Lee–Yang eigencombinations;
- the principle that unitarity of the S-matrix builds DH functions.

## Side result (off-tree, routine; novelty not checked beyond arXiv abstracts; folklore-level, not a discovery)
- Hikami's strange-identity L-series of a torus knot T(s,t) is a Dirichlet L-function only for the trefoil. Every other torus knot and every Brieskorn sphere gives a sum of two or more characters, with zeros in Re s > 1 (Saias–Weingartner).
- For T(2,5), the eigencombinations of the Lee–Yang S-matrix are Davenport–Heilbronn-type FE functions with zeros at 1.195270648642 + 5.539314643116i and 0.652400414407 + 122.846724764454i.

## Suggested NOGO line
- explore11 (Orca: quantum knot invariants, quantum modularity, Habiro ring) NO IDEA. L-functions enter only through special values (strange identities give L(−2k−1, χ); volumes give L(2, χ_{−3}); the Habiro ring gives Borel and p-adic regulators), the FE of unary theta vectors, and Stokes constants = Dirichlet coefficients. The field realizes DH itself: torus knots other than the trefoil, and Brieskorn spheres, give sums of ≥ 2 L-functions. T(2,5)'s Lee–Yang S-eigencombinations vanish at 1.19527+5.53931i and 0.65240+122.84672i, while its strange identity, integrality, unitary S and positive Fishburn-type coefficients all hold. [ExploreKnots](Riemann/Deep/ExploreKnots/REPORT.md)
- Possible FILTER line: a unitary (real orthogonal) modular S-matrix is what builds DH functions (±1 FE eigencombinations), so unitarity of quantum representations cannot imply RH.

## Files
`Riemann/Deep/ExploreKnots/`:
- `t1_strange.py` → `t1_strange_out.txt`
- `t2_sensitivity.py` → `t2_out.txt`
- `t3_zeros.py` → `t3_out.txt`, `t3_zeros_out.json`
- `t3b_check.py` → `t3b_out.txt`
- `t3c_smatrix.py` → `t3c_out.txt`
- `t3d_poincare.py` → `t3d_out.txt`, `t3d_out.json`
- `t3e_verify.py` → `t3e_out.txt`
- `t4_kashaev41.py` → `t4_out.txt`
- `t5_family.py` → `t5_out.txt`
- `t6_fishburn.py` → `t6_out.txt`
- `derivations.md` (D1–D5)

Each script ran alone and single-threaded under `perl -e 'alarm N'` (N ≤ 1700), with `uv run --no-project --with mpmath --with numpy` (plain python3 for t1 and t6). Longest run: 73 s (t3b). Peak memory: 73 MB. Scratch: `.scratch/explore11/STATE.md`. No downloads.

RH relevance: N1 and N13 (positivity, unitarity and self-adjoint structure from quantum topology: barrier, since the field contains DH functions); N12 (Habiro ring: special values and regulators only); N9 and N5 (twisted moments and the NB Gram matrix: QMF reciprocity is a known tool); N7 not touched.
