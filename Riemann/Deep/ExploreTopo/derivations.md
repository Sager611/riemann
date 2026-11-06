# explore18 derivations: topological phases of matter -> RH tree

Notation. F entire (or meromorphic in a strip) with F(s) = F(1-s) and F(conj s) = conj F(s); then A(t) := F(1/2+it) is real.
Zeros come in quadruples {rho, conj rho, 1-rho, 1-conj rho}; on-line zeros = real zeros of A. Test functions: xi, mv
(xi with its 6th/7th zeros 37.5862, 40.9187 replaced by 1/2 +- 0.15 + 39.2524i, FE and reality kept), DH = completed
Davenport-Heilbronn Lambda_f(s) = (5/pi)^{(s+1)/2} Gamma((s+1)/2) f(s).

## D1. Class D: the Z2 invariant is Hardy's sign change
- Physics: in class D (particle-hole symmetry C, C^2 = +1) the spectrum is symmetric E <-> -E; a single level at E = 0
  is protected, two can gap out. 0D invariant: sign Pf; 1D Kitaev invariant: sign Pf H(0) * sign Pf H(pi).
- Transfer: the FE + reality make the zero set symmetric under rho <-> 1 - conj rho (reflection in the line). Along any
  continuous FE-preserving family F_u, a simple on-line zero stays on the line until it meets another one: an off-line
  pair {rho, 1 - conj rho} has equal heights and is born only from a double real zero of A_u (Hurwitz). Hence the window
  invariant nu(I) = sign A(a) sign A(b) = (-1)^{N_0(I)} is the Kitaev invariant with the two window ends as the two
  TRIM points; N_0 >= #{windows with nu = -1}.
- This is Hardy 1914 (Gram-type windows) and Selberg 1942 (windows of length c/log T, parity detected through
  |int_I A| < int_I |A| with mollified moments). Pairs in one window are invisible, whether close on-line pairs or
  off-line pairs; DH has the same symmetry and the same invariant. On a finite window Z2 (edge) plus the argument
  principle N(T) (bulk) is Turing's method, i.e. how RH is verified to 3*10^12; it proves nothing asymptotically.

## D2. Class AIII/BDI: the integer winding is the argument principle (Levinson)
- Physics: chiral symmetry labels zero modes by chirality; only opposite chiralities annihilate; edge zero modes
  >= |winding of the bulk symbol| (Toeplitz index = -winding, Gohberg-Krein).
- Transfer: a companion is E holomorphic near the line with Re E(1/2+it) = w(t) A(t), w > 0. With E = |E| e^{i psi},
  A = 0 iff cos psi = 0, so (edge bound)  N_0(I) >= |W| - 1,  W := (1/pi) Delta_I psi.
  Bulk side: on the rectangle R = [1/2, X] x I, 2 pi kappa = winding of E around dR (kappa = zeros of E in R, the
  bulk defects), so W = N_Weyl - 2 kappa + b, b = horizontal-edge (eta-like) term, O(log T).
- Exact companion E_c = F + c F' (c > 0). On the line F'(1/2+it) = -i A'(t), so E_c = A - i c A' vanishes on the line
  only at a double zero of A. Hence kappa(c) cannot change with c (no zero crosses the line), and as c -> 0+ the zeros
  of E_c tend to those of F: an on-line zero rho goes to rho - c (left, since E_c ~ F'(rho)(s - rho + c)), an off-line
  zero stays. So kappa(c) = #{zeros of F with Re > 1/2} for every c > 0: the winding is exactly the argument principle,
  N_0 = N - 2 kappa. Polynomial check: A = z^2 + 1 in the z-plane gives E = z^2 + 1 - 2icz with zeros
  i(c +- sqrt(c^2+1)), one on each side for every c (kappa = 1 = number of non-real pairs).
- Limits: kappa(c) = 0 for all c is Lagarias' Re xi'/xi > 0 on Re s > 1/2 (Acta Arith. 1999; R20). c -> infinity gives
  zeros of F' and Laguerre-type defects (Speiser, Levinson-Montgomery; R9, R14, R47, R72). c = 2/L, L = log(t/2pi):
  xi + (2/L) xi' = h[(1 + (2/L) h'/h) zeta + (2/L) zeta'] with h'/h = (1/2) log(t/2 pi) + i pi/4 + O(1/t), so
  E_{2/L} = h (2 + O(1/L)) [zeta + zeta'/L]: Levinson's G = zeta + zeta'/L is the exact companion up to O(1/L).
- Where Levinson loses: kappa is bounded by Littlewood's lemma at sigma_0 = 1/2 - R/L,
  kappa <= (1/(1/2 - sigma_0)) sum_{beta > sigma_0} (beta - sigma_0); (beta - sigma_0)_+/(1/2 - sigma_0) is the
  least convex majorant of 1_{beta >= 1/2} supported on [sigma_0, inf) (Jensen/Littlewood weights are convex because
  log|G| is subharmonic), so the 'good' zeros of G near rho - 1/L (inside [sigma_0, 1/2]) are charged (R - 1)/R each,
  and the mollified mean value of |G M|^2 on sigma_0 does the rest. The topological step (edge >= |winding|) is
  lossless; all loss is analytic. Nothing in the topological dictionary improves the convex weight or the mean value.

## D3. Integer protection needs a definite Krein metric (= RH)
- Zeros that cannot annihilate pairwise need equal 'chirality' for all of them: psi' > 0 at every zero, i.e. E
  Hermite-Biehler with respect to the line (|E*(s)| < |E(s)| on Re s > 1/2, E*(s) = conj E(1 - conj s)). Then A has
  only real zeros: de Branges / Hermite-Biehler (R29, R34, R51), equivalent to RH. With kappa negative squares of the
  de Branges kernel, off-line pairs <= kappa (Krein-Langer; Pontryagin index, 1706.09103, R315, R325). For a real
  function alone the chirality sign A'(gamma) alternates, so only parity (D1) is protected: Xi on the line is class D,
  not class AIII.

## D4. The Euler product as a topologically trivial bulk (N7)
- Model F_B = G(s)B(s) + G(1-s)B(1-s), G a Gamma factor, B a Dirichlet series with B -> 1 as Re s -> inf. Companion
  E = G B: W = (theta_G/pi)|_I + (1/pi) Delta_I arg B(1/2+it), and Delta arg B along the line = -2 pi kappa_B +
  (horizontal edges of B), kappa_B = zeros of B in Re s > 1/2 (the bulk defects).
- B = P_X = prod_{p <= X} (1 - p^{-s})^{-1}: kappa_B = 0 and |arg P_X(sigma + it)| <= sum_{p<=X} arcsin p^{-sigma} on
  sigma >= 1/2, so every window has W = theta/pi + O(C_X): the line carries the full zeta-like count up to O(C_X)
  (compare Gonek 0704.3448: a related, non-analytic family built from finite Euler products satisfies RH with at most
  finitely many exceptions). The truncated product is not invertible on Re s = 0 (poles at 2 pi i k/log p); those
  'bulk gap closings' put one extra zero of F_B next to each pole for large t (distance ~ (t/2 pi)^{-1/2}), near
  Re s = 0 and 1, away from the line (test T3: none with 0.505 < Re s < 0.8 for X = 3, 29).
- B = Dirichlet polynomial (zeta or DH partial sums) or B = DH's combination of two Euler products: B has zeros in
  Re s > 1/2 with positive density (nonzero mean motion of arg B; Bohr-Jessen), each removing two units of winding.
- For zeta itself (X = infinity) the product diverges in the strip; 'B zero-free on Re s > 1/2 with bounded arg on the
  line' is RH plus S(t) = O(1), which is false (S(t) is unbounded). Circular, and the limit is exactly where it fails.
- Fredholm view (explore13 T3): Toeplitz operators with symbol zeta(sigma + it), 1/2 < sigma <= 1, are not Fredholm
  (zeta is unbounded above and below on such lines); the regularised Witten index (2501.15436) of non-Fredholm Toeplitz
  operators counts boundary zeros with weight 1/2, i.e. the argument principle with half residues.

## D5. Density invariants (Prodan-Schulz-Baldes, Bellissard gap labels) see only positive densities
- Index pairings of aperiodic media are traces per unit volume (IDS, Chern densities, Johnson-Moser rotation numbers,
  math-ph/0506016). The zero analogue on vertical lines is the Bohr-Jessen mean motion
  c(sigma) = lim (1/T) Delta arg F(sigma + it) = -2 pi (density of zeros right of sigma); on lines where an almost
  periodic F is bounded away from 0 it lies in the frequency module (Bohr), the gap-labelling statement.
- zeta: c(sigma) = 0 for all sigma > 1/2 (Bohr-Landau 1914, N(sigma, T) = o(T)); DH: c(sigma) != 0 on part of
  (1/2, 1.x) (Davenport-Heilbronn 1936; Voronin). So density invariants separate zeta from DH, but only by
  Bohr-Landau, and they are blind to o(T) sets of off-line zeros: they stop at zero-density estimates (N7 barrier).

## D6. Symmetry class of a Hilbert-Polya operator (N13)
- Spectrum {+-gamma}: H needs an operator anticommuting with it: chiral (unitary) or particle-hole (antiunitary).
  Berry-Keating xp with T = complex conjugation: T xp T^{-1} = -xp, T^2 = +1: class D. Bulk statistics are GUE
  (class A) because the pairing acts only near gamma = 0, where zeta has a gap (14.13). The 0D class-D invariant
  (sign Pf) is the parity of zero modes: the forced central zero for root number -1, i.e. Katz-Sarnak SO(odd)
  (AZ classes C, D, B for USp, SO(even), SO(odd)); for zeta (root number +1) trivial.
- No tenfold-way symmetry makes a spectrum real; reality comes from self-adjointness, or from an unbroken PT symmetry
  with a positive metric (Bender-Brody-Mueller 1608.03679 and Bellissard's comment 1704.02644), which is the definite
  Krein metric of D3. DH has the same symmetries (Lambda_f even and real on the line), so any symmetry-protection
  argument applies verbatim to DH, whose zeros do leave the line (in pairs, e.g. 0.808517 + 85.699348i).

## D7. Gap closing for Weil positivity (N1)
- Protection needs a uniform gap along the deformation. Weil's window form has lambda_min -> 0 as the window grows,
  even under RH (R84, about e^{-4 pi c}); the negative index is the Pontryagin index of D3 (= number of off-line pairs
  seen by the window). Nothing topological bounds lambda_min from below.

## D8. Isolated zeros are symmetry-protected zero modes (N8, already in use)
- A disc D symmetric about the line with exactly one zero of F contains an on-line zero (off-line zeros come in
  pairs rho, 1 - conj rho inside D); by Rouche this survives perturbations g with |g| < |F| on dD (gap = min |F| on
  dD). This is the mechanism of Bombieri-Hejhal (1995) and Lamzouri 2610.26841 for linear combinations, and of our
  isozero results (0.6724 c-isolated zeros, 0.3449 for combinations). Odd clusters (3 zeros in a symmetric disc) also
  force one on-line zero, pairs do not; the gain over singletons is the share of 3-clusters, negligible at small c.
