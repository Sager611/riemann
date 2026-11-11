# explore25 derivations

## D1 Chandrasekhar-Darboux = Maass weight shift (T1)
- Schwarzschild: V_RW, V_Z = W^2 -/+ W' + beta (superpartners); intertwiners A_+/- = -/+ d/dr* + W.
  Reflection amplitudes: R_Z(w) = R_RW(w) * (mu^2(mu^2+2) - 12iwM)/(mu^2(mu^2+2) + 12iwM), mu^2 = (l-1)(l+2):
  unimodular rational factor for real w; same transmission and same QNMs except the algebraically special
  frequency w_as = -i mu^2(mu^2+2)/(12M) (Chandrasekhar 1975; Glampedakis-Johnson-Kennefick 1702.06459).
- Modular surface, weight k (even): E_k(z,s) = sum_{Gamma_inf \ Gamma} j_g(z)^{-k} Im(gz)^s, j = (cz+d)/|cz+d|.
  Constant term y^s + phi_k(s) y^{1-s}, phi_k(s) = I_k(s) zeta(2s-1)/zeta(2s),
  I_k(s) = int_R ((u+i)/|u+i|)^{-k} (1+u^2)^{-s} du = (-1)^{k/2} sqrt(pi) Gamma(s-1/2) Gamma(s) / (Gamma(s+k/2) Gamma(s-k/2)).
- Maass raising R_k = iy d_x + y d_y + k/2 maps E_k -> (s+k/2) E_{k+2}; on the constant term
  R_k y^s = (s+k/2) y^s, R_k y^{1-s} = (1-s+k/2) y^{1-s}, so
  phi_{k+2}(s) = phi_k(s) C_k(s),  C_k(s) = (1-s+k/2)/(s+k/2),  |C_k(1/2+it)| = 1.
  Zero of C_k at s = 1+k/2 (k >= 2: holomorphic Eisenstein series of weight k+2, a lowest-weight point = the
  "algebraically special mode"; k = 0: cancels the pole of phi_0 at s = 1, leaving E_2*), pole at s = -k/2 = 1 - (1+k/2).
  Same structure as Chandrasekhar's factor.
- The arithmetic factor zeta(2s-1)/zeta(2s) is untouched: poles rho/2 are common to all weights ("isospectral").
- Crum/Darboux chains multiply the S-matrix by finite Blaschke-type factors. Moved-pair model:
  xi_mod = xi * R, R rational, so phi_mod(s) = phi(s) R(2s-1)/R(2s): a finite Darboux chain maps zeta's scattering
  matrix to the moved-pair one and back. Hence no Darboux/Chandrasekhar invariant separates RH-true from RH-false.
- DH: the weight shift acts only on the Gamma factor, so a formal phi_DH(s) = Lambda_DH(2s-1)/Lambda_DH(2s) carries the
  same C_k and the same unitarity |phi(1/2+it)| = 1 (unitarity needs only FE + reality).

## D2 Structured pseudospectrum and EP births (T2)
- Class: entire F with F(s) = F(1-s), real on Re s = 1/2, gamma factor (13/pi)^{s/2} Gamma(s/2) (even, q = 13).
  Real basis on the line: Z_12 (quadratic chi, an L-function), Z_A = Z_3 + Z_9, Z_B = Z_4 + Z_10 (DH-type blocks;
  checked Z_9(t) = Z_3(-t), Z_10(t) = Z_4(-t), so Z_A, Z_B are FE-symmetric with root number +1).
- Family F_c = cos c L_12 + sin c G_A. Zeros on the line: r(t) := Z_12/Z_A = -tan c.
  Real zeros stay between consecutive zeros of Z_A while on the line (Z_c = cos c Z_12 != 0 there).
  At a critical point t* of r on the line, Im(L/G) = 0 also along a perpendicular curve; when -tan c passes a local
  MINIMUM value r* < 0 two real zeros collide (square-root EP) and leave the line along that curve; at a local MAXIMUM
  r* < 0 a pair re-enters. So the structured epsilon-pseudospectrum {s : L/G(s) in [-eps, 0]} is read off the critical
  values of one real function, and the first off-line height is T_1(c) = min{t* : negative local min of r, |r*| < tan c}
  for c > 0 (for c < 0: positive local maxima of r with r* < |tan c|). Near a birth, beta - 1/2 ~ sqrt(c - c*) (checked: ratio 6.7
  vs sqrt(47) = 6.9).
- Structured pseudospectrum := union over 0 <= tan c' <= eps of zero sets = {s : (L/G)(s) in [-eps, 0]}, the scalar
  analogue of {z : ||(A - z)^{-1}|| >= 1/eps} restricted to admissible (FE-preserving, real) perturbations (2611.00853).

## D3 Box / mirror (T3)
- Mirror (Dirichlet wall) at cusp height a: normal modes = zeros of B_a(s) = a^s X(2s) + a^{1-s} X(2s-1) (constant term of
  E(., s) times X(2s); Colin de Verdiere's pseudo-Laplacian; Lagarias-Suzuki math/0412039 for X = Lambda: all zeros on
  Re s = 1/2 for a >= 1 except two real ones iff a > 4 pi e^{-gamma}).
- Hermite-Biehler: E(s) = a^s xi_X(2s), E#(s) = conj E(1 - conj s); B_a = E + E#, and on the line B_a = 2 Re E.
  |E#/E| = a^{1-2 sigma} |xi_X(2s-1)/xi_X(2s)|. Group each zero rho = beta + i gamma with its FE partner 1 - conj(rho)
  = 1 - beta + i gamma; with w = 2s = 1 + x + i tau (x = 2 sigma - 1 > 0), d = tau - gamma:
  |w-1-rho| |w-1-rho'| / (|w-rho| |w-rho'|) = (|x-beta+id| / |x+beta+id|) (|x-(1-beta)+id| / |x+(1-beta)+id|),
  < 1 for every x > 0 when 0 < beta < 1; for beta > 1 the second ratio exceeds 1 (block A mod 13: HB fails only near
  its zero 1.183 + 77.19i, box mode 0.67836243 + 38.547714i off the line). So the box spectrum is real whenever every zero has 0 < Re rho < 1 (PNT-level), equivalently
  Re(xi'/xi)(1+iu) = sum_rho (1-beta)/|1+iu-rho|^2 > 0 (Lagarias positivity at sigma = 1, unconditional).
  Moved pair (any delta < 1/2) and DH below its first zero with beta >= 1 pass: the mirror cannot see RH.
