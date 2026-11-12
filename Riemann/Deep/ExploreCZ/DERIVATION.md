# explore46 derivations (append-only, <= 60 lines per append)

## 1. Setting (node N7, large values)
D(t) = sum_{N<n<=2N} a_n n^{-it}, |a_n| <= 1, G = sum |a_n|^2 <= N, V = N^sigma, T = N^tau.
W subset [0,T], 1-separated, |D(t)| >= V on W; R = |W|. MVT (Montgomery-Vaughan): for every interval I,
  int_I |D|^2 <= |I| G + c_MV sum_n n |a_n|^2,  c_MV <= 3pi (MV 1974), >= 2pi (a_n = 1, I centred at t = 0).
Halasz-Montgomery: R^2 V^2 <= G sum_{r,s} |S(t_r - t_s)|, S(u) = sum_{N<n<=2N} n^{iu}.

## 2. Localization lemma (CZ at height V^2, Montgomery scale)
If avg_I |D|^2 > V^2 > G then |I| < c_MV * 2N * G/(V^2 - G) (use sum n|a_n|^2 <= 2N G).
So every Calderon-Zygmund stopping interval of |D|^2 at height V^2 has length < L* = 2 c_MV N G/(V^2 - G)
(~ 6pi N^{2-2sigma} when G = N); the weak-(1,1) count of stopping intervals is <= 3 (T + 2 c_MV N) G / V^2 = MVT.
Proof: one line from the local MVT. Classical content (the local MVT); no new bound.

## 3. The class of real-variable arguments (CZ, maximal, Carleson embedding, sparse, stopping time)
Data N, T, V, G_1..G_K (G_k >= coefficient norm of D^k; G_k <= N^{k+o(1)} for |a_n| <= 1 by the divisor bound).
w >= 0 on R is ADMISSIBLE if (A1) supp w^ in [-log 2, log 2]; (A2) for all intervals I and k <= K:
  int_I w^k <= |I| G_k + c 2^k N^k G_k   (c = 3pi).
Every |D|^2 is admissible: D^k has frequencies log m, m in (N^k, 2^k N^k], and Montgomery-Vaughan's local MVT gives (A2).
CZ stopping intervals, M f, weak (1,1), good-lambda, John-Nirenberg, Carleson embedding of PW measures (large sieve with
measures), sparse bounds of averaging operators: all consume w only through positivity, local averages of w^k and (A1).
R(w) := max number of 1-separated t in [0,T] with w(t) >= V^2.

## 4. Theorem P (class optimum = k-th power mean value bound)
F(t) = (sin(t log2/2)/(t log2/2))^2: F >= 0, F(0) = 1, F^ supported in [-log 2, log 2], int F = 2pi/log 2 = 9.0647,
int F^2 = 6.0431, int F^3 = 4.9856. Comb overlap: (period mass of (sum_r F(.-rL))^k)/int F^k <= 1.13 for L >= 20 (t3_out).
Put rho_k = 1.13 int F^k, C_k = c 2^k N^k, B_k = (T + C_k) G_k/(rho_k V^{2k}), B = min_{k<=K} B_k, L' = 2T/B.
If L' >= 20 and 4 rho_k V^{2k} <= C_k G_k (k <= K) then w = V^2 sum_{r <= B/2} F(t - rL') is admissible, R(w) >= floor(B/2).
Proof: w <= the infinite comb, whose k-th power has period mass <= rho_k V^{2k}; a window of length H <= T meets <= H/L' + 2 periods:
  mass <= (H/L' + 2) rho_k V^{2k} <= (H/2) G_k (T + C_k)/T + 2 rho_k V^{2k} <= H G_k/2 + C_k G_k/2 + C_k G_k/2 <= H G_k + C_k G_k;
  a window with H > T: mass <= (B/2 + 2) rho_k V^{2k} <= (T + C_k) G_k/2 + C_k G_k/2 <= H G_k + C_k G_k. (A1) holds termwise.
  (Correction 20:15: the first draft used B spikes; the H -> T windows need the factor 2. t4 used B spikes and still passed, ratio <= 0.893.)
Conversely every admissible w has R(w) <= C_0 B, C_0 absolute: w = |g|^2 with g of type (log 2)/2 (Akhiezer-Krein), g^k = g^k * psi_k
(psi_k Schwartz, psi_k^ = 1 on the band) gives w(t_r)^k <= C int w^k |psi_k(t_r - u)| du; sum over 1-separated t_r and use (A2).
So the best large-values bound any class argument can prove is the k-th power MVT bound min_k (T + N^k) G_k V^{-2k}, up to
absolute constants. With G_k = N^{k+o(1)}: R_class = N^{LV_c + o(1)}, LV_c(sigma,tau) = min_k [max(tau,k) + k(1 - 2 sigma)].

## 5. Corollaries (t3_out)
(a) No class argument proves Huxley or Guth-Maynard where they beat LV_c: e.g. (sigma,tau) = (0.8,1.2): LV_c = 0.6, GM = Huxley
    = 0.4; (0.75,1.5): 1.0 vs GM 0.9; (0.9,1.5): 0.4 vs 0.2. At tau >= 2 the listed bounds do not beat LV_c.
(b) Through the GM sec. 13 / Ingham zero detector (ExploreKakeya bookkeeping, copied), LV_c gives A(sigma) = 3/(2 - sigma)
    exactly (2.0689, 2.1420, 2.2211, 2.3067, 2.4000, 2.4980 at 0.55..0.80 vs Ingham 2.0690 ... 2.5000; grid effect).
    So the class reproduces Ingham, which is the record on [1/2, 7/10]: there it cannot improve; beyond 7/10 every
    improvement (Huxley, Heath-Brown, Jutila, Bourgain, GM, TTY) uses the phase of S(u) = sum n^{iu} (Halasz-Montgomery).
(c) The HM kernel: |S(u)| <= 12 N/|u| on [1, 0.9*2piN] (a Calderon-Zygmund/Hilbert kernel, already handled losslessly by
    Montgomery-Vaughan), while for u ~ 2piN m, max|S| = (1.26..2.52) sqrt N and rms|S| = (0.79..1.00) sqrt N (N = 300, 1000,
    m <= 12): the CZ size ratio max|S| u/N reaches 6000 at N = 1000, m = 12. Not a CZ kernel; sparse domination cannot act.

## 6. Smallest cases of Theorem P (t4_profile.py, t4_out.txt)
Minimal data G_k = N^k (hardest for the profile), K = 3, sigma = 0.8, c = 3pi; sliding windows of 60 lengths in [1, T+400]:
  N = 10, tau = 3: B = 29.2, L' = 34.3, max window ratio 0.893 <= 1, R(w) = 53;  N = 30, tau = 2: B = 18.6, ratio 0.863;
  N = 100, tau = 1.5: B = 17.8, ratio 0.852;  N = 300, tau = 1.2: B = 21.0, ratio 0.881. All admissible (ratio < 1).
At fixed N the constants dominate; the barrier is the exponent statement of sec. 4-5.

## 7. Other suggested bridges, one line each
- N8 (Levinson bad set): Jensen's step <log|f|> <= (1/2) log <|f|^2> is lost only through small values of |f|. CZ, maximal and
  weak-type tools bound upper tails |{|f| > lambda}|; the constant |f| = ||f||_2 meets every Chebyshev-type upper-tail bound and
  has Jensen gap 0, so no upper-tail information improves the Jensen step. Needs lower-tail (negative moment) input. Record (67.36%)
  is not Levinson-based anyway.
- N10 (psi in short intervals): sup_{h >= h0} |psi(x+h) - psi(x) - h|/h <= M(F_{h0})(x) with F the explicit-formula derivative
  smoothed at h0; Hardy-Littlewood in L^2 costs a factor 4 against the fixed-h0 mean square (no log loss). Routine; the exponent
  (almost all x, h >= x^{2/15+eps} from Guth-Maynard) is unchanged.

## 8. Side corollary candidate (off-tree, harvest type; not reviewed)
Lacey-Mousavi-Rahimi 2101.10401 (A_N f = N^{-1} sum_{n<=N} Lambda(n) f(x-n)): sharp l^p-improving bound (Thm t:first, l.192),
sparse bound for sup_N A_N and Cor c:logell with t = 1 under GRH, t = 2 unconditionally. TeX main.tex: GRH enters only via
(e:primecount) l.185, psi(x;q,a) = x/phi(q) + O(x^{1/2} log^2 q), used once, in PNTlemma (l.414, step l.452: error Q N^{-1/2+eps}),
which l:QQ (l.497) applies only for q <= N_0 = N~^{1/400} with N~ = N^{1/5} (l.484) and needs error << N~^{-1/1000}; and through
the absence of exceptional characters. t = 2 comes only from the trivial range <f><g> << e^{-c sqrt n} forced by N~ = e^{c sqrt n/4}
(l.954-960). The Lean-checked half-plane (no zero of any L(s,chi), Re s > B_2 = 0.87338) gives psi(x;q,a) = x/phi(q) +
O(x^{B_2} log^2(qx)) uniformly, no exceptional zeros, and PNTlemma error Q N^{B_2-1+eps} <= N^{1/2000 + B_2 - 1 + eps} = N^{-0.126}
<< N~^{-1/1000} = N^{-1/5000}. Any uniform half-plane Re s > theta with theta < 1 - 1/2000 - 1/5000 = 0.9993 suffices. Reviewer (verdict-only, timed out with a summary, 12 Nov 20:12 UTC): INCORRECT AS STATED, not refuted after repair.
(1) PNTlemma needs Q >= max(q, N|xi - a/q|); small q does not give small Q, so the arcs must be narrowed to
|xi - a/q| <= N^{-theta-delta} (l.415, 514-516, 545-553), and the gap between the arcs needs Vinogradov's bound (q' in (N_0, N^theta]).
(2) Summing over residues costs phi(q)(1 + N|beta|) (l.439-461). (3) Separation inequality l.535 fails for q_0 = 1 < q.
(4) In Thm t:first the Err bound N~^{-1/1000} gives the needed (loglog Q)/Q only for Q <= N~^{1/1000}; the trivial range
<f><g> < cN^{-1/4} (l.995) must shrink to cN^{-eta}, eta ~ 1/2500 (constant 1/eta), under GRH too (l.970-998).
(5) Page-exceptional zeros beta_q > 1 - c/log q exceed theta for large q, so the Ex term vanishes.
Status: the t = 1 bound under a uniform half-plane theta < 1 is plausible after standard repairs (narrow arcs + minor-arc bound
+ eta-threshold), but LMR's own GRH proof needs the same repairs. No verbatim substitution. Do not file without a repaired proof.
