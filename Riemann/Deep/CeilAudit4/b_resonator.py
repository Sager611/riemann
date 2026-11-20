# B: 2606.15635 (long resonator, L(1/2, chi_{8q}), q prime) under a zero-free half-plane Re s > theta
import math
theta = 0.8733458561067185768
# (i) term-by-term: Lemma logd error X^{theta+eps} log(c+2) summed over AFE l <= X^{1/2+eps} with weight l^{-1/2}: X^{theta+1/4}; times |R| <= N = X^nu
print("(i) AFE-summed error exponent theta+1/4 =", theta+0.25, "> 1: fails already for |R| = 1 (nu = 0); GRH: 3/4")
# unbalanced AFE: lengths L1*L2 = X (root number +1), error X^theta (sqrt L1 + sqrt L2) >= 2 X^{theta+1/4}
for y in [0.0, 0.1, 0.25]:
    print(f"   unbalanced AFE, L1 = X^(1/2+{y}), L2 = X^(1/2-{y}): error exponent {theta + (0.5+y)/2:.4f}")
# (ii) twist-uniform quadratic large sieve (Heath-Brown 1995): sum_{l~L1} |E_l|^2 << X^eps (L1 + X) X log X
for a in [0.1, 0.25, 0.5]:
    # dyadic block l ~ X^a: L1^{-1/2} * sum |E_l| <= L1^{-1/2} * (L1 * (L1 + X) X)^{1/2}
    L1 = a; ex = -L1/2 + 0.5*(L1 + max(L1, 1.0) + 1.0)
    print(f"(ii) large sieve, l ~ X^{a}: block error exponent {ex:.4f} (needs < 1 - nu)")
# (iii) untwisted first moment (m = n = 1): Heath-Brown's family zero density sum_{l<=L} N(s,T,chi_l) << (LT)^eps L^{3(1-s)/(2-s)} T^{(3-2s)(1-s)/(2-s)}
best = max((s - 0.25 + 0.5*3*(1-s)/(2-s), s) for s in [0.5 + i*(theta-0.5)/100000 for i in range(100001)])
print(f"(iii) untwisted first moment, l <= X^(1/2): max_s [s - 1/4 + 3(1-s)/(2(2-s))] = {best[0]:.5f} at s = {best[1]:.4f} (< 1: transfers, error X^0.80)")
# size of resonator moduli: log m <= sum_k a log N/(k^2 log3 N) * (k + 1 + log(log N log2 N))
for logN in [1e2, 1e4, 1e8, 1e16]:
    l2 = math.log(logN); l3 = math.log(l2); a = 1.0; delta = 0.9
    K = int(l2**delta)
    lm = sum(a*logN/(k*k*l3)*(k + 1 + math.log(logN*l2)) for k in range(1, K+1))
    print(f"   log N = {logN:.0e}: max log m / log N = {lm/logN:.2f} (grows like log2 N / log3 N): twist conductor m n = N^(that)")
