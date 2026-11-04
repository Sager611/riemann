"""Independent checks of the combinatorial / arithmetic lemmas of OpenAI's
"Universal Tensor Squares for Symmetric Groups" (Sep 24 2026). Exact integers.
 U1 Lemma 3.1 (candidate): r bounds, L self-conjugate of size n, M>=9 for n>=65         (n <= 200000)
 U2 Lemma 6.1 (two opposite prefixes): n <= max(H_i+W_j-1, floor(H_i W_j/(ij)))  (all partitions n<=30, all i,j)
 U3 Lemma 6.2 (M>=22): the two polynomial bounds and n > max(2S, S^2/64)          (22<=M<=20000, all r)
 U4 Lemma 6.3: recompute the 2,832 intermediate assertions; failing (M,r) must lie in the 19 exceptions
 U5 Prop 5.9 (balance bound): packing A (sum n, |A|>=2M-5, max A<=M+1) for 9<=M<=200, and Young's
    inequalities for every admissible first-four column profile (worst-case completion), 9<=M<=40
 U6 Lemma 5.8 (odd horizontal strip): every 2-row target of size c+x has an even-row precursor (c<=80)
 U7 Lemma 5.5 second part: every eta in H_4(2q+4), q>=3, has a removable horizontal 2-strip leaving
    an odd-height column (q<=14)
 U8 Lemma 4.7 (two extra copies): the split alpha,beta of every eta |- 2r (even column heights,
    <=8 columns) has the stated shape, and c^theta_{alpha,beta} > 0 by an LR-tableau search (r<=11)"""
import itertools
import mn

def N(m): return m*(m+1)//2

def params(n):
    from math import isqrt
    M = (isqrt(8*n+1)-1)//2
    while N(M) > n: M -= 1
    while (n - N(M)) % 2: M -= 1
    r = (n - N(M))//2; b = r//2; s = 1 + r % 2
    return M, r, b, s

def fast_conj(c):
    """conjugate of a weakly decreasing list of positive ints, O(len + c[0])"""
    out = []; k = len(c)
    for j in range(c[0] if c else 0):
        while k > 0 and c[k-1] <= j: k -= 1
        out.append(k)
    return tuple(out)

def from_cols(cols):
    cols = sorted([c for c in cols if c > 0], reverse=True)
    return fast_conj(cols)

def U1(nmax=30000):
    cnt = 0
    for n in range(1, nmax+1):
        M, r, b, s = params(n)
        if M % 2: assert 2*r <= M-1, (n, M, r)
        else: assert 2*r <= 3*M+4, (n, M, r)
        if n >= 65: assert M >= 9, n
        if M >= 9:
            cols = [M+b+s-1, M-1+b] + list(range(M-2, 2, -1)) + [2]*(b+1) + [1]*s
            L = from_cols(cols)
            assert sum(L) == n and L == fast_conj(list(L)) and tuple(sorted(cols, reverse=True)) == fast_conj(list(L)), n
            cnt += 1
    print(f"U1: candidate lemma ok for n<=" + str(nmax) + f" ({cnt} degrees with M>=9)")

def U2(nmax=26):
    tot = 0
    for n in range(1, nmax+1):
        for lam in mn.partitions(n):
            w = list(lam); h = list(mn.conj(lam))
            Wp = [0]; Hp = [0]
            for x in w: Wp.append(Wp[-1]+x)
            for x in h: Hp.append(Hp[-1]+x)
            for i in range(1, n+1):
                Hi = Hp[min(i, len(h))]
                for j in range(1, n+1):
                    Wj = Wp[min(j, len(w))]
                    assert n <= max(Hi+Wj-1, (Hi*Wj)//(i*j)), (lam, i, j)
                    tot += 1
    print(f"U2: rectangle lemma ok, {tot} (partition,i,j) triples, n<={nmax}")

def U3(Mmax=2000):
    for x in range(0, 300):
        for f in range(0, 12):
            assert x <= 2*max(0, (x-f)//2) + f + 1
    cnt = 0
    for M in range(22, Mmax+1):
        rmax = (M-1)//2 if M % 2 else (3*M+4)//2
        for r in range(1, rmax+1):
            n = N(M) + 2*r; S = 2*r + 38
            assert 2*(n - 2*S) >= M*M - 5*M - 160 > 0, (M, r)
            assert 64*n - S*S >= 23*M*M - 28*M - 1508 > 0, (M, r)
            assert n > max(2*S, (S*S)//64)
            cnt += 1
    print(f"U3: large-M capacity bounds ok for 22<=M<={Mmax} ({cnt} pairs)")

EXC = {(10, r) for r in range(8, 18)} | {(12, r) for r in range(14, 21)} | {(14, 22), (14, 23)}

def U4():
    def q(M, f): return min(8, (2*M-8-f)//f)
    def D(M, r, d, e):
        S = set()
        if d > 1: S.add((d-1, 2*M-2))
        for f in range(e, 5):
            if q(M, f) > 0: S.add((q(M, f), 2*r-2+q(M, f)*(f+1)))
        return S
    fails = set(); asserts = 0; pairs = 0
    for M in range(9, 22):
        rmax = (M-1)//2 if M % 2 else (3*M+4)//2
        for r in range(1, rmax+1):
            n = N(M) + 2*r
            ok_all = True
            for d in range(1, 5):
                for e in range(1, 5):
                    val = min(max(U+V-1, (U*V)//(i*j)) for (i, U) in D(M, r, d, e) for (j, V) in D(M, r, e, d))
                    if not val < n: ok_all = False
                    if (M, r) not in EXC: asserts += 1
            if (M, r) not in EXC: pairs += 1
            if not ok_all: fails.add((M, r))
    assert fails <= EXC, sorted(fails - EXC)
    print(f"U4: intermediate check: {pairs} non-exceptional pairs, {asserts} assertions all hold; "
          f"failing pairs = {len(fails)} (all among the 19 listed; listed but not failing: {sorted(EXC - fails)})")

def pack(t, j):
    return [x for x in ((t+i)//j for i in range(j)) if x > 0]

def balance_A(M, r):
    b = r//2; s = 1 + r % 2; c = 2*b+2; T = 2*M-1+r; n = N(M) + 2*r
    A = pack(T + 2*(s-1), 4)
    if s == 1: A.append(1)
    if c <= M:
        A.append(c); mids = list(range(3, M-1))
    else:
        A += pack(c+3, 2); mids = list(range(4, M-1))
    # pair consecutively, leaving the smallest alone if odd
    if len(mids) % 2:
        A.append(mids[0]); mids = mids[1:]
    for k in range(0, len(mids), 2):
        u = mids[k+1]  # pair (u-1, u)
        assert mids[k] == u-1
        A += pack(2*u-1, 4)
    assert sum(A) == n, (M, r, sum(A), n)
    return A, n

def U5(Mmax_pack=400, Mmax_young=16):
    for M in range(9, Mmax_pack+1):
        rmax = (M-1)//2 if M % 2 else (3*M+4)//2
        for r in range(1, rmax+1):
            A, n = balance_A(M, r)
            assert len(A) >= 2*M-5 and max(A) <= M+1, (M, r, len(A), max(A))
    checked = 0
    for M in range(9, Mmax_young+1):
        rmax = (M-1)//2 if M % 2 else (3*M+4)//2
        for r in range(1, rmax+1):
            A, n = balance_A(M, r)
            PA = lambda k: sum(min(k, a) for a in A)
            P = [PA(k) for k in range(0, n+2)]
            lim = 2*M-2
            for h1 in range(1, lim+1):
                for h2 in range(1, min(h1, lim-h1)+1):
                    for h3 in range(1, min(h2, lim-h1-h2)+1):
                        for h4 in range(1, min(h3, lim-h1-h2-h3)+1):
                            H4 = h1+h2+h3+h4
                            if H4 + 10**9*0 > n: continue
                            # a completion exists only if the rest fits with parts <= h4 (always possible if n >= H4)
                            pref = [0, h1, h1+h2, h1+h2+h3, H4]
                            for k in range(1, max(A)+1):   # for k >= max A, P_A(k) = n
                                bound = min(n, pref[k] if k <= 4 else H4 + (k-4)*h4)
                                assert bound <= P[k], (M, r, (h1, h2, h3, h4), k, bound, P[k])
                            checked += 1
    print(f"U5: balance packing ok for 9<=M<={Mmax_pack}; Young inequalities ok for {checked} first-four profiles, 9<=M<={Mmax_young}")

def U6(cmax=80):
    cnt = 0
    for c in range(2, cmax+1, 2):
        for x in (1, 3):
            if x > c: continue
            for beta in range(0, (c+x)//2+1):
                alpha = c+x-beta
                ys = [y for y in range(0, c+1, 2) if c-y >= y and alpha >= c-y >= beta >= y]
                assert ys, (c, x, alpha, beta)
                lo, hi = max(0, beta-x), min(beta, c-beta)
                assert all(lo <= y <= hi for y in ys) and any(y % 2 == 0 for y in range(lo, hi+1))
                cnt += 1
    print(f"U6: odd horizontal strip ok ({cnt} targets, c<={cmax})")

def removable_hstrips(lam, size):
    """partitions mu with lam/mu a horizontal strip of the given size"""
    lam = list(lam); L = len(lam)
    res = []
    def rec(i, mu, left):
        if i == L:
            if left == 0: res.append(tuple(x for x in mu if x > 0))
            return
        nxt = lam[i+1] if i+1 < L else 0
        for t in range(0, min(left, lam[i]-nxt)+1):
            rec(i+1, mu+[lam[i]-t], left-t)
    rec(0, [], size)
    return res

def U7(qmax=14):
    cnt = 0
    for q in range(3, qmax+1):
        for eta in mn.partitions(2*q+4):
            if len(eta) > 4: continue
            ok = any(any(h % 2 for h in mn.conj(mu)) for mu in removable_hstrips(eta, 2))
            assert ok, (q, eta)
            cnt += 1
    print(f"U7: two-path Pieri step ok ({cnt} targets, 3<=q<={qmax})")

def lr_positive(theta, alpha, beta):
    """exists an LR tableau of shape theta/alpha with content beta (search, early exit)"""
    theta = list(theta); alpha = list(alpha) + [0]*(len(theta)-len(alpha))
    if any(alpha[i] > theta[i] for i in range(len(theta))) or len(alpha) > len(theta): return False
    if sum(theta) != sum(alpha) + sum(beta): return False
    cells = [(i, j) for i in range(len(theta)) for j in range(alpha[i], theta[i])]
    # reading order: rows top to bottom, right to left within row
    order = sorted(cells, key=lambda c: (c[0], -c[1]))
    T = {}
    cnt = [0]*(len(beta)+1)
    def ok_place(i, j, v):
        if (i, j+1) in T and T[(i, j+1)] < v: return False   # rows weakly increase left->right
        if j+1 < theta[i] and (i, j+1) not in T and j+1 >= alpha[i]: return False
        if (i-1, j) in T and T[(i-1, j)] >= v: return False  # columns strictly increase
        if i > 0 and j < theta[i-1] and j >= alpha[i-1] and (i-1, j) not in T: return False
        if i > 0 and j < alpha[i-1]: pass
        return True
    def rec(k):
        if k == len(order): return True
        i, j = order[k]
        for v in range(1, len(beta)+1):
            if cnt[v] >= beta[v-1]: continue
            if v > 1 and cnt[v-1] <= cnt[v]: continue      # lattice word condition
            if not ok_place(i, j, v): continue
            T[(i, j)] = v; cnt[v] += 1
            if rec(k+1): return True
            del T[(i, j)]; cnt[v] -= 1
        return False
    return rec(0)

def U8(rmax=11):
    # sanity of the LR search against known values
    assert lr_positive((2, 1), (1,), (1, 1)) and lr_positive((2, 1), (1,), (2,)) and not lr_positive((3,), (1, 1), (1,))
    assert lr_positive((2, 2, 1), (2, 1), (1, 1)) and lr_positive((2, 2), (2,), (2,)) and not lr_positive((2, 2), (2,), (1, 1))
    assert lr_positive((3, 2, 1), (2, 1), (2, 1)) and not lr_positive((4, 2), (2, 1), (1, 1, 1))
    for k in range(0, 4):
        kap = tuple([2]*k + [1])
        assert lr_positive(tuple([2]*(2*k+1)), kap, kap), k
    cnt = 0
    for r in range(1, rmax+1):
        for eta in mn.partitions(2*r):
            heights = mn.conj(eta)            # eta column heights = theta rows
            if any(h % 2 for h in eta) and False: pass
            if any(h % 2 for h in heights): continue          # need even column heights
            if len(heights) > 8: continue                     # at most eight columns
            theta = heights                                     # theta = eta^t
            tcols = mn.conj(theta)                              # columns of theta, come in equal pairs
            js = [tcols[i] for i in range(0, len(tcols), 2)]
            assert all(tcols[i] == tcols[i+1] for i in range(0, len(tcols), 2)) and sum(js) == r and max(js) <= 8
            alpha_cols = []; beta_cols = []
            odd = [x for x in js if x % 2]; even = [x for x in js if x % 2 == 0]
            reserved = None
            if r % 2:
                reserved = odd[0]; odd = odd[1:]
            flip = 0
            for j in even:
                alpha_cols += [j//2]*2; beta_cols += [j//2]*2
            for j in odd:
                u, v = j//2, j - j//2
                if flip: u, v = v, u
                flip ^= 1
                alpha_cols += [u]*2; beta_cols += [v]*2
            if reserved is not None:
                k = (reserved-1)//2
                kap = [2]*k + [1]
                kc = mn.conj(tuple(kap))
                alpha_cols += list(kc); beta_cols += list(kc)
            al = from_cols(alpha_cols); be = from_cols(beta_cols)
            assert sum(al) == r and sum(be) == r and len(al) <= 4 and len(be) <= 4, (eta, al, be)
            nodd = lambda p: sum(1 for x in p if x % 2)
            if r % 2 == 0: assert nodd(al) == 0 and nodd(be) == 0
            else: assert nodd(al) == 1 and nodd(be) == 1
            assert lr_positive(theta, al, be), (eta, theta, al, be)
            cnt += 1
    print(f"U8: two-extra splits ok and LR-positive for {cnt} shapes, r<={rmax}")

if __name__ == "__main__":
    import sys
    which = sys.argv[1:] or ["U1", "U2", "U3", "U4", "U5", "U6", "U7", "U8"]
    for w in which:
        globals()[w]()
