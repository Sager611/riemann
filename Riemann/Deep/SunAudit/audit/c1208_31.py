import sys, math
def first_fail(A, B, NM=80):
    u = [0, 1]
    for n in range(1, NM + 3): u.append(A*u[n] - B*u[n-1])
    if any(x <= 0 for x in u[1:]): return 'nonpos'
    l = [None] + [math.log(u[n]) / n for n in range(1, NM + 3)]
    for n in range(4, NM):
        d2 = l[n+2] - 2*l[n+1] + l[n]
        if not d2 < -1e-12: return n
    return None
def exact_fail(A, B, n):
    u = [0, 1]
    for k in range(1, n + 3): u.append(A*u[k] - B*u[k-1])
    # r_n <= r_{n+1}  <=>  u_{n+1}^{2 n (n+2)} <= u_n^{(n+1)(n+2)} u_{n+2}^{n (n+1)}
    return u[n+1]**(2*n*(n+2)) <= u[n]**((n+1)*(n+2)) * u[n+2]**(n*(n+1))
res = {}
for A in range(2, 31):
    worst = None
    Bmin = -9 if A == 2 else -1500
    for B in range(-1, Bmin - 1, -1):
        f = first_fail(A, B)
        if f is not None:
            ok = exact_fail(A, B, f) if isinstance(f, int) else None
            worst = (B, f, 'exact-confirmed' if ok else 'NOT-confirmed'); break
    res[A] = worst
print('largest failing B per A (B in [-1500,-1]):')
for A, w in res.items(): print(' A=%d:' % A, w)
pos = [(A, B, first_fail(A, B)) for A in range(2, 31) for B in range(1, (A*A + 3)//4) if A*A > 4*B and first_fail(A, B) is not None]
print('failures with B>0:', pos[:10], len(pos))
u = [0, 1]
for k in range(1, 9): u.append(3*u[k] + 41*u[k-1])
print('A=3,B=-41: u_1..u_8 =', u[1:9])
