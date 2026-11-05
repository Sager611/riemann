def isprime(n):
    if n < 2: return False
    if n % 2 == 0: return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0: return False
        d += 2
    return True
C = [x for x in range(1, 4006) if isprime(2*x - 3) and isprime(2*x + 3)]
Cs = set(C)
reps = [(x, 4006 - x) for x in C if (4006 - x) in Cs]
print('C has', len(C), 'elements below 4006; first:', C[:12])
print('representations of 4006 as c1+c2 (c1,c2 in C):', reps)
for n in (4005, 4007, 4008):
    print(n, [(x, n - x) for x in C if (n - x) in Cs][:3])
# which n>=4006 up to 20000 are missing?
C2 = [x for x in range(1, 20001) if isprime(2*x - 3) and isprime(2*x + 3)]
S = set(C2)
miss = [n for n in range(4000, 20001) if not any((n - x) in S for x in C2 if x < n)]
print('missing n in [4000,20000]:', miss)
