import sys
sys.path.insert(0, '.')
mode = sys.argv[1]; n = int(sys.argv[2]); r = int(sys.argv[3]) if len(sys.argv) > 3 else 1
NODE_LIMIT = int(sys.argv[4]) if len(sys.argv) > 4 else 3 * 10**7
N = n * n
if r == 1:
    BITS = N
    FULL = (1 << N) - 1
    elems = [x for x in range(N) if x % n]
    def shift(B, e):
        return ((B << e) | (B >> (N - e))) & FULL if e else B
    GOOD = 0
    for x in range(n, N, n): GOOD |= 1 << x
else:
    BITS = N * N
    FULL = (1 << BITS) - 1
    elems = [(a, b) for a in range(N) for b in range(N) if a % n or b % n]
    LO = {}; HI = {}
    for s in range(N):
        lo = 0; hi = 0
        for x in range(N):
            for y in range(N):
                if y < N - s: lo |= 1 << (x * N + y)
                else: hi |= 1 << (x * N + y)
        LO[s] = lo; HI[s] = hi
    def shift(B, e):
        a, b = e
        if b:
            B = ((B & LO[b]) << b) | ((B & HI[b]) >> (N - b))
        if a:
            B = ((B << (a * N)) | (B >> (BITS - a * N))) & FULL
        return B
    GOOD = 0
    for a in range(0, N, n):
        for b in range(0, N, n):
            if a or b: GOOD |= 1 << (a * N + b)

# modes: 'w' any nonempty subset; 'i' subsets of size n or n+1; 's' subsets of size exactly n
if mode == 'w':
    J = None
else:
    J = n + 1 if mode == 'i' else n
best = [0, None]
nodes = [0]
E = len(elems)
def is_good(state):
    if mode == 'w': return state & GOOD
    if mode == 'i': return (state[n] | state[n + 1]) & GOOD
    return state[n] & GOOD
def add(state, e):
    if mode == 'w':
        return state | shift(state, e)
    new = list(state)
    for j in range(J, 0, -1):
        new[j] = state[j] | shift(state[j - 1], e)
    return new
seq = []
def dfs(start, state, depth):
    nodes[0] += 1
    if nodes[0] > NODE_LIMIT: raise TimeoutError
    if depth > best[0]:
        best[0] = depth; best[1] = list(seq)
    for i in range(start, E):
        e = elems[i]
        st2 = add(state, e)
        if is_good(st2): continue
        seq.append(e)
        dfs(i, st2, depth + 1)
        seq.pop()
init = 1 if mode == 'w' else [1] + [0] * J
try:
    dfs(0, init, 0); status = 'exhaustive'
except TimeoutError:
    status = 'node-limit'
print('mode', mode, 'n', n, 'r', r, status, 'nodes', nodes[0], 'max bad size', best[0], 'example', best[1])
