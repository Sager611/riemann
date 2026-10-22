#!/usr/bin/env python3
"""Counterexample to Sun's circular-permutation conjecture (as stated in arXiv 2610.06289, Conjecture 1.1):
q = 11, a0 = 2. G_{q,a0}: vertices F_q^*, x ~ y iff x != y and a0 + x*y is a primitive root mod q.
Claim: G_{11,2} has no Hamilton cycle. Prints the neighbourhoods and runs an exhaustive search."""
q, a0 = 11, 2
phi = q - 1
P = {g for g in range(1, q) if all(pow(g, phi // p, q) != 1 for p in (2, 5))}
assert P == {2, 6, 7, 8}
V = list(range(1, q))
adj = {x: {y for y in V if y != x and (a0 + x * y) % q in P} for x in V}
print("primitive roots mod 11:", sorted(P))
for x in V:
    print(x, sorted(adj[x]))
assert adj[2] == {3, 8} and adj[9] == {3, 8}
def ham():
    path, used = [1], {1}
    def dfs():
        if len(path) == len(V):
            return 1 in adj[path[-1]]
        for y in adj[path[-1]]:
            if y not in used:
                used.add(y); path.append(y)
                if dfs():
                    return True
                path.pop(); used.discard(y)
        return False
    return dfs()
assert not ham()
print("no Hamilton cycle in G_{11,2}: verified exhaustively")
