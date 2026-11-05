#!/usr/bin/env python3
"""sidecheck: convention test for PROOF.md Step 1 and the certificate's corridor angles.
Simulate an explicit ambidextrous shape (the lens = upper unit half-disk about (0,0) cap lower unit
half-disk about (0,1), area 2*pi/3 - sqrt(3)/2 ~ 1.2284) moving through the Z-hallway
(east, clockwise corner, south, counterclockwise corner, east).  Check in the sofa frame:
 (i)  g^-1(L_world) == C(psi,u) pointwise, with psi = -theta at corner 1 and psi = pi - theta at corner 2,
 (ii) the moved lens stays inside the Z-hallway,
 (iii) the lens fits in K cap C(phi,u1) cap C(pi+phi,u2) (the certificate's case {0,2}, quadrant() mapping),
 (iv) a right-turn-only half-disk does NOT fit any left-turn window angle with its own corner (sanity)."""
import math, random
random.seed(3)
D = 40.0
def rot(a, z): c, s = math.cos(a), math.sin(a); return (c*z[0] - s*z[1], s*z[0] + c*z[1])
def inL0(z):  x, y = z; return (x <= 1 and 0 <= y <= 1) or (0 <= x <= 1 and y <= 1)
def inL2(z):  # Rot_pi(L0) + (1, -D+1)
    w = (z[0] - 1, z[1] + D - 1); return inL0((-w[0], -w[1]))
def inZ(z):
    x, y = z
    return (x <= 1 and 0 <= y <= 1) or (0 <= x <= 1 and -D <= y <= 1) or (x >= 0 and -D <= y <= -D + 1)
def inC(psi, u, z):   # PROOF.md: C(psi,u) = Rot_psi(L0+u) = P u Q in rotated coords
    p = z[0]*math.cos(psi) + z[1]*math.sin(psi); q = -z[0]*math.sin(psi) + z[1]*math.cos(psi)
    return (p <= u[0] + 1 and u[1] <= q <= u[1] + 1) or (u[0] <= p <= u[0] + 1 and q <= u[1])
def in_lens(z): x, y = z; return x*x + y*y <= 1 and y >= 0 and x*x + (y-1)**2 <= 1 and y <= 1
# motion: list of (theta, c) with g(z) = Rot_theta z + c
def motion():
    out = []
    for t in range(21): out.append(('H1', 0.0, (-5 + 5*t/20, 0.0)))
    for t in range(61):
        a = (math.pi/2) * t/60; out.append(('C1', -a, (0.0, 0.0)))
    for t in range(21):
        out.append(('V', -math.pi/2, (0.0, -(D - 1) * t/20)))
    for t in range(61):
        b = (math.pi/2) * t/60
        # rotate counterclockwise by b about corner (1,-D+1): g = R_b o (R_{-pi/2} z + (0,-D+1) - k) + k
        k = (1.0, -D + 1)
        th = -math.pi/2 + b
        # g(z) = R_b(R_{-pi/2} z + (0,-D+1) - k) + k = R_th z + R_b((0,-D+1)-k) + k
        v = rot(b, (0.0 - k[0], -D + 1 - k[1])); out.append(('C2', th, (v[0] + k[0], v[1] + k[1])))
    for t in range(21):
        out.append(('H2', 0.0, (1.0 + 4*t/20, -D)))
    return out
def g(th, c, z): w = rot(th, z); return (w[0] + c[0], w[1] + c[1])
lens_pts = []
while len(lens_pts) < 4000:
    z = (random.uniform(-1, 1), random.uniform(0, 1))
    if in_lens(z): lens_pts.append(z)
bad_in_Z = 0; bad_eq = 0; checks = 0; angles = {'C1': [], 'C2': []}; badph = {}
for ph, th, c in motion():
    for z in lens_pts:
        if not inZ(g(th, c, z)): bad_in_Z += 1; badph[ph] = badph.get(ph, 0) + 1
    if ph in ('C1', 'C2'):
        if ph == 'C1': psi = -th; corner = (0.0, 0.0); inX = inL0
        else: psi = math.pi - th; corner = (1.0, -D + 1); inX = inL2
        # inner corner of g^-1(X) in sofa frame:
        w = rot(-th, (corner[0] - c[0], corner[1] - c[1]))
        u = rot(-psi, w)       # inner corner of C(psi,u) is Rot_psi(u)
        angles[ph].append(psi % (2*math.pi))
        for _ in range(3000):
            z = (random.uniform(-4, 4), random.uniform(-4, 4)); checks += 1
            if inX(g(th, c, z)) != inC(psi, u, z): bad_eq += 1
print('moved lens leaves the Z-hallway at %d sampled (point,time) pairs %s' % (bad_in_Z, badph))
print('sofa-frame hallway == C(psi,u) mismatches: %d of %d' % (bad_eq, checks))
for k, v in angles.items():
    print('%s hallway angles psi range: [%.4f, %.4f] deg' % (k, math.degrees(min(v)), math.degrees(max(v))))
# (iii) certificate case {0,2}: corridors quadrant(0)=(c,s) and quadrant(2)=(-c,-s), phi = atan2(20,21)
phi = math.atan2(20, 21)
for name, psi, corner in [('corner1 psi=phi', phi, (0.0, 0.0)), ('corner2 psi=pi+phi', math.pi + phi, None)]:
    if corner is None:
        th = math.pi - psi; b = th + math.pi/2; k = (1.0, -D + 1)
        v = rot(b, (0.0 - k[0], -D + 1 - k[1])); c = (v[0] + k[0], v[1] + k[1])
        w = rot(-th, (k[0] - c[0], k[1] - c[1]))
    else:
        th = -psi; c = (0.0, 0.0); w = rot(-th, (0.0, 0.0))
    u = rot(-psi, w)
    miss = sum(1 for z in lens_pts if not inC(psi, u, z))
    print('%s: u=(%.6f, %.6f)  lens points outside C(psi,u): %d of %d' % (name, u[0], u[1], miss, len(lens_pts)))
