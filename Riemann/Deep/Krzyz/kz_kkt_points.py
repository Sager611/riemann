#!/usr/bin/env python3
"""Dump numerically found KKT points (solutions of the reduced moment system), normalized
(sorted angles in [0,2pi), conjugated if theta_1+theta_N > 2pi), with |a_n|.  usage: n tries out.json"""
import sys, json
import kz_soundness as S
n, tries, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
sols = S.find(n, tries, 777 + n)
json.dump(sols, open(out, 'w'))
for s in sols: print(round(s['an_direct'], 9), round(s['t'], 5), round(s['q'], 6), [round(x, 4) for x in s['theta']], [round(x, 4) for x in s['lam']])
