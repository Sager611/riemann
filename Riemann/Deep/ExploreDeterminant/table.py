#!/usr/bin/env python3
"""Parse spacing2 logs into a markdown table. usage: table.py log [log ...]"""
import sys, re
rows = []
for fn in sys.argv[1:]:
    cur = {}
    for line in open(fn):
        m = re.match(r'alpha=([\d.]+) log10t=([\d.]+) nu=([\d.]+).*?H=([\d.]+).*?kept\(odd\)=(\d+)', line)
        if m: cur = dict(alpha=m[1], t=m[2], nu=m[3], H=float(m[4]), I=int(m[5])); continue
        m = re.match(r'TRUE\s+diag=(\d+) off: x3=(\d+) x2x3=(\d+) x2x3x4=(\d+) x1x2x3=(\d+) all4=(\d+)', line)
        if m: cur.update(x23=int(m[3]), all4=int(m[6])); continue
        m = re.match(r'CONIC .*max_per_conic=(\d+)', line)
        if m: cur['mpc'] = int(m[1]); continue
        m = re.match(r'SHAPE V=([\d.]+) x2x3/HWshape=([\d.]+) all4/\(HWshape\*H\^-2\)=([\d.]+)', line)
        if m: cur.update(s23=m[2], s4=m[3]); continue
        m = re.match(r'RATIO off_all4/diag=([\d.e+-]+|nan|0)\s+true/rand: x2x3=([\d.]+|nan) x2x3x4=([\d.]+|nan) all4=([\d.]+|nan|inf)', line)
        if m: cur.update(od=m[1], r23=m[2], r4=m[4]); rows.append(cur); cur = {}
print("| log10 t | ν | H | arcs | (x2,x3) pairs | ÷[H-W] shape | all-4 pairs | ÷(shape·H⁻²) | true/random (x2x3, all4) | off/diag | max per generic conic |")
print("|---|---|---|---|---|---|---|---|---|---|---|")
for r in rows:
    print(f"| {r['t']} | {r['nu']} | {r['H']:.2f} | {r['I']} | {r['x23']} | {r.get('s23','')} | {r['all4']} | {r.get('s4','')} | {r['r23']}, {r['r4']} | {float(r['od']):.3g} | {r.get('mpc','')} |")
