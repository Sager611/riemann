#!/usr/bin/env python3
"""q3b_intervals.py -- explore16: the q-set where Cherednik's RH holds for T(6,7) and T(7,8) (q3_daha.py found it is not an
interval). Prints the maximal q-intervals in (0, 5] (grid step 0.005) on which all 2 delta zeros of P_q lie on |u| = 1.
Run: uv run --no-project --with numpy python3 -u q3b_intervals.py
"""
from q3_daha import qt_catalan, unit_circle_count

def main():
    for n in (6, 7):
        C = qt_catalan(n); d = n * (n - 1) // 2
        coef = {}
        for (a, b), c in C.items():
            key = (d - a + b, d - a - b); coef[key] = coef.get(key, 0) + c
        grid = [k / 200 for k in range(1, 1001)]
        ok = [unit_circle_count(coef, d, q, 30000) == 2 * d for q in grid]
        iv = []; start = None
        for q, o in zip(grid, ok):
            if o and start is None: start = q
            if not o and start is not None: iv.append((start, q - 0.005)); start = None
        if start is not None: iv.append((start, grid[-1]))
        print(f'T({n},{n + 1}): RH on q-grid intervals {iv}')

if __name__ == '__main__':
    main()
