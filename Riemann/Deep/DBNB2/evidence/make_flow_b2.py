#!/usr/bin/env python3
"""dbnb2 producer: Raymond's 67-row comparison flow (Zenodo 23752103) restarted from
the Lean-checked half-plane Re s > b0 = 0.87338196661320752 instead of 7/8.

Unchanged: the 67 rows (plan M8_LATE_PLAN.json, source MAIN_CLOCK_SOURCE.json, hash-checked),
step 1/2048000, cushion 1e-12, eta 1e-12, startup tau 1e-14, 2^-96 dyadic outward rounding,
stop at the first arc end with q <= STOP (author: 1e-4).
Changed: q_init = (2 b0 - 1)^2 + 1e-12 (exact rational).
Exact FLINT rationals throughout. The author's verifier (with only its hard-coded
old-run constants replaced) checks every arc of the output.
"""
from pathlib import Path
import argparse, gzip, hashlib, json, sys, time
sys.dont_write_bytecode = True
from flint import fmpq as Q

B0 = Q(87338196661320752, 10**17)
PLAN = ('lambda_round14/material/m8_late/M8_LATE_PLAN.json', 'ec346b7c1136eda38021eb6323033b3b48ed21e6f7458eeeac7cc9579db2497a')
SOURCE = ('lambda_round15/flow_main_v2/MAIN_CLOCK_SOURCE.json', 'b59905d7b994d3424ba14f9347be75a2429ae4d4b0de2c2404971ffda1d4d536')
D = 2**96
STEP = Q(1, 2048000)
ETA = Q(1, 10**12)

def fl(x):
    return int(x.numerator)/int(x.denominator)

def floor_dy(x):
    return Q(int((x * D).floor()), D)

def ceil_dy(x):
    return Q(int((x * D).ceil()), D)

def own_via_newton(q, row):
    # identical algebra to the author's verify_arcs_fresh.own_via_newton
    result = Q(0)
    for h, v, w in row['rational_nodes']:
        hh = h*h; vv = v*v
        den = (hh-q)*(hh-q) + vv*(2*hh+2*q+vv)
        assert hh > q and den > 0
        result += 2*w*h*(hh+vv-q)/den
    if row['overlay'] is None:
        return result
    alpha, centres = row['overlay']
    result *= alpha
    for p, coefficients in centres:
        pole = p*p - q
        assert pole > 0
        root_sum = 2*p/pole; root_product = 1/pole
        previous = Q(2); current = root_sum; scale = 1/p
        for m, b in enumerate(coefficients):
            result += 2*b*scale*current*(-1 if m % 2 else 1)
            previous, current = current, root_sum*current - root_product*previous
            scale *= p
    return result

def pade_negative_upper(x):
    assert x >= 0
    positive = Q(1); negative = Q(1)
    for c in (42, 840, 10080, 75600, 332640, 665280):
        positive = positive*x + c
        negative = -negative*x + c
    assert positive > 0
    return negative/positive

def prepare_rows(plan, source):
    rows = []
    prev_end = Q(0)
    for i, (p, s) in enumerate(zip(plan['rows'], source['records'])):
        assert p['row'] == s['row'] == i and p['mode'] == s['mode']
        start, end = map(Q, p['time'])
        cap0, cap = map(Q, p['q_interval'])
        assert start == prev_end < end and cap0 == 0 < cap
        main, debit, first, last = (Q(s[k]) for k in ('main_lo', 'native_debit', 'first_hi', 'last_hi'))
        assert first >= last > 0
        intercept = main - debit - first
        slope = (first - last)/(end - start)
        nodes = [tuple(Q(n[k]) for k in ('height', 'shift', 'weight')) for n in p['nodes']]
        overlay = None
        if p.get('joint_jet_overlay') is not None:
            o = p['joint_jet_overlay']
            overlay = (Q(o['alpha']), [(Q(c['p']), list(map(Q, c['b']))) for c in o['centres']])
        rows.append(dict(start=start, end=end, cap=cap, intercept=intercept, slope=slope,
                         rational_nodes=nodes, overlay=overlay))
        prev_end = end
    assert prev_end == Q(87, 800)
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workspace', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--stop', default='1/10000')
    ap.add_argument('--b0', default=None, help='override b0 as p/q (sensitivity runs only)')
    ap.add_argument('--tube', action='store_true', help='keep the retained row-55 tube [289/10000, 361/10000] by outward (upward) resets')
    a = ap.parse_args()
    b0 = Q(a.b0) if a.b0 else B0
    stop = Q(a.stop)
    t0 = time.monotonic()
    data = {}
    for k, (rel, dig) in (('plan', PLAN), ('source', SOURCE)):
        raw = (a.workspace/rel).read_bytes()
        assert hashlib.sha256(raw).hexdigest() == dig, k
        data[k] = json.loads(raw)
    rows = prepare_rows(data['plan'], data['source'])
    q0 = (2*b0 - 1)**2
    initial = q0 + ETA
    t = Q(0); q = initial; i = 0
    recs = 0; min_cap = None; first_rate = None; resets = 0
    tmp = a.out.with_suffix('.records.tmp')
    with gzip.open(tmp, 'wt', compresslevel=6, encoding='utf-8') as stream:
        while True:
            while t + STEP > rows[i]['end']:
                i += 1
                assert i < 67, 'ran past the 67-row horizon'
            row = rows[i]
            assert 0 < q <= row['cap'], ('cap', recs, i)
            cc = row['cap'] - q
            if min_cap is None or cc < min_cap: min_cap = cc
            sv = row['intercept'] + row['slope']*(t - row['start'])
            lower = floor_dy(sv)
            upper = ceil_dy(own_via_newton(q, row))
            r = floor_dy(2*(lower - upper) - ETA)
            assert r > 0 and 2*(lower - upper) - r >= ETA
            if first_rate is None: first_rate = r
            positive = q - (2 + r*q)*STEP
            assert positive > 0, ('positivity', recs)
            ideal_up = (q + 2/r)*pade_negative_upper(r*STEP) - 2/r
            qout = ceil_dy(ideal_up)
            if a.tube and i == 55:
                L = Q(289, 10000)
                m = int(((row['end'] - (t + STEP))/STEP).floor())
                c = ceil_dy((L + 2*STEP)/(1 - 40*STEP))
                need = c + Q(m, D) if m > 0 else L
                if qout < need:
                    qout = need; resets += 1
                assert positive >= L and qout >= L and q <= Q(361, 10000) and qout <= Q(361, 10000), ('tube', recs)
            assert 0 < qout < q
            rec = {'row': i, 't': str(t), 'dt': str(STEP), 'qin': str(q), 'qout': str(qout),
                   'source_lower': str(lower), 'own_upper': str(upper), 'rate': str(r),
                   'linear_positive_lower': str(positive)}
            stream.write((',' if recs else '') + json.dumps(rec, separators=(',', ':')))
            recs += 1
            t += STEP; q = qout
            if recs % 20000 == 0:
                print('ARCS', recs, 'row', i, 't', fl(t), 'q', fl(q), 'el', round(time.monotonic()-t0, 1), flush=True)
            if q <= stop:
                break
    tau = Q(1, 10**14)
    startup = ETA - first_rate*initial*tau
    endpoint = t + q/2
    # upward 12-digit decimal ceiling, as the author's headline
    head = Q(int((endpoint*10**12).floor()) + 1, 10**12)
    assert endpoint < head < endpoint + ETA
    header = {'schema': 'MULTIPLICATIVE_PACKED_FLOW_V1', 'initial_sigma': str(b0), 'initial_q': str(initial),
              'steps': recs, 'terminal_t': str(t), 'terminal_q': str(q), 'endpoint': str(endpoint),
              'headline': str(head), 'headline_decimal': fl(head), 'minimum_kernel_cap_margin': str(min_cap),
              'startup_tau': str(tau), 'startup_gap_lower': str(startup), 'step_denominator': 2048000,
              'plan_sha256': PLAN[1], 'source_sha256': SOURCE[1],
              'producer_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    with gzip.open(a.out, 'wt', compresslevel=6, encoding='utf-8') as out:
        h = json.dumps(header, separators=(',', ':'))
        out.write(h[:-1] + ',"records":[')
        with gzip.open(tmp, 'rt', encoding='utf-8') as src:
            for chunk in iter(lambda: src.read(1 << 20), ''):
                out.write(chunk)
        out.write(']}')
    tmp.unlink()
    summary = dict(header); summary['endpoint_decimal'] = fl(endpoint); summary['terminal_t_decimal'] = fl(t)
    summary['terminal_q_decimal'] = fl(q); summary['startup_gap_decimal'] = fl(startup)
    summary['first_rate_decimal'] = fl(first_rate); summary['initial_q_decimal'] = fl(initial)
    summary['seconds'] = time.monotonic() - t0
    summary['tube_upward_resets'] = resets
    (a.out.parent/(a.out.name.replace('.json.gz', '') + '_SUMMARY.json')).write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({k: summary[k] for k in ('steps', 'terminal_t', 'terminal_t_decimal', 'terminal_q_decimal', 'endpoint_decimal', 'headline', 'startup_gap_decimal', 'first_rate_decimal', 'seconds')}, indent=2))

if __name__ == '__main__':
    main()
