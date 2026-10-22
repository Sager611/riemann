"""Step 7 (dbn5): Planat's UNCHANGED producer and auditor on sources whose floors are raised to our new
certified paired-probe floors (max of Planat's refined floor and ours, per provider, cell and shift).
Inputs : .scratch/dbn4v/evidence/SOURCE_*.json (dbn4 evidence sources: Planat refined + 786 kernels; sha checked)
         EVIDENCE2/floors_dbn5.json (our certified floors, one per (cell, height, shift); see run5.py / check5.py)
Outputs: EVIDENCE2/07_planat_run_dbn5.txt, EVIDENCE2/classical_dbn5.json.gz, EVIDENCE2/classical_dbn5_AUDIT.json
Scratch: .scratch/dbn5/src (rebuilt sources), .scratch/dbn5/hist"""
import gzip
import hashlib
import json
import shutil
import signal
import sys
import time
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
E2 = HERE.parent
REPO = E2.parents[3]
PLANAT = REPO / '.scratch/dbn4v/arxiv/anc/certificate'
EV4 = REPO / '.scratch/dbn4v/evidence'
WORK = REPO / '.scratch/dbn5'
SRC = WORK / 'src'
HIST = WORK / 'hist'
NAMES = ['SOURCE_OLD_ENHANCED.json', 'SOURCE_RETUNED.json']
EVIDENCE_SHA256 = ['235c34ae1c0efb03458a39687685451f9b2e0febbaab3d2d6ef38cf9d9cad714',
                   '29840261b1fee53867a61ab32c4ab1f4b81ecc588c93463b0ef28d02fcb76058']
SH = ['9/20', '19/40', '49/100', '1/2', '51/100', '21/40', '107/200', '113/200', '23/40', '3/5']
LINES = []


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    LINES.append(s)


def require(cond, msg):
    if not cond:
        raise SystemExit('FAILED: ' + msg)


def read_bytes(p):
    try:
        return Path(p).read_bytes()
    except OSError as e:
        raise SystemExit('cannot read %s: %s' % (p, e))


def load_json(p):
    try:
        if str(p).endswith('.gz'):
            with gzip.open(p, 'rt') as f:
                return json.load(f)
        return json.loads(read_bytes(p))
    except (OSError, ValueError) as e:
        raise SystemExit('cannot parse %s: %s' % (p, e))


def write_bytes(p, b):
    try:
        Path(p).write_bytes(b)
    except OSError as e:
        raise SystemExit('cannot write %s: %s' % (p, e))


def sha(p):
    return hashlib.sha256(read_bytes(p)).hexdigest()


def canon(x):
    return (json.dumps(x, separators=(',', ':'), sort_keys=True) + '\n').encode()


def planat_module(name):
    rel = 'code/' + name + '.py'
    want = load_json(PLANAT / 'SHA256SUMS.json')[rel]['sha256']
    got = sha(PLANAT / rel)
    require(got == want, 'Planat module changed: ' + rel)
    sys.dont_write_bytecode = True
    if str(PLANAT / 'code') not in sys.path:
        sys.path.insert(0, str(PLANAT / 'code'))
    return __import__(name), got


def build_sources(ours):
    for pi, n in enumerate(NAMES):
        b = read_bytes(EV4 / n)
        require(hashlib.sha256(b).hexdigest() == EVIDENCE_SHA256[pi], 'dbn4 evidence source hash ' + n)
        e = load_json(EV4 / n)
        require(canon(e) == b, 'evidence source is canonical JSON ' + n)
        gains = []
        for i in range(608):
            try:
                h = F(e['plan']['rows'][i]['height'])
            except (KeyError, ValueError, ZeroDivisionError) as err:
                raise SystemExit('bad plan height in %s row %d: %r' % (n, i, err))
            require((i, h) in ours, 'floor record for cell %d height %s' % (i, h))
            new = ours[(i, h)]
            for s in SH:
                k = 'r_' + s.replace('/', '_')
                old = F(e['rows'][i][k])
                if new[s] > old:
                    e['rows'][i][k] = str(new[s])
                    gains.append((i, float(new[s] - old)))
        nb = canon(e)
        write_bytes(SRC / n, nb)
        log('provider', pi, n, ': floors raised', len(gains), 'of 6080; new source sha256', hashlib.sha256(nb).hexdigest())
        for (a_, b_) in ((0, 100), (100, 200), (200, 300), (300, 400), (400, 500), (500, 608)):
            g = [x[1] for x in gains if a_ <= x[0] < b_]
            log('   cells %3d-%3d: raised %4d, mean gain %.5f, max %.5f' % (a_, b_ - 1, len(g), (sum(g) / len(g)) if g else 0.0, max(g) if g else 0.0))


def main():
    signal.alarm(3600)
    log('# dbn5: Planat unchanged producer + auditor with the new certified paired-probe floors (s7_new_floor_run.py)')
    prod, hp = planat_module('produce_classical_portfolio')
    aud, ha = planat_module('audit_classical_portfolio')
    log('producer sha256', hp, '; auditor sha256', ha, '(both equal Planat SHA256SUMS.json entries)')
    FL = load_json(E2 / 'floors_dbn5.json')
    log('floors_dbn5.json sha256', sha(E2 / 'floors_dbn5.json'), '; records', len(FL['floors']))
    try:
        ours = {(int(r['cell']), F(r['h'])): {s: F(v) for s, v in r['floor'].items()} for r in FL['floors']}
    except (KeyError, ValueError, TypeError, ZeroDivisionError) as err:
        raise SystemExit('bad floors_dbn5.json: %r' % (err,))
    SRC.mkdir(parents=True, exist_ok=True)
    HIST.mkdir(parents=True, exist_ok=True)
    build_sources(ours)
    t = time.time()
    res = prod.produce(SRC, HIST, [0, 1], 843, 32, 'classical_dbn5')
    hist = HIST / 'classical_dbn5.json.gz'
    log('\n## producer: U =', res['upper_exact'], '=', res['upper_decimal'], 'steps', res['steps'], '(%.0f s)' % (time.time() - t))
    log('history sha256', sha(hist))
    S = [load_json(SRC / n) for n in NAMES]
    hashes = [aud.sha(SRC / n) for n in NAMES]
    H = load_json(hist)
    t = time.time()
    A = aud.audit(H, S, hashes)
    log('## auditor: PASS =', A['PASS'], 'U =', A['upper_exact'], 'steps', A['steps'], '(%.0f s)' % (time.time() - t))
    require(A['PASS'] and A['upper_exact'] == res['upper_exact'], 'auditor result')
    neg = aud.negative_controls(H, S, hashes)
    for c in neg['cases']:
        log('  negative control rejected: %-30s %s' % (c['case'], c['reason']))
    log('Planat negative controls rejected: %d' % len(neg['cases']))
    kuse = {}
    for rec in H['records']:
        kuse[rec['kernel']] = kuse.get(rec['kernel'], 0) + 1
    log('records', len(H['records']), '; distinct kernels', len(kuse), '; appended (dbn4) kernels used in',
        sum(v for k, v in kuse.items() if k >= 57), 'records')
    try:
        shutil.copy(hist, E2 / 'classical_dbn5.json.gz')
    except OSError as e:
        raise SystemExit('cannot copy history: %s' % e)
    write_bytes(E2 / 'classical_dbn5_AUDIT.json', (json.dumps(A, indent=1, sort_keys=True) + '\n').encode())
    write_bytes(E2 / '07_planat_run_dbn5.txt', ('\n'.join(LINES) + '\n').encode())


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError, IndexError, TypeError, AssertionError, ZeroDivisionError) as err:
        raise SystemExit('s7_new_floor_run failed: %r' % (err,))
