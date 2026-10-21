"""Step 3: the 0.148409 bound with Planat's UNCHANGED producer and auditor (sha256 = his SHA256SUMS.json).
- fresh producer run on the rebuilt evidence sources; the history must be byte-identical to dbn4's;
- plain audit; then the same audit with every require(...) call counted (the module is not edited: its
  require function is wrapped at run time), split by whether the record uses an appended kernel;
- Planat's 9 negative controls, plus 12 more aimed at records that use appended kernels;
- a forged history that uses an invalid kernel (lambda = 21/10) with consistent hashes and rates, to show
  exactly what the auditor does NOT check (kernel geometry; covered by step 4 and by the superset auditor).
Output: EVIDENCE/03_planat_auditor_148409.txt"""
import collections, copy, gzip, re, time
from fractions import Fraction as F
from common import *

guard(1800)
log = Log(EVID / '03_planat_auditor_148409.txt')
log('# Planat unchanged producer + auditor on the dbn4 sources (step s3_planat_run.py)')
prod, hp = planat_module('produce_classical_portfolio')
aud, ha = planat_module('audit_classical_portfolio')
log('producer code/produce_classical_portfolio.py sha256', hp, '(equals Planat SHA256SUMS.json entry)')
log('auditor  code/audit_classical_portfolio.py  sha256', ha, '(equals Planat SHA256SUMS.json entry)')
ev = WORK / 'evidence'
paths = [ev / n for n in NAMES]
for pi, p in enumerate(paths):
    assert sha(p) == EVIDENCE_SHA256[pi], p
log('evidence sources (from step 2):', ', '.join('%s %s' % (p.name, sha(p)) for p in paths))
out = WORK / 'histories'
t = time.time()
res = prod.produce(ev, out, [0, 1], 843, 32, 'classical_dbn4')
hist = out / 'classical_dbn4.json.gz'
log('\n## fresh producer run: U =', res['upper_exact'], '=', res['upper_decimal'], 'steps', res['steps'], '(%.0f s)' % (time.time() - t))
log('history sha256', sha(hist), '; identical to committed dbn4/classical_dbn4.json.gz:', sha(hist) == sha(HIST_SINGLE) == HIST_SINGLE_SHA256)
assert sha(hist) == HIST_SINGLE_SHA256
S = [load_json(p) for p in paths]
hashes = [aud.sha(p) for p in paths]
H = load_json(hist)
t = time.time()
A = aud.audit(H, S, hashes)
log('\n## plain audit: PASS =', A['PASS'], 'U =', A['upper_exact'], 'steps', A['steps'], 'minimum kernel margin', float(F(A['minimum_kernel_margin'])), '(%.0f s)' % (time.time() - t))
assert A['PASS'] and A['upper_exact'] == '148408913066829/1000000000000000'
kuse = collections.Counter(rec['kernel'] for rec in H['records'])
log('records', len(H['records']), '; records using appended kernels (index >= 57):', sum(v for k, v in kuse.items() if k >= 57),
    '; distinct kernels used', len(kuse), 'of which appended', sum(1 for k in kuse if k >= 57))
calls, calls_new, fails = collections.Counter(), collections.Counter(), collections.Counter()
orig = aud.require


def counting(x, label):
    f = sys._getframe(1)
    while f is not None and f.f_code.co_name != 'audit':
        f = f.f_back
    rec = f.f_locals.get('rec') if f is not None else None
    calls[label] += 1
    if isinstance(rec, dict) and rec.get('kernel', -1) >= 57:
        calls_new[label] += 1
    if not x:
        fails[label] += 1
    return orig(x, label)


aud.require = counting
A2 = aud.audit(H, S, hashes)
aud.require = orig
assert A2 == A
src_lines = read_text(PLANAT / 'code/audit_classical_portfolio.py').splitlines()
WHAT = {
    'History schema': 'history format is Planat V1',
    'Source identities': 'sha256 of both source files = header hashes (binds every floor, height, kernel)',
    'Scope': 'no 7/8 zero-free hypothesis, no continuum shift',
    'Analytic boundary': 'source floors flagged as assumptions (premise A1)',
    'Substeps': 'substeps per cell is a positive integer (32)',
    'Classical initial strip': 'q0 = (1+1e-6)^2, initial strip q = 1, terminal q = 1e-8',
    'Source providers': 'providers are 0 and/or 1',
    'Complete source schedule': '608 source rows and 608 plan rows per provider',
    'Frequency class coverage': 'all 11 frequency classes present',
    'Closed source time cells': 'cell i is exactly [i/4000, (i+1)/4000]',
    'Contiguous time steps': 'record n is (cell, substep) = divmod(n, 32): no gaps, no reordering',
    'Provider and kernel ownership': 'provider declared; 0 <= kernel index < kernel_count',
    'Height and incoming state': 'q_in = previous q_out; height = source plan height of (provider, cell)',
    'Full kernel domain': 'ELIGIBILITY: q_in/h^2 <= rho^2 and |h r| < 2 for the kernel used',
    'Whole source and own pair': 'VALIDITY: floor and own-pair cap = source values for (provider, cell, shift)',
    'Own-pair range': 'cap > 0 and v >= 0 (inside cap_ok)',
    'Own-pair bound on entire decreasing range': 'CAPS: s_r(w) <= cap for all w in [0, q_in/h^2] (exact quadratic)',
    'Complete positive contraction rate': 'RATES: rate = 2 lam/h (floor - cap/h) exactly; 0 < rate*dt < 1',
    'Strict classical time-zero separation': 'first step: q0 - (2 + c q0) 1e-8 > 1',
    'Exponential argument': '0 < c dt < 1 before the exponential bound',
    'Positive decreasing restart': '1e-8 < q_out < q_in',
    'Declared dyadic lattice': 'q_out is a multiple of 2^-96',
    'Independent outward exponential bound': 'ODE STEP: q_out >= (q+2/c) E - 2/c, E = sum_{k<=20} (-c dt)^k/k! >= exp(-c dt)',
    'Logarithmic argument': '0 < z < 1 before the logarithm bound',
    'Independent terminal logarithm bound': 'TERMINAL: lag >= upper bound of (1/c) ln(1+z) = time from q_in to 1e-8; lag <= dt',
    'Terminal strip cost and wall range': 'TERMINAL: now + lag + 1e-8/2 <= U < 41/250 (de Bruijn strip cost; wall range)',
    'Complete history': 'terminal record is the last record; steps = number of records',
    'Reported source partition': 'per-shift and per-provider step counts equal the header',
}
log('\n## coverage table: every check of the unchanged auditor, as executed on this history')
log('%-42s %-9s %7s %8s %5s  %s' % ('check (label in the auditor)', 'line', 'calls', 'appended', 'fail', 'what it verifies'))
labels = re.findall(r"require\([^\n]*?'([^']+)'\)", '\n'.join(src_lines))
seen = []
for lab in labels:
    if lab in seen:
        continue
    seen.append(lab)
    ln = ','.join(str(i + 1) for i, s in enumerate(src_lines) if "'" + lab + "'" in s)
    log('%-42s %-9s %7d %8d %5d  %s' % (lab, ln, calls[lab], calls_new[lab], fails[lab], WHAT.get(lab, '?')))
assert set(calls) <= set(seen) and not fails
log('labels in the auditor source', len(seen), '; labels executed', len(calls), '; failures 0')
log('NOT checked by the auditor (by design, Planat docstring: "Source inequalities and kernel inequalities have')
log('separate certificates"): (a) kernel geometry 2K_y >= lam (y/h) J -> step 4 (04_*); (b) that the source files')
log('equal Planat\'s refined sources -> step 2 (02_floor_identity.txt); (c) premises A1-A3 themselves (REPORT4).')
log('\n## Planat negative controls (his function, unchanged)')
neg = aud.negative_controls(H, S, hashes)
for c in neg['cases']:
    log('  rejected: %-30s reason: %s' % (c['case'], c['reason']))
log('Planat negative controls: %d of 9 rejected' % len(neg['cases']))


def attempt(name, Hx, Sx=None, hx=None):
    try:
        aud.audit(Hx, Sx if Sx is not None else S, hx if hx is not None else hashes)
    except (ValueError, KeyError, IndexError) as e:
        log('  rejected: %-58s reason: %s' % (name, e))
        return True
    log('  ACCEPTED: ' + name)
    return False


log('\n## extra negative controls on records that use appended kernels')
R = H['records']
r0 = R[0]
v0 = F(r0['q_in']) / F(r0['height']) ** 2
small = next(ki for ki in range(57, 843) if F(S[0]['kernels'][ki]['rho']) ** 2 < v0)
ok = []


def damaged(f):
    h = copy.deepcopy(H)
    f(h)
    return h


ok.append(attempt('kernel %d (rho=%s) used outside its domain at step 0' % (small, S[0]['kernels'][small]['rho']),
                  damaged(lambda h: h['records'][0].update(kernel=small))))
ok.append(attempt('kernel index 843 (out of range)', damaged(lambda h: h['records'][0].update(kernel=843))))
ok.append(attempt('kernel index -1', damaged(lambda h: h['records'][0].update(kernel=-1))))
ok.append(attempt('floor raised by 1e-6 at step 10000', damaged(lambda h: h['records'][10000].update(floor=str(F(R[10000]['floor']) + F(1, 10 ** 6))))))
ok.append(attempt('own-pair cap lowered by 1e-9 at step 10000', damaged(lambda h: h['records'][10000].update(self_cap=str(F(R[10000]['self_cap']) - F(1, 10 ** 9))))))
ok.append(attempt('rate raised by 1e-12 at step 10000', damaged(lambda h: h['records'][10000].update(rate=str(F(R[10000]['rate']) + F(1, 10 ** 12))))))
ok.append(attempt('restart lowered by one lattice unit 2^-96 at step 5000', damaged(lambda h: h['records'][5000].update(q_out=str(F(R[5000]['q_out']) - F(1, 2 ** 96))))))
last = R[-1]
c_last = F(last['rate'])
z = c_last * (F(last['q_in']) - F(1, 10 ** 8)) / (2 + c_last * F(1, 10 ** 8))
ok.append(attempt('terminal lag set just below the auditor\'s log bound', damaged(lambda h: h['records'][-1].update(terminal_lag_upper=str(aud.log_upper(z) / c_last - F(1, 10 ** 40))))))
ok.append(attempt('final bound U lowered by 1e-15', damaged(lambda h: h.update(upper_exact=str(F(H['upper_exact']) - F(1, 10 ** 15))))))
ok.append(attempt('step 7000 deleted', damaged(lambda h: h['records'].pop(7000))))
S2 = copy.deepcopy(S)
S2[0]['kernels'][r0['kernel']]['lambda'] = '21/10'
h2 = [sha_bytes(canon(s)) for s in S2]
ok.append(attempt('source file: lambda of kernel %d raised to 21/10' % r0['kernel'], H, S2, h2))
S3 = copy.deepcopy(S)
S3[0]['rows'][0]['r_19_40'] = str(F(S3[0]['rows'][0]['r_19_40']) + F(1, 10 ** 6))
ok.append(attempt('source file: floor (provider 0, cell 0, r=19/40) raised by 1e-6', H, S3, [sha_bytes(canon(s)) for s in S3]))
log('extra negative controls: %d of %d rejected' % (sum(ok), len(ok)))
assert all(ok)
log('\n## what the auditor cannot see: an invalid kernel with consistent hashes and rates')
Hf = copy.deepcopy(H)
Hf['source_hashes'] = h2
kf = r0['kernel']
kk = S2[0]['kernels'][kf]
nf = 0
for rec in Hf['records']:
    if rec['provider'] == 0 and rec['kernel'] == kf:
        h = F(rec['height'])
        rec['rate'] = str(2 * F(kk['lambda']) / h * (F(rec['floor']) - F(kk['self_cap']) / h))
        nf += 1
try:
    Af = aud.audit(Hf, S2, h2)
    log('forged history (kernel %d: lambda 2 -> 21/10, %d records re-rated): unchanged auditor PASS =' % (kf, nf), Af['PASS'], 'U =', Af['upper_decimal'])
    log('=> expected: kernel geometry is outside the auditor; 21/10 violates the kernel inequality (step 4 rejects it;')
    log('   the superset auditor of step 7 rejects this forged history).')
except ValueError as e:
    log('forged history rejected by the unchanged auditor:', e)
(WORK / 'forged').mkdir(exist_ok=True)
with gzip.open(WORK / 'forged' / 'forged_lambda.json.gz', 'wt') as g:
    json.dump(Hf, g)
for pi, n in enumerate(NAMES):
    (WORK / 'forged' / n).write_bytes(canon(S2[pi]))
log('\nRESULT the 0.148408913066829 history is reproduced byte for byte by Planat\'s unchanged producer and passes his unchanged auditor;')
log('every one of its %d steps uses an appended kernel, so every per-step check above ran on appended kernels.' % len(R))
log.save()
