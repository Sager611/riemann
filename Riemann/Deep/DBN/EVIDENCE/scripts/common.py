"""Shared paths and helpers for the EVIDENCE scripts (standard library only).

Layout (all paths derived from this file's location):
  REPO/Riemann/Deep/DBN/EVIDENCE/scripts/   these scripts
  REPO/Riemann/Deep/DBN/dbn4/               dbn4 catalogue, mixtures, histories, producer code
  WORK = $DBN4V_WORK or REPO/.scratch/dbn4v  scratch (arXiv download + extraction, rebuilt inputs, new histories)
  PLANAT = WORK/arxiv/anc/certificate        Planat's certificate folder, extracted from the verified arXiv tarball
"""
import hashlib, json, os, signal, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EVID = HERE.parent
DBN = EVID.parent
REPO = DBN.parents[2]
DBN4 = DBN / 'dbn4'
WORK = Path(os.environ.get('DBN4V_WORK', str(REPO / '.scratch' / 'dbn4v')))
PLANAT = WORK / 'arxiv' / 'anc' / 'certificate'
TARBALL = WORK / 'dl' / 'e-print_2610.18214v1'
TARBALL_URL = 'https://arxiv.org/e-print/2610.18214v1'
TARBALL_SHA256 = '4599b663fe3c0913b78c445834b93ff2721e33e819601185a43bb88dd226ed5b'
NAMES = ['SOURCE_OLD_ENHANCED.json', 'SOURCE_RETUNED.json']
# hashes written into the dbn4 history headers (the evidence files Planat's auditor binds to)
EVIDENCE_SHA256 = ['235c34ae1c0efb03458a39687685451f9b2e0febbaab3d2d6ef38cf9d9cad714',
                   '29840261b1fee53867a61ab32c4ab1f4b81ecc588c93463b0ef28d02fcb76058']
HIST_SINGLE = DBN4 / 'classical_dbn4.json.gz'
HIST_MIX = DBN4 / 'classical_dbn4_mix.json.gz'
HIST_SINGLE_SHA256 = 'e8b13ddcc455a5ed18cc60d4a022f0d63b249f19af12d7cbb50e8e421e5ec2fc'


def read_bytes(p):
    try:
        return Path(p).read_bytes()
    except OSError as e:
        raise SystemExit('cannot read %s: %s' % (p, e))


def read_text(p):
    return read_bytes(p).decode()


def loads(s):
    try:
        return json.loads(s)
    except ValueError as e:
        raise SystemExit('bad JSON: %s' % e)


def sha(p):
    return hashlib.sha256(read_bytes(p)).hexdigest()


def sha_bytes(b):
    return hashlib.sha256(b).hexdigest()


def canon(x):
    """The serialisation Planat's sources and the dbn4 evidence files use."""
    return (json.dumps(x, separators=(',', ':'), sort_keys=True) + '\n').encode()


def guard(seconds):
    """Hard wall-clock limit: SIGALRM terminates the process."""
    signal.alarm(int(seconds))


class Log:
    def __init__(self, path):
        self.path = Path(path)
        self.lines = []

    def __call__(self, *a):
        s = ' '.join(str(x) for x in a)
        print(s, flush=True)
        self.lines.append(s)

    def save(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text('\n'.join(self.lines) + '\n')


def load_json(p):
    """Read a JSON (or .json.gz) file; fail closed with a clear message."""
    import gzip
    try:
        if str(p).endswith('.gz'):
            with gzip.open(p, 'rt') as f:
                return json.load(f)
        return json.loads(Path(p).read_text())
    except (OSError, ValueError) as e:
        raise SystemExit('cannot read %s: %s' % (p, e))


def planat_sums():
    return load_json(PLANAT / 'SHA256SUMS.json')


def planat_module(name):
    """Import one of Planat's code/*.py modules after checking its bytes against his SHA256SUMS.json."""
    rel = 'code/' + name + '.py'
    want = planat_sums()[rel]['sha256']
    got = sha(PLANAT / rel)
    if got != want:
        raise SystemExit('Planat module changed: %s %s != %s' % (rel, got, want))
    sys.dont_write_bytecode = True
    if str(PLANAT / 'code') not in sys.path:
        sys.path.insert(0, str(PLANAT / 'code'))
    mod = __import__(name)
    return mod, got
