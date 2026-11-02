
import json, lzma, gzip, statistics as st, collections, os, sys
cert = sys.argv[1]
man = json.load(open(cert + '/manifest.json'))
nk, nd, tot_digits, maxd = [], [], 0, collections.Counter()
raw_nums = raw_keys = 0
for name in man['files']:
    with lzma.open(cert + '/' + name, 'rt') as f:
        for line in f:
            idx, den, dk, vals = json.loads(line)
            nk.append(len(vals))
            ds = [len(str(abs(v))) for v in vals]
            nd.extend(ds); maxd[max(ds)] += 1
            raw_nums += len(','.join(map(str, vals)))
            raw_keys += len(','.join(map(str, dk)))
size = sum(os.path.getsize(cert + '/' + n) for n in man['files'])
L = len(nk)
print('leaves', L, 'mult bytes', size, 'B/leaf %.0f' % (size / L))
print('keys/leaf mean %.1f median %d max %d' % (st.mean(nk), st.median(nk), max(nk)))
print('digits/num mean %.2f median %d' % (st.mean(nd), st.median(nd)))
print('max digits per leaf hist', sorted(maxd.items()))
print('raw chars: nums %.0f/leaf, keys %.0f/leaf' % (raw_nums / L, raw_keys / L))
# entropy estimate of numerators if digits random: sum(digits)*log2(10)/8
print('entropy-ish numerator bytes/leaf %.0f' % (sum(nd) * 3.3219 / 8 / L))
