import json
print("| certificate | θ (7/8) | C′ all x (7/8) | θ (B′) | C′ all x (B′) | C, x ≥ x_a (B′) |")
print("|---|---|---|---|---|---|")
for r in json.load(open('logs/bprime_results.json')):
    print("| %s | %.4f | 10^%.2f | %.4f | 10^%.2f | 10^%.2f (x ≥ 10^%d) |" % (r['name'].replace('_Bp', ''), r['old_theta'], r['old_log10Cprime'],
          r['theta'], r['allx']['Cprime_log10'], r['best']['log10C'], r['best']['log10xa']))
for r in json.load(open('logs/bprime_results.json')):
    if 'C1C2_recheck' in r: print("recheck", r['name'], r['C1C2_recheck'])
