
import re, sys, collections
pats = {k: re.compile(v, re.I) for k, v in {
 'stepanov': r'stepanov', 'multiplic': r'multiplicit', 'iterLift': r'iterat\w*lift|lift\w*iterat|liftchain|chainlift|lift_chain|chain_lift|liftIter|iterLift|nfoldLift|mFoldLift|mfold', 'seesaw': r'seesaw',
 'besselPeriod': r'besselperiod|bessel_period', 'ggp': r'ggp|gangrossprasad|gan_gross', 'rallis': r'rallis', 'siegelweil': r'siegelweil|siegel_weil',
 'mass': r'massformula|siegelmass|mass_formula|siegel_mass', 'chowla': r'chowla', 'weilfunctional': r'weilfunctional|guinand|explicitformula|explicit_formula',
 'determinant_aux': r'auxdet|auxiliarydet|wronsk', 'cost': r'liftcost|costof|cost_of|_cost|Cost\b', 'interlac': r'interlac', 'boxplus': r'boxplus|⊞'}.items()}
cnt = {k: collections.Counter() for k in pats}
ex = {k: [] for k in pats}
with open('/Users/adriansager/playground/riemann/.scratch/alpcore/decls.txt', encoding='utf-8', errors='replace') as fh:
    for line in fh:
        path, _, rest = line.partition(':')
        for k, rx in pats.items():
            if rx.search(rest):
                d = '/'.join(path.split('/')[1:4]); cnt[k][d] += 1
                if len(ex[k]) < 4: ex[k].append(line.strip()[:200])
for k in pats:
    tot = sum(cnt[k].values())
    print('==', k, tot, '|', '; '.join(f'{d}:{n}' for d, n in cnt[k].most_common(5)))
    for e in ex[k][:3]: print('   ', e)
