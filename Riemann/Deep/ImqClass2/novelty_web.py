# novelty_web.py: GitHub (iwaokimura/cln-imq-2k tags and latest commits) and HAL (newest "class number" imaginary quadratic).
import json, urllib.request, urllib.parse
UA = {'User-Agent': 'Mozilla/5.0 (novelty scan; riemann repo)'}
def get(u):
    with urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=60) as r:
        return json.load(r)
try:
    tags = get('https://api.github.com/repos/iwaokimura/cln-imq-2k/tags')
    print('GitHub tags:', [t['name'] for t in tags])
    for c in get('https://api.github.com/repos/iwaokimura/cln-imq-2k/commits?per_page=6'):
        print('  commit', c['commit']['committer']['date'], c['sha'][:8], c['commit']['message'].splitlines()[0][:110])
except Exception as e:
    print('GITHUB FAIL', str(e)[:200])
try:
    q = urllib.parse.urlencode({'q': '"class number" imaginary quadratic', 'sort': 'producedDate_tdate desc', 'rows': 20, 'fl': 'halId_s,title_s,producedDate_tdate', 'wt': 'json'})
    d = get('https://api.archives-ouvertes.fr/search/?' + q)
    print('HAL (20 newest for "class number" imaginary quadratic):')
    for x in d['response']['docs']:
        print('  ', x.get('producedDate_tdate', '?')[:10], x.get('halId_s'), (x.get('title_s') or ['?'])[0][:110])
except Exception as e:
    print('HAL FAIL', str(e)[:200])
