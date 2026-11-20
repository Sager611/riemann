import tarfile, gzip, os, io, sys
for id in open('ids.txt').read().split():
    d=f'src/{id}'; os.makedirs(d,exist_ok=True)
    raw=open(f'{id}.tgz','rb').read()
    try:
        with tarfile.open(fileobj=io.BytesIO(raw),mode='r:gz') as t:
            t.extractall(d, filter='data')
    except tarfile.ReadError:
        open(f'{d}/main.tex','wb').write(gzip.decompress(raw))
    tex=[os.path.join(r,f) for r,_,fs in os.walk(d) for f in fs if f.endswith('.tex')]
    print(id, [ (os.path.relpath(t,d), sum(1 for _ in open(t,errors='replace'))) for t in tex])
