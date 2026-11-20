"""Constant-memory row reader for the string arrays of Weisz's zeta_zeros_14400.npz (no full np.load)."""
import zipfile
import numpy as np

def _hdr(f, ver):
    rd = np.lib.format.read_array_header_1_0 if ver == (1, 0) else np.lib.format.read_array_header_2_0
    return rd(f, max_header_size=10**6) if 'max_header_size' in rd.__code__.co_varnames else rd(f)

class NpzRows:
    """Iterate rows of one 1-D array inside an .npz without loading the array."""
    def __init__(self, npz_path, name):
        self.zf = zipfile.ZipFile(npz_path)
        self.f = self.zf.open(name + '.npy')
        ver = np.lib.format.read_magic(self.f)
        shape, fortran, dtype = _hdr(self.f, ver)
        assert len(shape) == 1 and not fortran, (name, shape, fortran)
        self.n, self.dtype = shape[0], dtype

    def __iter__(self):
        it = self.dtype.itemsize
        for _ in range(self.n):
            buf = self.f.read(it)
            assert len(buf) == it
            v = np.frombuffer(buf, dtype=self.dtype)[0]
            yield v.item() if hasattr(v, 'item') else v

def members(npz_path):
    with zipfile.ZipFile(npz_path) as zf:
        out = {}
        for nm in zf.namelist():
            f = zf.open(nm)
            ver = np.lib.format.read_magic(f)
            out[nm] = _hdr(f, ver)
        return out
