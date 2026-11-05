"""fxkern.py -- ctypes wrapper of libfxkern (dense midpoint kernels on fixed-point blocks; see fxkern.c)."""
import ctypes, os
import numpy as np

_lib = None
def lib(path=None):
    global _lib
    if _lib is None:
        path = path or os.environ.get('FXKERN_LIB', os.path.join(os.getcwd(), 'libfxkern.dylib'))
        _lib = ctypes.CDLL(path)
        u16p = np.ctypeslib.ndpointer(dtype=np.uint16, flags='C_CONTIGUOUS')
        L, I = ctypes.c_long, ctypes.c_int
        _lib.fx_chol.argtypes = [L, I, L, L, u16p, u16p, L]; _lib.fx_chol.restype = I
        _lib.fx_trsm.argtypes = [L, L, I, L, L, u16p, u16p, u16p, L]; _lib.fx_trsm.restype = I
        _lib.fx_trsm_left.argtypes = [L, L, I, L, L, u16p, u16p, u16p, L]; _lib.fx_trsm_left.restype = I
        _lib.fx_trsm_left_t.argtypes = [L, L, I, L, L, u16p, u16p, u16p, L]; _lib.fx_trsm_left_t.restype = I
    return _lib

def chol(fx, A, prec):
    n = A.shape[0]; A = np.ascontiguousarray(A); out = np.zeros_like(A)
    r = lib().fx_chol(n, fx.ND, fx.F, fx.T, A, out, prec)
    return out, r

def trsm(fx, S, Ld, prec):
    """X with X L^T = S."""
    m, n = S.shape[0], S.shape[1]; S = np.ascontiguousarray(S); Ld = np.ascontiguousarray(Ld); out = np.zeros_like(S)
    r = lib().fx_trsm(m, n, fx.ND, fx.F, fx.T, S, Ld, out, prec)
    if r: raise OverflowError("fx_trsm: value out of range")
    return out

def trsm_left(fx, Ld, B, prec, trans=False):
    """X with L X = B (trans: L^T X = B)."""
    n, m = B.shape[0], B.shape[1]; B = np.ascontiguousarray(B); Ld = np.ascontiguousarray(Ld); out = np.zeros_like(B)
    f = lib().fx_trsm_left_t if trans else lib().fx_trsm_left
    r = f(n, m, fx.ND, fx.F, fx.T, Ld, B, out, prec)
    if r: raise OverflowError("fx_trsm_left: value out of range")
    return out
