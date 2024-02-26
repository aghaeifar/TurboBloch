import ctypes
import pathlib
import numpy as np


libname = pathlib.Path().absolute() / ".." / "lib" / "libbloch.so"

def bloch_simulation(b1:np.ndarray, gr:np.ndarray, pr:np.ndarray, b0:np.ndarray=None, td=1e-3, T1=10000., T2=10000., m0:np.ndarray=None):
    '''
    b1  : (Tesla)      nTime x [X,Y,Z]
    gr  : (Tesla/m)    3 x nTime
    pr :  (meter)      3 x [X,Y,Z]   
    b0  : (Tesla)      [X,Y,Z] x 1
    td  : (second)    
    T1  : (second)
    T2  : (second)
    m0  : ()           3 x [X,Y,Z]
    '''

    if pr.ndim==1: # must be at least 2D
        pr = np.expand_dims(pr, -1)
    if b1.ndim==1: # if is 1D, take spatial dimensions from pr
        for i in range(1, pr.ndim):
            b1 = np.expand_dims(b1, -1).repeat(pr.shape[i], -1)
    pr_shape = b1.shape[1:]     # spatial dimension shape
    nTime = b1.shape[0]         # number of time points
    nPos  = np.prod(pr_shape)   # number of spatial points
    
    if gr.shape != (3, nTime) or pr.shape != (3,)+pr_shape:
        raise ValueError(f'gr shape {gr.shape} or pr shape {[pr.shape]} does not match b1 shape {b1.shape}.')

    if b0 is None:
        b0 = np.zeros(pr_shape, dtype=np.float32, order='F')    
    if m0 is None:
        m0 = np.zeros((3,)+pr_shape, dtype=np.float32, order='F')
        m0[-1,:] = 1

    if b0.shape != pr_shape or m0.shape != (3,)+pr_shape:
        raise ValueError(f'Error in b0 shape {b0.shape} or m0 shape {[m0.shape]}. Spatial dimension size must be {pr_shape}')
    
    T1 = np.ones(b0.shape, dtype=np.float32, order='F') * T1
    T2 = np.ones(b0.shape, dtype=np.float32, order='F') * T2
    # simulator expects column-major layout (a.k.a Fortran-contiguous)
    b1 = np.asfortranarray(b1) if not b1.flags['F_CONTIGUOUS'] else b1
    gr = np.asfortranarray(gr) if not gr.flags['F_CONTIGUOUS'] else gr
    pr = np.asfortranarray(pr) if not pr.flags['F_CONTIGUOUS'] else pr
    b0 = np.asfortranarray(b0) if not b0.flags['F_CONTIGUOUS'] else b0
    m0 = np.asfortranarray(m0) if not m0.flags['F_CONTIGUOUS'] else m0
    # simulator expects single floating point input
    b1 = b1.astype(np.complex64) if b1.dtype != np.complex64 else b1
    gr = gr.astype(np.float32) if gr.dtype != np.float32 else gr
    pr = pr.astype(np.float32) if pr.dtype != np.float32 else pr
    b0 = b0.astype(np.float32) if b0.dtype != np.float32 else b0
    m0 = m0.astype(np.float32) if m0.dtype != np.float32 else m0
    # output 
    result = np.zeros_like(m0, dtype=np.float32, order='F')
    # load shared library and simulate
    handle  = ctypes.CDLL(libname)
    handle.bloch_sim.argtypes = [np.ctypeslib.ndpointer(np.complex64, ndim=b1.ndim, flags='F'), # B1
                                 np.ctypeslib.ndpointer(np.float32, ndim=gr.ndim, flags='F'),   # gradients
                                 ctypes.c_float,                                                # Dwell-time
                                 np.ctypeslib.ndpointer(np.float32, ndim=b0.ndim, flags='F'),   # off-resonance
                                 np.ctypeslib.ndpointer(np.float32, ndim=pr.ndim, flags='F'),   # spatial pritions
                                 np.ctypeslib.ndpointer(np.float32, ndim=T1.ndim, flags='F'),   # T1
                                 np.ctypeslib.ndpointer(np.float32, ndim=T2.ndim, flags='F'),   # T2
                                 np.ctypeslib.ndpointer(np.float32, ndim=m0.ndim, flags='F'),   # initial magnetization
                                 ctypes.c_int32,                                                # number of spatial pritions
                                 ctypes.c_int32,                                                # number of time points
                                 np.ctypeslib.ndpointer(np.float32, ndim=result.ndim, flags='F'), # output
                                 ctypes.c_bool,                                                 # save all time points
                                 ctypes.c_bool]                                                 # is constant T1 and T2

    handle.bloch_sim(b1, gr, td, b0, pr, T1, T2, m0, nPos, nTime, result, False, False)
    return result
