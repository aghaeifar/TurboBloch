import ctypes
import pathlib
import numpy as np


libname = '/usr/local/lib/libbloch.so'

def simulate(b1:np.ndarray, gr:np.ndarray, pr:np.ndarray, b0:np.ndarray=None, dt=1e-3, T1=10000., T2=10000., m0:np.ndarray=None):
    '''
    b1  : (Tesla)      nTime x [X,Y,Z]
    gr  : (Tesla/m)    3 x nTime
    pr :  (meter)      3 x [X,Y,Z]   
    b0  : (Tesla)      [X,Y,Z] x 1
    dt  : (second)     nTime x 1 or 1 x 1
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
    # check and adjust dwell-time
    if isinstance(dt, np.ndarray) == False:
        dt = np.array(dt)
    if dt.size == 1:
        dt = np.ones(nTime, dtype=np.float32, order='F') * dt
    if dt.size != nTime:
        raise ValueError(f'dt must be 1D with size 1 or size nTime={nTime}')
    

    T1 = np.ones(b0.shape, dtype=np.float32, order='F') * T1
    T2 = np.ones(b0.shape, dtype=np.float32, order='F') * T2
    # simulator expects column-major layout (a.k.a Fortran-contiguous)
    b1 = np.asfortranarray(b1) if not b1.flags['F_CONTIGUOUS'] else b1
    gr = np.asfortranarray(gr) if not gr.flags['F_CONTIGUOUS'] else gr
    dt = np.asfortranarray(dt) if not dt.flags['F_CONTIGUOUS'] else dt
    pr = np.asfortranarray(pr) if not pr.flags['F_CONTIGUOUS'] else pr
    b0 = np.asfortranarray(b0) if not b0.flags['F_CONTIGUOUS'] else b0
    m0 = np.asfortranarray(m0) if not m0.flags['F_CONTIGUOUS'] else m0
    # simulator expects single floating point input
    b1 = b1.astype(np.complex64) if b1.dtype != np.complex64 else b1
    gr = gr.astype(np.float32) if gr.dtype != np.float32 else gr
    dt = dt.astype(np.float32) if dt.dtype != np.float32 else dt
    pr = pr.astype(np.float32) if pr.dtype != np.float32 else pr
    b0 = b0.astype(np.float32) if b0.dtype != np.float32 else b0
    m0 = m0.astype(np.float32) if m0.dtype != np.float32 else m0
    # output 
    result = np.zeros_like(m0, dtype=np.float32, order='F')
    # load shared library and simulate
    handle  = ctypes.CDLL(libname)
    handle.bloch_sim.argtypes = [np.ctypeslib.ndpointer(np.complex64, ndim=b1.ndim, flags='F'), # B1
                                 np.ctypeslib.ndpointer(np.float32, ndim=gr.ndim, flags='F'),   # gradients
                                 np.ctypeslib.ndpointer(np.float32, ndim=dt.ndim, flags='F'),   # Dwell-time
                                 np.ctypeslib.ndpointer(np.float32, ndim=b0.ndim, flags='F'),   # off-resonance
                                 np.ctypeslib.ndpointer(np.float32, ndim=pr.ndim, flags='F'),   # spatial pritions
                                 np.ctypeslib.ndpointer(np.float32, ndim=T1.ndim, flags='F'),   # T1
                                 np.ctypeslib.ndpointer(np.float32, ndim=T2.ndim, flags='F'),   # T2
                                 np.ctypeslib.ndpointer(np.float32, ndim=m0.ndim, flags='F'),   # initial magnetization
                                 ctypes.c_int32,                                                # number of spatial pritions
                                 ctypes.c_int32,                                                # number of time points                                 
                                 ctypes.c_bool,                                                 # save all time points
                                 ctypes.c_bool,                                                 # is constant T1 and T2
                                 ctypes.c_bool,                                                 # is constant dwell-time
                                 np.ctypeslib.ndpointer(np.float32, ndim=result.ndim, flags='F') # output
                                ]                                                 

    handle.bloch_sim(b1, gr, dt, b0, pr, T1, T2, m0, nPos, nTime, False, False, False, result)
    return result
