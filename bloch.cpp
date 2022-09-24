/* Bloch simulation
 *
 * Author : Ali Aghaeifar <ali.aghaeifar.mri [at] gmail [dot] com>
 *
 */
#include "bloch.h"
#include <iostream>
#include <vector>
#include <execution>
#include <numeric> // std::inner_product, std::iota 
#include <algorithm>

#ifndef __NOPTX__
#include "mkl.h"
#ifdef __SINGLE_PRECISION__
typedef MKL_Complex8  _MKL_COMPLEX;
#else
typedef MKL_Complex16 _MKL_COMPLEX;
#endif
#endif


#ifdef __SEQUENTIAL__
#define __MODE__  (std::execution::seq)
#else
#define __MODE__  (std::execution::par_unseq)
#endif



#define GAMMA_T 267522187.44

#ifndef __NOPTX__
void (*p_cblas_Xgemm)(const CBLAS_LAYOUT, const CBLAS_TRANSPOSE, const CBLAS_TRANSPOSE, const MKL_INT, const MKL_INT, const MKL_INT, const void *, const void *, const MKL_INT, const void *, const MKL_INT, const void *, void *, const MKL_INT);
#endif

// q = {x, y, z, w}
void apply_rot_quaternion(_T q[4], _T *m0, _T *m1)
{
    _T t[3];
    // see https://blog.molecular-matters.com/2013/05/24/a-faster-quaternion-vector-multiplication/
    if(q[0] == 0 && q[1] == 0) // no RF pulse case
    {
        t[0] =-2*q[2]*m0[1];
        t[1] = 2*q[2]*m0[0];

        m1[0] = m0[0] + q[3]*t[0] - q[2]*t[1];
        m1[1] = m0[1] + q[3]*t[1] + q[2]*t[0];
        m1[2] = m0[2];
    }
    else if(q[2] == 0) // no gradient and off-resonance
    {
        t[0] = 2*q[1]*m0[2];
        t[1] =-2*q[0]*m0[2];
        t[2] = 2*(q[0]*m0[1] - q[1]*m0[0]);

        m1[0] = m0[0] + q[3]*t[0] + q[1]*t[2];
        m1[1] = m0[1] + q[3]*t[1] - q[0]*t[2];
        m1[2] = m0[2] + q[3]*t[2] + q[0]*t[1] - q[1]*t[0];
    }
    else
    {
        t[0] = 2*(q[1]*m0[2] - q[2]*m0[1]);
        t[1] = 2*(q[2]*m0[0] - q[0]*m0[2]);
        t[2] = 2*(q[0]*m0[1] - q[1]*m0[0]);

        m1[0] = m0[0] + q[3]*t[0] + q[1]*t[2] - q[2]*t[1];
        m1[1] = m0[1] + q[3]*t[1] + q[2]*t[0] - q[0]*t[2];
        m1[2] = m0[2] + q[3]*t[2] + q[0]*t[1] - q[1]*t[0];
    }
}

void create_quaternion(_T nx, _T ny, _T nz, _T q[4])
{ 
    // see https://paroj.github.io/gltut/Positioning/Tut08%20Quaternions.html
    _T phi = 0., s = 0., c = 0.;
    int n = 1;
    // q = sin(θ/2)(xi+yj+zk) + cos(θ/2)
    if(nx == 0 && ny == 0)
    {  // Only gradients and off-resonance, no RF
        phi  = abs(nz);
        s    = sin(0.5 * phi);
        q[0] = 0.;
        q[1] = 0.;
        q[2] = s * (nz>0 ? 1:-1); // equal to nz/phi == nz/abs(nz) == sign(nz)
    }
    else if(nz == 0)
    {  // Only RF, no gradients and off-resonance
        phi  = sqrt(nx*nx + ny*ny);
        s    = sin(0.5 * phi);
        _T sp= s / phi;
        q[0] = nx * sp;
        q[1] = ny * sp;
        q[2] = 0.; 
    }
    else
    {
        phi  = sqrt(nx*nx + ny*ny + nz*nz);
        s    = sin(0.5 * phi);
        _T sp= s / phi; // /phi because [nx, ny, nz] is unit length in definition. This will be effective in the next lines, where nx, ny, nz * sp 
        q[0] = nx * sp;
        q[1] = ny * sp;
        q[2] = nz * sp;        
    }
    q[3] = sqrt(1. - s*s);// cos(0.5 * phi); // sqrt(1. - s*s) is faster
}

// ----------------------------------------------- //

void bloch::timekernel( std::complex<_T> *b1xy, 
                        _T *gr,
                        _T *pr, 
                        _T b0, 
                        _T td_gamma, 
                        _T *m0,
                        _T e1, _T e2, 
                        _T *output)
{
    _T m1[3], q[4];
    _T rotx, roty, rotz;
    _T e1_1 = e1 - 1;  
    std::copy(m0, m0+3, output);  // setting starting magnetization
    
    if(e1 > 0 && e2 > 0) // including relaxations
    {        
        for (int ct=0; ct<m_lNTime; ct++, gr += 3)
        {            
            // rotations are right handed, thus all are negated.
            rotx = -b1xy[ct].real() * td_gamma;
            roty = -b1xy[ct].imag() * td_gamma;            
            rotz = -std::inner_product(gr, gr+3, pr, b0) * td_gamma; // -(gx*px + gy*py + gz*pz + b0) * dT * gamma
            
            create_quaternion(-rotx, -roty, -rotz, q); // quaternion needs additional sign reverse because looking down the axis of rotation, positive rotations appears clockwise
            apply_rot_quaternion(q, output, m1);

            output += m_lStepTime;
            output[0] = m1[0] * e2;
            output[1] = m1[1] * e2;
            output[2] = m1[2] * e1 - e1_1;            
        }
    }
    else // excluding relaxations, can we make combine all rotations into one and and improve run time? 
    {   
        /*
        std::complex<_T> al0(1,0), bt0(0,0), al1, bt1, al2, bt2;
        for (int ct=0; ct<m_lNTime; ct++, gr += 3)
        {            
            // rotations are right handed, thus all are negated.
            rotx = -b1xy[ct].real() * td_gamma;
            roty = -b1xy[ct].imag() * td_gamma;            
            rotz = -std::inner_product(gr, gr+3, pr, b0) * td_gamma; // -(gx*px + gy*py + gz*pz + b0) * dT * gamma
           
            create_CayleyKlein(rotx, roty, rotz, al1, bt1);             
                                
            al2 = al1 * al0 - std::conj(bt1) * bt0;
            bt2 = bt1 * al0 + std::conj(al1) * bt0;
            al0 = al2;
            bt0 = bt2; 
                               
        }   
        // consider https://www.intel.com/content/www/us/en/developer/articles/technical/onemkl-improved-small-matrix-performance-using-just-in-time-jit-code.html
        apply_rot_CayleyKlein(al2.real(), al2.imag(), bt2.real(), bt2.imag(), m0, output);  
        */ 
    } 
}

// ----------------------------------------------- //

bloch::bloch(long nPosition, long nTime, long nCoil, bool saveAll, bool isPTx)
{
    m_lNPos   = nPosition;
    m_lNTime  = nTime;
    m_lNCoil  = nCoil;
    m_lStepPos    = saveAll ? 3 * (nTime+1) : 3;
    m_lStepTime   = saveAll ? 3 : 0;
    m_lStepB1     = 0;
    m_pB1combined = NULL;
#ifndef __NOPTX__
    if(isPTx)
    {
        m_pB1combined = new std::complex<_T>[nTime*nPosition];
        m_lStepB1 = nTime;
    }
        
#ifdef __SINGLE_PRECISION__
    p_cblas_Xgemm = &cblas_cgemm;
#else
    p_cblas_Xgemm = &cblas_zgemm;
#endif
#endif
}

bloch::~bloch()
{
    if(m_pB1combined != NULL)
        delete[] m_pB1combined;
}


// ----------------------------------------------- //

bool bloch::run(std::complex<_T> *pB1,   // m_lNTime x m_lNCoil [Volt]: column-major order {t0c0, t1c0, t2c0,...,t0c1, t1c1, t2c1,...}
                _T *pGr,                 // 3 x m_lNTime [Tesla/m] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
                _T td,                   // [second]
                _T *pB0,                 // m_lNPos x 1  [Tesla]
                _T *pPos,                // 3 x m_lNPos  [meter] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
                std::complex<_T> *pSens, // m_lNCoil x m_lNPos [Tesla/Volt]: column-major order {c0p0, c1p0, c2p0,...,c0p1, c1p1, c2p1,...}
                _T T1, _T T2,            // [second]
                _T *pM0,                 // 3 x m_lNPos : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
                _T *pResult)             // 3 x (m_lNTime+1) x m_lNPos : column-major order {x1t0,y1t0,z1t0,...,x1tn,y1tn,z1tn,x2t0,y2t0,z2t0,...}, result equals m0 at t0
{
    // Calculate the E1 and E2 values.
    _T e1 = T1 <= 0. ? -1.0 : exp(-td/T1);
    _T e2 = T2 <= 0. ? -1.0 : exp(-td/T2);
    _T td_gamma = td * GAMMA_T;
    std::complex<_T> *pB1combined = pB1;
 
#ifndef __NOPTX__
    if(m_pB1combined != NULL)
    {
        _MKL_COMPLEX alpha, beta;
        alpha.real = 1.0; alpha.imag = 0.0;
        beta.real  = 0.0; beta.imag  = 0.0;
        // consider gemm3m for faster calculation but higher numerical rounding errors
        p_cblas_Xgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m_lNTime, m_lNPos, m_lNCoil, &alpha, pB1, m_lNTime, pSens, m_lNCoil, &beta, m_pB1combined, m_lNTime);
        pB1combined = m_pB1combined;
    }
#endif 

    // =================== Do The Simulation! =================== 
    try
    {
        std::vector<int> a(m_lNPos);
        std::iota (a.begin(), a.end(),0);
        std::for_each (__MODE__, std::begin(a), std::end(a), [&](int cpos){
            timekernel(pB1combined+cpos*m_lStepB1, pGr, pPos+cpos*3, *(pB0+cpos), td_gamma, pM0+cpos*3, e1, e2, pResult+cpos*m_lStepPos);
        });      
    }
    catch( std::exception &ex )
    {
        std::cout << "Simulation failed." << std::endl;
        std::cout << ex.what() << std::endl;
        return false;
    } 
    return true;
}


extern "C" {
        bool bloch_sim(
        std::complex<_T> *pB1,   // m_lNTime x m_lNCoils [Volt] : column-major order {t0c0, t1c0, t2c0,...,t0c1, t1c1, t2c1,...}
        _T *pGr,                 // 3 x m_lNTime [Tesla/m] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
        _T td,                   // [second]
        _T *pB0,                 // m_lNPos x 1  [Tesla]
        _T *pPos,                // 3 x m_lNPos  [meter] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}                    
        std::complex<_T> *pSens, // m_lNCoils x m_lNPos [Tesla/Volt]: column-major order {c0p0, c1p0, c1p0,...,c0p1, c1p1, c1p1,...}
        _T T1, _T T2,            // [second]
        _T *pM0,                 // 3 x m_lNPos : column-maj
        long nPosition, 
        long nTime, 
        long nCoil,
        _T *pResult,             // 3 x (m_lNTime+1) x m_lNPos : column-major order {x1t0,y1t0,z1t0,...,x1tn,y1tn,z1tn,x2t0,y2t0,z2t0,...}, result equals m0 at t0
        bool saveAll,            // return all time-points or only the final magnetization         
        bool isPTx)
{
    bloch bloch_obj(nPosition, nTime, nCoil, saveAll, isPTx);
    return bloch_obj.run(pB1, pGr, td, pB0, pPos, pSens, T1, T2, pM0, pResult);
}
} // extern