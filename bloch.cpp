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

#ifdef __FASTER__
#include "sincos_table.h" 
#endif


#ifdef __SEQUENTIAL__
#define __MODE__  (std::execution::seq)
#else
#define __MODE__  (std::execution::par_unseq)
#endif
     

 
#define GAMMA_T 267522187.44  

_T (*mysin)(_T in); 
_T (*mycos)(_T in);

// q = {x, y, z, w}
void apply_rot_quaternion(_T q[4], _T *m0, _T *m1)
{
    _T t[3];
    // see https://blog.molecular-matters.com/2013/05/24/a-faster-quaternion-vector-multiplication/
    if(q[0] == 0 && q[1] == 0) 
    {   // Only gradients and off-resonance, no RF
        t[0] =-2*q[2]*m0[1];
        t[1] = 2*q[2]*m0[0];

        m1[0] = m0[0] + q[3]*t[0] - q[2]*t[1];
        m1[1] = m0[1] + q[3]*t[1] + q[2]*t[0];
        m1[2] = m0[2];
    }
    else if(q[2] == 0)
    {   // Only RF, no gradients and off-resonance
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
    //q = sin(θ/2)(xi+yj+zk) + cos(θ/2)
    if(nx == 0 && ny == 0)
    {   // Only gradients and off-resonance, no RF
        phi  = abs(nz); 
        s    = nz == 0? 0:mysin(0.5 * phi);
        q[0] = 0.;
        q[1] = 0.;
        q[2] = s * (nz>0 ? 1:-1); // equal to nz/phi == nz/abs(nz) == sign(nz)
    }
    else if(nz == 0)
    {  // Only RF, no gradients and off-resonance
        phi  = sqrt(nx*nx + ny*ny);
        s    = mysin(0.5 * phi);
        _T sp= s / phi; 
        q[0] = nx * sp;
        q[1] = ny * sp;
        q[2] = 0.; 
    }
    else
    {
        phi  = sqrt(nx*nx + ny*ny + nz*nz);
        s    = mysin(0.5 * phi);
        _T sp= s / phi; // /phi because [nx, ny, nz] is unit length in definition. This will be effective in the next lines, where nx, ny, nz * sp 
        q[0] = nx * sp;
        q[1] = ny * sp;
        q[2] = nz * sp;        
    }
    q[3] = mycos(0.5 * phi); // sqrt(1. - s*s) is faster
}

// ----------------------------------------------- //

void bloch::timekernel( const std::complex<_T> *b1xy, 
                        const _T *gr,
                        const _T *pr, 
                        const _T *b0, 
                        const _T td_gamma, 
                        const _T *m0,
                        const _T e1, const _T e2, 
                        _T *output)
{
    _T m1[3], q[4];
    _T rotx, roty, rotz;
    _T e1_1 = e1 - 1;  
    std::copy(m0, m0+3, output);  // setting starting magnetization

    for (int ct = 0; ct < m_lNTime; ct++, gr += 3)
    {
        // rotations are right handed, thus all are negated.
        rotx = -b1xy[ct].real() * td_gamma;
        roty = -b1xy[ct].imag() * td_gamma;
        rotz = -std::inner_product(gr, gr + 3, pr, *b0) * td_gamma; // -(gx*px + gy*py + gz*pz + b0) * dT * gamma

        create_quaternion(-rotx, -roty, -rotz, q); // quaternion needs additional sign reverse because looking down the axis of rotation, positive rotations appears clockwise
        apply_rot_quaternion(q, output, m1);

        output += m_lStepTime;
        output[0] = m1[0] * e2;
        output[1] = m1[1] * e2;
        output[2] = m1[2] * e1 - e1_1;

        if(m_bHasDiffusion)
        {
            pr += 3;
            b0 += 1;
        }
    }
}

// ----------------------------------------------- //

bloch::bloch(long nPosition, long nTime, bool hasDiffusion, bool saveAll)
{
    m_lNPos   = nPosition;
    m_lNTime  = nTime;
    m_lStepPos   = saveAll ? (nTime+1) : 1;
    m_lStepTime  = saveAll ? 3 : 0;
    m_bHasDiffusion = hasDiffusion;

#ifdef __FASTER__
mysin = &fast_sin;
mycos = &fast_cos;
#else
#ifdef __SINGLE_PRECISION__
mysin = &sinf;
mycos = &cosf;
#else
mysin = &sin;
mycos = &cos;
#endif
#endif
}

bloch::~bloch()
{

}


// ----------------------------------------------- //

bool bloch::run(const std::complex<_T> *pB1,   // RF; m_lNTime x 1 [Volt]
                const _T *pGr,                 // gradients; 3 x m_lNTime [Tesla/m] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
                const _T td,                   // dwell-time; [second]
                const _T *pB0,                 // off-resonance; m_lNPos x 1  [Tesla]
                const _T *pPos,                // spatial positions; 3 x m_lNPos  [meter] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
                const _T T1, const _T T2,      // relaxations; [second]
                const _T *pM0,                 // initial magnetization; 3 x m_lNPos : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
                 _T *pResult)                  // 3 x (m_lNTime+1) x m_lNPos : column-major order {x1t0,y1t0,z1t0,...,x1tn,y1tn,z1tn,x2t0,y2t0,z2t0,...}, result equals m0 at t0
{
    // Calculate the E1 and E2 values.
    _T e1 = T1 <= 0. ? 1.0 : exp(-td/T1);
    _T e2 = T2 <= 0. ? 1.0 : exp(-td/T2);
    _T td_gamma = td * GAMMA_T;
 
    // =================== Do The Simulation! =================== 
    try
    {
        int step_diffusion = m_bHasDiffusion?m_lNPos:1;
        std::vector<int> a(m_lNPos);
        std::iota (a.begin(), a.end(),0);
        std::for_each (__MODE__, std::begin(a), std::end(a), [&](int cpos){
            timekernel(pB1, pGr, pPos+3*cpos*step_diffusion, pB0+cpos*step_diffusion, td_gamma, pM0+3*cpos, e1, e2, pResult+3*cpos*m_lStepPos);
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
        std::complex<_T> *pB1,   // m_lNTime x 1 [Volt] 
        _T *pGr,                 // 3 x m_lNTime [Tesla/m] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
        _T td,                   // [second]
        _T *pB0,                 // m_lNPos x 1  [Tesla]
        _T *pPos,                // 3 x m_lNPos  [meter] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}                    
        _T T1, _T T2,            // [second]
        _T *pM0,                 // 3 x m_lNPos : column-maj
        int nPosition, 
        int nTime, 
        _T *pResult,            // 3 x (m_lNTime+1) x m_lNPos : column-major order {x1t0,y1t0,z1t0,...,x1tn,y1tn,z1tn,x2t0,y2t0,z2t0,...}, result equals m0 at t0
        bool hasDiffusion,
        bool saveAll            // return all time-points or only the final magnetization         
        )
{ 
    bloch bloch_obj(nPosition, nTime, hasDiffusion, saveAll);
    return bloch_obj.run(pB1, pGr, td, pB0, pPos, T1, T2, pM0, pResult);
}

} // extern

