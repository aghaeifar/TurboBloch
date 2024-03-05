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


#ifdef __SEQUENTIAL__
#define __MODE__  (std::execution::seq)
#else
#define __MODE__  (std::execution::par_unseq)
#endif


const long GAMMA = 267522187; 
const _T half = 0.5;

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
        s    = nz == 0? 0:sin(half * phi);
        q[0] = 0.;
        q[1] = 0.;
        q[2] = s * (nz>0 ? 1:-1); // equal to nz/phi == nz/abs(nz) == sign(nz)
    }
    else if(nz == 0)
    {  // Only RF, no gradients and off-resonance
        phi  = sqrt(nx*nx + ny*ny);
        s    = sin(half * phi);
        _T sp= s / phi; 
        q[0] = nx * sp;
        q[1] = ny * sp;
        q[2] = 0.; 
    }
    else
    {
        phi  = sqrt(nx*nx + ny*ny + nz*nz);
        s    = sin(half * phi);
        _T sp= s / phi; // /phi because [nx, ny, nz] is unit length in definition. This will be effective in the next lines, where nx, ny, nz * sp 
        q[0] = nx * sp;
        q[1] = ny * sp;
        q[2] = nz * sp;        
    }
    q[3] = cos(half * phi); // sqrt(1. - s*s) is faster
}

// ----------------------------------------------- //

void bloch::timekernel( const std::complex<_T> *b1, 
                        const _T *gr,
                        const _T *pr, 
                        const _T *b0, 
                        const _T *m0,
                        const _T T1, 
                        const _T T2, 
                        _T *output)
{
    _T m1[3], q[4];
    _T rotx, roty, rotz;
    _T e1, e2, dt;  
    std::copy(m0, m0+3, output);  // setting starting magnetization

    for (int ct = 0; ct < m_lNTime; ct++, gr += 3)
    {
        dt = m_dt[ct];
        // =================== p r e c e s s i o n ===================
        // rotations are right handed, thus all are negated.
        rotx = -b1[ct].real() * dt * GAMMA;
        roty = -b1[ct].imag() * dt * GAMMA;
        rotz = -std::inner_product(gr, gr + 3, pr, *b0) * dt * GAMMA; // -(gx*px + gy*py + gz*pz + b0) * dT * gamma

        create_quaternion(-rotx, -roty, -rotz, q); // quaternion needs additional sign reverse because looking down the axis of rotation, positive rotations appears clockwise
        apply_rot_quaternion(q, output, m1); // m1 = q*output 
        // =================== r e l a x a t i o n ===================
        e1 = exp(-dt/T1); e2 = exp(-dt/T2); 
        output += m_lStepTime;
        output[0] = m1[0] * e2;
        output[1] = m1[1] * e2;
        output[2] = m1[2] * e1 - e1 + 1;
    }
}

// ----------------------------------------------- //

bloch::bloch(size_t n_spatial_position,  // number of spatial positions or spins 
             size_t n_timepoints,        // number of time points 
             bool save_all_timepoints, // save final magnetization or whole evolution of spins
             bool is_relaxation_constant,
             bool is_dwelltime_constant
            )
{
    m_lNPos                 = n_spatial_position;
    m_lNTime                = n_timepoints;
    m_lStepPos              = save_all_timepoints ? (n_timepoints+1) : 1;
    m_lStepTime             = save_all_timepoints ? 3 : 0;
    m_bConstantT2T2         = is_relaxation_constant;
    m_bConstantDwellTime    = is_dwelltime_constant;
}

bloch::~bloch()
{

}


// ----------------------------------------------- //

bool bloch::run(const std::complex<_T> *pB1,// RF pulse [T]          ; n_timepoints x n_spatial_position
                const _T *pGr,              // gradients [T/m]       ; 3 x n_timepoints: column-major order {gx1,gy1,gz1,gx2,gy2,gz2,...,gxm,gym,gzm}
                const _T *dt,               // time step [Sec]       ; 1 x n_timepoints
                const _T *pB0,              // off-resonance [T]     ; 1 x n_spatial_position
                const _T *pPos,             // spatial positions [m] ; 3 x n_spatial_position: column-major order {x1,y1,z1,...,xm,ym,zm}
                const _T *T1,               // relaxations T1 [Sec]  ; 1 or 1 x n_spatial_position depends on "is_relaxation_constant"
                const _T *T2,               // relaxations T2 [Sec]  ; 1 or 1 x n_spatial_position depends on "is_relaxation_constant"                
                const _T *pM0,              // initial magnetization ; 3 x n_spatial_position : column-major order {x1,y1,z1,x2,y2,z2,...,xm,ym,zm}
                _T *pResult                 // output                ; 3 x n_spatial_position or 3 x (n_timepoints+1) x n_spatial_position, depends on "save_all_timepoints": column-major order {x1t1,y1t1,z1t1,...,x1tn,y1tn,z1tn,x2t1,y2t1,z2t1,...,x2tn,y2tn,z2tn,...,xmt1,ymt1,zmt1,...,xmtn,ymtn,zmtn}, result equals m0 at t0
                )                
{ 
    m_dt = const_cast<_T*>(dt);
    // =================== Do The Simulation! =================== 
    try
    {
        std::vector<int> a(m_lNPos);
        std::iota (a.begin(), a.end(),0);
        std::for_each (__MODE__, std::begin(a), std::end(a), [&](int cpos){
            timekernel( pB1 + cpos*m_lNTime, 
                        pGr, 
                        pPos + 3*cpos, 
                        pB0 + cpos, 
                        pM0 + 3*cpos, 
                        m_bConstantT2T2 ? *T1 : T1[cpos], 
                        m_bConstantT2T2 ? *T2 : T2[cpos], 
                        pResult + 3*cpos*m_lStepPos);
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
        std::complex<_T> *pB1,  // m_lNTime x n_spatial_position [Tesla] 
        _T *pGr,                // 3 x m_lNTime [Tesla/m] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
        _T *dt,                 // m_lNTime x 1 [second]
        _T *pB0,                // m_lNPos x 1  [Tesla]
        _T *pPos,               // 3 x m_lNPos  [meter] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}                    
        _T *T1, 
        _T *T2,                 // [second]
        _T *pM0,                // 3 x m_lNPos : column-maj
        int n_spatial_position, 
        int n_timepoints,         
        bool save_all_timepoints, // return all time-points or only the final magnetization  
        bool is_relaxation_constant,
        bool is_dwelltime_constant,
        _T *pResult             // 3 x (m_lNTime+1) x m_lNPos : column-major order {x1t0,y1t0,z1t0,...,x1tn,y1tn,z1tn,x2t0,y2t0,z2t0,...}, result equals m0 at t0       
        )
{ 
    bloch bloch_obj(n_spatial_position, n_timepoints, save_all_timepoints, is_relaxation_constant, is_dwelltime_constant);
    return bloch_obj.run(pB1, pGr, dt, pB0, pPos, T1, T2, pM0, pResult);
}

} // extern
