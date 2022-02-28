/* Bloch simulation
 *
 * Author : Ali Aghaeifar <ali.aghaeifar.mri [at] gmail [dot] com>
 *
 */
#include <iostream>
#include <numeric>
#include <chrono>
#include <thread>
#include <mkl.h>

#ifdef _TBB
#include "tbb/parallel_for.h" // Intel Threading Building Blocks
#elif _WIN32
#include <ppl.h>
#include <windows.h> // without this ppl will be single thread
#endif

#include "bloch.h"

#define GAMMA_T 267522187.44

using namespace std;


// Find the rotation matrix that rotates |n| radians about the vector given by nx,ny,nz
void apply_rot_CayleyKlein(float ar, float ai, float br, float bi, float *m0, float *m1)
{
    float arar, aiai, arai2, brbr, bibi, brbi2, arbi2, aibr2, arbr2, aibi2;
    float  rmat[9];

    if (ar == 1 && ai == 0 && br == 0 &&  bi == 0 )
    {
        rmat[0] = 1; rmat[1] = 0; rmat[2] = 0;
        rmat[3] = 0; rmat[4] = 1; rmat[5] = 0;
        rmat[6] = 0; rmat[7] = 0; rmat[8] = 1;
    }
    else
    {
        // Auxiliary variables to speed this up
        arar  = ar*ar;
        aiai  = ai*ai;
        arai2 = 2*ar*ai;
        brbr  = br*br;
        bibi  = bi*bi;
        brbi2 = 2*br*bi;
        arbi2 = 2*ar*bi;
        aibr2 = 2*ai*br;
        arbr2 = 2*ar*br;
        aibi2 = 2*ai*bi;

        // Make rotation matrix.
        rmat[0] =  arar  -aiai -brbr +bibi;
        rmat[1] = -arai2 -brbi2; // arai2 + brbi2
        rmat[2] = -arbr2 +aibi2; // arbr2 +aibi2
        rmat[3] =  arai2 -brbi2;
        rmat[4] =  arar  -aiai +brbr -bibi;
        rmat[5] = -aibr2 -arbi2; // -aibr2 +arbi2
        rmat[6] =  arbr2 +aibi2;
        rmat[7] =  arbi2 -aibr2;
        rmat[8] =  arar  +aiai -brbr -bibi;
    }

    m1[0] = std::inner_product(rmat+0, rmat+3, m0, 0.0);
    m1[1] = std::inner_product(rmat+3, rmat+6, m0, 0.0);
    m1[2] = std::inner_product(rmat+6, rmat+9, m0, 0.0);
}

void create_CayleyKlein(float nx, float ny, float nz, complex<float> &alpha, complex<float> &beta)
{
    // Cayley-Klein parameters - See Paulies paper "Parameter Relations for the Shinnar-Le Roux Selective Excitation Pulse Design Algorithm"
    float phi = sqrt(nx*nx + ny*ny + nz*nz);
    if(phi == 0) // this if is just for speedup 
    {
        alpha.real(1);
        alpha.imag(0);
        beta.real(0);
        beta.imag(0);
    }
    else
    {
        float hp = phi / 2;
        float sp = sin(hp) / phi; // /phi because [nx, ny, nz] is unit length in defs.

        alpha.real(cos(hp));
        alpha.imag(-nz * sp);
        beta.real(+ny * sp); // this is +
        beta.imag(-nx * sp);
    }
}

// q = {x, y, z, w}
void apply_rot_quaternion(float q[4], float *m0, float *m1)
{
    float t[3];
    // see https://blog.molecular-matters.com/2013/05/24/a-faster-quaternion-vector-multiplication/
    if(q[0] == 0 && q[1] == 0) // no RF pulse case
    {
        t[0] =-2*q[2]*m0[1];
        t[1] = 2*q[2]*m0[0];

        m1[0] = m0[0] + q[3]*t[0] - q[2]*t[1];
        m1[1] = m0[1] + q[3]*t[1] + q[2]*t[0];
        m1[2] = m0[2];
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

void create_quaternion(float nx, float ny, float nz, float q[4])
{
    float phi = 0.;
    // q = sin(θ/2)(xi+yj+zk) + cos(θ/2)
    if(nx == 0 && ny == 0)
    {  // Only gradients and off-resonance, no RF
        phi  = abs(nz);
        q[0] = 0;
        q[1] = 0;
        q[2] = sin(phi/2) * (nz>0 ? 1:-1); // equal to nz/phi == nz/abs(nz) == sign(nz)
    }
    else
    {
        phi = sqrt(nx*nx + ny*ny + nz*nz);
        float sp = sin(phi/2) / phi; // /phi because [nx, ny, nz] is unit length in definition. This will be effective in the next lines, where nx, ny, nz * sp 
        q[0] = nx * sp;
        q[1] = ny * sp;
        q[2] = nz * sp;        
    }
    q[3] = cos(phi/2);
}

// ----------------------------------------------- //

void timekernel(std::complex<float> *b1xy, 
                float *gr,
                float *pr, 
                float b0, 
                float td_gamma, 
                float *m0,
                float e1, float e2, 
                long nNTime, 
                float *output)
{
    float m1[3]={0.,0.,1.}, q[4];
    float e1_1 = e1 - 1;    

    if(e1 >= 0 && e2 >= 0) // including relaxations
    {
        float rotx, roty, rotz;
        std::complex<float> alpha0, beta0;
        for (int ct=0; ct<nNTime; ct++)
        {            
            // rotations are right handed, thus all are negated.
            rotx = -b1xy[ct].real() * td_gamma;
            roty = -b1xy[ct].imag() * td_gamma;
            
            rotz = -std::inner_product(gr, gr+3, pr, b0) * td_gamma; // -(gx*px + gy*py + gz*pz + b0) * dT * gamma
            gr += 3; // move to the next time point
            
            //quaternion method is faster
            //create_CayleyKlein(rotx, roty, rotz, alpha0, beta0);
            //apply_rot_CayleyKlein(alpha0.real(), alpha0.imag(), beta0.real(), beta0.imag(), m0, m1);
            create_quaternion(-rotx, -roty, -rotz, q); // quaternion needs additional sign reverse because looking down the axis of rotation, positive rotations appears clockwise
            apply_rot_quaternion(q, m0, m1);
            m1[0] *= e2;
            m1[1] *= e2;
            m1[2]  = m1[2] * e1 - e1_1;

            std::copy(m1, m1+3, m0); // set magnetization for the next iteration
        }
    }
    else // excluding relaxations, faster calculation
    {
        complex<float> alpha0(1, 0), beta0(0, 0), alpha1, beta1, alpha2, beta2;
        // consider https://www.intel.com/content/www/us/en/developer/articles/technical/onemkl-improved-small-matrix-performance-using-just-in-time-jit-code.html
        //apply_rot_CayleyKlein(ar2, ai2, br2, bi2, m0, m1);  
    }

    std::copy(m1, m1+3, output);    
}

// ----------------------------------------------- //

bloch::bloch(long nPosition, long nTime, long nCoil)
{
    m_lNPos   = nPosition;
    m_lNTime  = nTime;
    m_lNCoil = nCoil;
    m_dMagnetization = new float[3*m_lNPos];    
    m_b1combined = new std::complex<float>[nTime*nPosition];
}

bloch::~bloch()
{
    delete[] m_dMagnetization;
    delete[] m_b1combined;
}


// ----------------------------------------------- //

bool bloch::run(std::complex<float> *b1,       // m_lNTime x m_lNCoil [Volt] : column-major order {t0c0, t1c0, t2c0,...,t0c1, t1c1, t2c1,...}
                    float *gr,                 // 3 x m_lNTime [Tesla/m] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
                    float td,                  // m_lNTime x 1 [second]
                    float *b0,                 // m_lNPos x 1  [Tesla]
                    float *pr,                 // 3 x m_lNPos  [meter] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}                    
                    std::complex<float> *sens, // m_lNCoil x m_lNPos [Tesla/Volt]: column-major order {c0p0, c1p0, c1p0,...,c0p1, c1p1, c1p1,...}
                    float T1, float T2,        // [second]
                    float *m0)                 // 3 x m_lNPos : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
{
    if (b1==NULL || gr==NULL || b0==NULL || pr==NULL || m0==NULL)
    {
        std::cout << "At least one input is NULL. Terminate simulation..."<<std::endl;
        return false;
    }
    
    // Calculate the E1 and E2 values at each time step.
    float e1 = T1 < 0. ? -1.0 : exp(-td/T1);
    float e2 = T2 < 0. ? -1.0 : exp(-td/T2);
    float td_gamma = td * GAMMA_T;
    
    if(sens != NULL)
    {
        MKL_Complex8 alpha, beta; // float precision complex values. MKL_Complex16 for double precision.
        alpha.real = 1.0; alpha.imag = 0.0;
        beta.real = 0.0; beta.imag = 0.0;
        // consider gemm3m for faster calculation but higher numerical rounding errors
        cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m_lNTime, m_lNPos, m_lNCoil, &alpha, b1, m_lNTime, sens, m_lNCoil, &beta, m_b1combined, m_lNTime);
    }
    else
    {
        for(int cpos=0; cpos<m_lNPos; cpos++)
            std::copy(b1, b1+m_lNTime, m_b1combined+cpos*m_lNTime);
    }

    // int nTimes_nPos = m_lNPos * m_lNTime;
    // float *b1xy = (float *) m_b1combined;
    // cblas_scopy(nTimes_nPos, m_b1x, 1, b1temp, 2);
    // =================== Do The Simulation! ===================
    auto start = std::chrono::system_clock::now();  
    try
    {
        // b1combined : {t0p0, t1p0, t2p0,... , t0p1, t1p1, t2p1, ...} 
#ifdef _TBB
    tbb::parallel_for(tbb::blocked_range<int>(0, m_lNPos), [&](tbb::blocked_range<int> r) {
        for (int cpos=r.begin(); cpos<r.end(); cpos++)
            timekernel(m_b1combined+cpos*m_lNTime, gr, pr+cpos*3, *(b0+cpos), td_gamma, m0+cpos*3, e1, e2, m_lNTime, m_dMagnetization+cpos*3);
    });
#elif _OPENMP     
    #pragma omp parallel for
    for(int cpos=0; cpos<m_lNPos; cpos++)        
        timekernel(m_b1combined+cpos*m_lNTime, gr, pr+cpos*3, *(b0+cpos), td_gamma, m0+cpos*3, e1, e2, m_lNTime, m_dMagnetization+cpos*3);
#elif _WIN32 
    concurrency::parallel_for(int(0), (int)m_lNPos, [&](int cpos){
        timekernel(m_b1combined +cpos*m_lNTime, gr, pr+cpos*3, *(b0+cpos), td_gamma, m0+cpos*3, e1, e2, m_lNTime, m_dMagnetization+cpos*3);
    });
#else
    for (int cpos=0; cpos<m_lNPos; cpos++)  // sequential   
            timekernel(m_b1combined+cpos*m_lNTime, gr, pr+cpos*3, *(b0+cpos), td_gamma, m0+cpos*3, e1, e2, m_lNTime, m_dMagnetization+cpos*3);
#endif

    }
    catch( std::exception &ex )
    {
        std::cout << "Simulation failed." << std::endl;
        std::cout << ex.what() << std::endl;
        return false;
    }

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
    std::cout<< "Simulation -> " << elapsed.count() << " millisecond" << std::endl;        
    return true;
}


// ----------------------------------------------- //
// 3 x m_lNPos : {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
bool bloch::getMagnetization(float result[])
{
    std::copy(m_dMagnetization, m_dMagnetization + m_lNPos*3, result);
    return true;
}



extern "C" {
        bool bloch_sim(
        std::complex<float> *pb1,   // m_lNTime x m_lNCoils [Volt] : column-major order {t0c0, t1c0, t2c0,...,t0c1, t1c1, t2c1,...}
        float *pgr,                 // 3 x m_lNTime [Tesla/m] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
        float td,                   // m_lNTime x 1 [second]
        float *pb0,                 // m_lNPos x 1  [Tesla]
        float *ppr,                 // 3 x m_lNPos  [meter] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}                    
        std::complex<float> *psens, // m_lNCoils x m_lNPos [Tesla/Volt]: column-major order {c0p0, c1p0, c1p0,...,c0p1, c1p1, c1p1,...}
        float T1, float T2,         // [second]
        float *pm0,                 // 3 x m_lNPos : column-maj
        long nPosition, 
        long nTime, 
        long nCoil,
        float *presult              // 3 x m_lNPos  [meter] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}  
        )            
{
    bloch bloch_obj(nPosition, nTime, nCoil);
    if(bloch_obj.run(pb1, pgr, td, pb0, ppr, psens, T1, T2, pm0) == false)
        return false;

    bloch_obj.getMagnetization(presult);
    return true;
}
} // extern