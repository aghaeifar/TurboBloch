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
void apply_rot_CayleyKlein(double ar, double ai, double br, double bi, double *m0, double *m1)
{
    double arar, aiai, arai2, brbr, bibi, brbi2, arbi2, aibr2, arbr2, aibi2;
    double  rmat[9];

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

void create_CayleyKlein(double nx, double ny, double nz, complex<double> &alpha, complex<double> &beta)
{
    // Cayley-Klein parameters - See Paulies paper "Parameter Relations for the Shinnar-Le Roux Selective Excitation Pulse Design Algorithm"
    double phi = sqrt(nx*nx + ny*ny + nz*nz);
    if(phi == 0) // this if is just for speedup 
    {
        alpha.real(1);
        alpha.imag(0);
        beta.real(0);
        beta.imag(0);
    }
    else
    {
        double hp = phi / 2;
        double sp = sin(hp) / phi; // /phi because [nx, ny, nz] is unit length in defs.

        alpha.real(cos(hp));
        alpha.imag(-nz * sp);
        beta.real(+ny * sp); // this is +
        beta.imag(-nx * sp);
    }
}

// q = {x, y, z, w}
void apply_rot_quaternion(double q[4], double *m0, double *m1)
{
    double t[3];
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

void create_quaternion(double nx, double ny, double nz, double q[4])
{
    double phi = 0.;
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
        double sp = sin(phi/2) / phi; // /phi because [nx, ny, nz] is unit length in definition. This will be effective in the next lines, where nx, ny, nz * sp 
        q[0] = nx * sp;
        q[1] = ny * sp;
        q[2] = nz * sp;        
    }
    q[3] = cos(phi/2);
}

// ----------------------------------------------- //

void bloch::timekernel( std::complex<double> *b1xy, 
                        double *gr,
                        double *pr, 
                        double b0, 
                        double td_gamma, 
                        double *m0,
                        double e1, double e2, 
                        double *output)
{
    double m1[3], q[4];
    double e1_1 = e1 - 1;  
    std::copy(m0, m0+3, output);  // setting starting magnetization
    
    if(e1 >= 0 && e2 >= 0) // including relaxations
    {
        double rotx, roty, rotz;
        std::complex<double> alpha0, beta0;
        for (int ct=0; ct<m_lNTime; ct++)
        {            
            // rotations are right handed, thus all are negated.
            rotx = -b1xy[ct].real() * td_gamma;
            roty = -b1xy[ct].imag() * td_gamma;
            
            rotz = -std::inner_product(gr, gr+3, pr, b0) * td_gamma; // -(gx*px + gy*py + gz*pz + b0) * dT * gamma
            gr += 3; // move to the next time point
            
            //quaternion method is faster
            // create_CayleyKlein(rotx, roty, rotz, alpha0, beta0);
            // apply_rot_CayleyKlein(alpha0.real(), alpha0.imag(), beta0.real(), beta0.imag(), output, m1);
            create_quaternion(-rotx, -roty, -rotz, q); // quaternion needs additional sign reverse because looking down the axis of rotation, positive rotations appears clockwise
            apply_rot_quaternion(q, output, m1);

            output += m_lStepTime;
            output[0] = m1[0] * e2;
            output[1] = m1[1] * e2;
            output[2] = m1[2] * e1 - e1_1;            
        }
    }
    else // excluding relaxations, faster calculation
    {
        complex<double> alpha0(1, 0), beta0(0, 0), alpha1, beta1, alpha2, beta2;
        // consider https://www.intel.com/content/www/us/en/developer/articles/technical/onemkl-improved-small-matrix-performance-using-just-in-time-jit-code.html
        //apply_rot_CayleyKlein(ar2, ai2, br2, bi2, m0, m1);  
    } 
}

// ----------------------------------------------- //

bloch::bloch(long nPosition, long nTime, long nCoil, bool saveAll)
{
    m_lNPos   = nPosition;
    m_lNTime  = nTime;
    m_lNCoil  = nCoil;
    m_lStepPos  = saveAll ? 3 * (nTime+1) : 3;
    m_lStepTime = saveAll ? 3 : 0;
    m_dB1combined = new std::complex<double>[nTime*nPosition];
}

bloch::~bloch()
{
    delete[] m_dB1combined;
}


// ----------------------------------------------- //

bool bloch::run(std::complex<double> *pB1,   // m_lNTime x m_lNCoil [Volt]: column-major order {t0c0, t1c0, t2c0,...,t0c1, t1c1, t2c1,...}
                double *pGr,                 // 3 x m_lNTime [Tesla/m] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
                double td,                   // [second]
                double *pB0,                 // m_lNPos x 1  [Tesla]
                double *pPos,                // 3 x m_lNPos  [meter] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
                std::complex<double> *pSens, // m_lNCoil x m_lNPos [Tesla/Volt]: column-major order {c0p0, c1p0, c2p0,...,c0p1, c1p1, c2p1,...}
                double T1, double T2,         // [second]
                double *pM0,                 // 3 x m_lNPos : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
                double *pResult)             // 3 x (m_lNTime+1) x m_lNPos : column-major order {x1t0,y1t0,z1t0,...,x1tn,y1tn,z1tn,x2t0,y2t0,z2t0,...}, result equals m0 at t0
{
    // Calculate the E1 and E2 values.
    double e1 = T1 < 0. ? -1.0 : exp(-td/T1);
    double e2 = T2 < 0. ? -1.0 : exp(-td/T2);
    double td_gamma = td * GAMMA_T;
    
    if(pSens != NULL)
    {
        MKL_Complex16 alpha, beta; // MKL_Complex8 for double precision. MKL_Complex16 for double precision.
        alpha.real = 1.0; alpha.imag = 0.0;
        beta.real = 0.0; beta.imag = 0.0;
        // consider gemm3m for faster calculation but higher numerical rounding errors
        // cblas_cgemm for complex float
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m_lNTime, m_lNPos, m_lNCoil, &alpha, pB1, m_lNTime, pSens, m_lNCoil, &beta, m_dB1combined, m_lNTime);
    }
    else
    {
        for(int cpos=0; cpos<m_lNPos; cpos++)
            std::copy(pB1, pB1+m_lNTime, m_dB1combined + cpos*m_lNTime);
    }
    
    // =================== Do The Simulation! ===================
    // auto start = std::chrono::system_clock::now();  
    try
    {
        // b1combined : {t0p0, t1p0, t2p0,... , t0p1, t1p1, t2p1, ...} 
#ifdef _TBB
    tbb::parallel_for(tbb::blocked_range<int>(0, m_lNPos), [&](tbb::blocked_range<int> r) {
        for (int cpos=r.begin(); cpos<r.end(); cpos++)
            timekernel(m_dB1combined+cpos*m_lNTime, pGr, pPos+cpos*3, *(pB0+cpos), td_gamma, pM0+cpos*3, e1, e2, pResult+cpos*m_lStepPos);
    });
#elif _OPENMP     
    #pragma omp parallel for
    for(int cpos=0; cpos<m_lNPos; cpos++)        
        timekernel(m_dB1combined+cpos*m_lNTime, pGr, pPos+cpos*3, *(pB0+cpos), td_gamma, pM0+cpos*3, e1, e2, pResult+cpos*m_lStepPos);
#elif _WIN32 
    concurrency::parallel_for(int(0), (int)m_lNPos, [&](int cpos){
        timekernel(m_dB1combined + cpos*m_lNTime, pGr, pPos + cpos*3, *(pB0+cpos), td_gamma, pM0 + cpos*3, e1, e2, pResult + cpos*m_lStepPos);
    });
#else
    for (int cpos=0; cpos<m_lNPos; cpos++)  // sequential   
            timekernel(m_dB1combined+cpos*m_lNTime, pGr, pPos+cpos*3, *(pB0+cpos), td_gamma, pM0+cpos*3, e1, e2, pResult+cpos*m_lStepPos);
#endif

    }
    catch( std::exception &ex )
    {
        std::cout << "Simulation failed." << std::endl;
        std::cout << ex.what() << std::endl;
        return false;
    }

    // auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
    // std::cout<< "Simulation -> " << elapsed.count() << " millisecond" << std::endl;        
    return true;
}


extern "C" {
        bool bloch_sim(
        std::complex<double> *pB1,   // m_lNTime x m_lNCoils [Volt] : column-major order {t0c0, t1c0, t2c0,...,t0c1, t1c1, t2c1,...}
        double *pGr,                 // 3 x m_lNTime [Tesla/m] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
        double td,                   // [second]
        double *pB0,                 // m_lNPos x 1  [Tesla]
        double *pPos,                // 3 x m_lNPos  [meter] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}                    
        std::complex<double> *pSens, // m_lNCoils x m_lNPos [Tesla/Volt]: column-major order {c0p0, c1p0, c1p0,...,c0p1, c1p1, c1p1,...}
        double T1, double T2,         // [second]
        double *pM0,                 // 3 x m_lNPos : column-maj
        long nPosition, 
        long nTime, 
        long nCoil,
        double *pResult,             // 3 x (m_lNTime+1) x m_lNPos : column-major order {x1t0,y1t0,z1t0,...,x1tn,y1tn,z1tn,x2t0,y2t0,z2t0,...}, result equals m0 at t0
        bool saveAll)               // return all time-points or only the final magnetization         
{
    bloch bloch_obj(nPosition, nTime, nCoil, saveAll);
    return bloch_obj.run(pB1, pGr, td, pB0, pPos, pSens, T1, T2, pM0, pResult);
}
} // extern