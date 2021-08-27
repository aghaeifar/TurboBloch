/* Bloch simulation
 *
 * Author : Ali Aghaeifar <ali.aghaeifar.mri [at] gmail [dot] com>
 *
 * Inspired by Brian Hargreaves's mex bloch simulator
 *
 */

#include <stdio.h>
#include <math.h>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <thread>
#include <vector>
#include <windows.h>
#include <ppl.h>
//#include <mkl.h>
#include "bloch_sim.h"

#ifdef USE_GPU
//#include <cuda_runtime.h>
//#include <device_launch_parameters.h>
//#include <cuda.h>
//#include <cublas_v2.h>
#include "./gpu_matrix_mul/gpu_matrix_mul.h"
#endif


#define GAMMA_T 267522187.44
#define TWOPI	6.283185307179586

// Find the rotation matrix that rotates |n| radians about the vector given by nx,ny,nz
void calcrotmat(double nx, double ny, double nz, double rmat[9])
{
    double ar, ai, br, bi, hp, cp, sp;
    double arar, aiai, arai2, brbr, bibi, brbi2, arbi2, aibr2, arbr2, aibi2;
    double phi;

    phi = sqrt(nx*nx + ny*ny + nz*nz);

    if (phi == 0.0)
    {
        rmat[0] = 1; rmat[1] = 0; rmat[2] = 0;
        rmat[3] = 0; rmat[4] = 1; rmat[5] = 0;
        rmat[6] = 0; rmat[7] = 0; rmat[8] = 1;
    }
    else
    {
        // Cayley-Klein parameters
        hp = phi/2;
        cp = cos(hp);
        sp = sin(hp)/phi;	// /phi because [nx, ny, nz] is unit length in defs.
        ar = cp;
        ai = -nz*sp;
        br = ny*sp;
        bi = -nx*sp;

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
}

void timekernel(std::complex<double> *b1, double *gr,
                double *pr, double b0, double dt_gamma, double *m0,
                double e1, double e2, long nNTime, double *output)
{
    double rotx, roty, rotz, rotmat[9], m1[3];
    double e1_1 = e1 - 1;
    std::copy(m0, m0+3, output);
    for (int ct=0; ct<nNTime; ct++)
    {
        rotx = b1[ct].real() * dt_gamma * -1.0;
        roty = b1[ct].imag() * dt_gamma; // Hao Sun has changed this sign to '-', but I beleive the original '+' is correct.
        rotz = std::inner_product(gr, gr+3, pr, b0) * dt_gamma * -1.0; // -(gx*px + gy*py + gz*pz + b0) * dT * gamma
        gr += 3; // move to the next position

        calcrotmat(rotx, roty, rotz, rotmat);
        m1[0] = std::inner_product(rotmat+0, rotmat+3, output, 0.0) * e2;
        m1[1] = std::inner_product(rotmat+3, rotmat+6, output, 0.0) * e2;
        m1[2] = std::inner_product(rotmat+6, rotmat+9, output, 0.0) * e1 - e1_1;
        std::copy(m1, m1+3, output); // set magnetization for the next iteration
    }
}


bloch_sim::bloch_sim(long nPositions, long nTime, long nCoils)
{
    m_lNPos = nPositions;
    m_lNTime = nTime;
    m_lNCoils = nCoils;
    m_magnetization = new double[m_lNPos*3]();
    m_cb1 = new std::complex<double>[nPositions*nTime](); // () init to zero
}

bloch_sim::~bloch_sim()
{
    delete[] m_magnetization;
    delete[] m_cb1;
}


// ----------------------------------------------- //

bool bloch_sim::run(std::complex<double> *b1,   // m_lNTime x m_lNCoils : {t0c0, t1c0, t2c0,...,t0c1, t1c1, t2c1,...}
                    double *gr,                 // m_lNTime x 3 [Tesla/m] : {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
                    double tp,                  // m_lNTime x 1 [second]
                    double *b0,                 // m_lNPos x 1  [Radian]
                    double *pr,                 // m_lNPos x 3  [meter] : {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
                    double T1, double T2,       // [second]
                    std::complex<double> *sens, // m_lNCoils x m_lNPos : {c0p0, c1p0, c1p0,...,c0p1, c1p1, c1p1,...}
                    double *m0)                 // {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
{
    auto start = std::chrono::system_clock::now();
    if (b1==NULL || gr==NULL || b0==NULL || pr==NULL || sens==NULL || m0==NULL)
    {
        std::cout << "Program found at least one input is NULL"<<std::endl;
        return false;
    }

    // Calculate the E1 and E2 values at each time step.
    double e1 = exp(-tp/T1);
    double e2 = exp(-tp/T2);
    double tp_gamma = tp * GAMMA_T;

#ifndef USE_GPU
    Eigen::MatrixXcd e_b1comb = b1_uc * e_sens; // m_lNTime * m_lNPos
#endif

#ifdef USE_GPU
    // cuBLAS library uses column-major storage
    gpu_matrix_mul myGPU;
    myGPU.mul(b1, sens, m_cb1, m_lNTime, m_lNCoils, m_lNPos);
#endif

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
    std::cout<< "preparation " << elapsed.count() << " millisecond" << std::endl;
    start = std::chrono::system_clock::now();
    // =================== Do The Simulation! ===================
    // m_cb1 : {t0p0, t1p0, t2p0,... , t0p1, t1p1, t2p1, ...}
    concurrency::parallel_for (int(0), (int)m_lNPos, [&](int cpos){
 //   for (int cpos=0; cpos<(int)m_lNPos; cpos++)
        timekernel(m_cb1+cpos*m_lNTime, gr, pr+cpos*3, *(b0+cpos), tp_gamma, m0+cpos*3, e1, e2, m_lNTime, m_magnetization+cpos*3);
    });

    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
    std::cout<< "Simulation1 " << elapsed.count() << " millisecond" << std::endl;

/* // was very slow, don't know why
    std::thread *threadarr = new std::thread[m_lNPos];
    for (int cpos=0; cpos<m_lNPos; cpos++)
        threadarr[cpos] = std::thread(timekernel, m_cb1+cpos*m_lNTime, gr, pr+cpos*3, *(b0+cpos), tp_gamma, m0+cpos*3, e1, e2, m_lNTime, m_magnetization+cpos*3);
    for(int cpos=0; cpos<m_lNPos; cpos++)
        threadarr[cpos].join();
    delete[] threadarr;
*/
    return true;
}


// ----------------------------------------------- //
bool bloch_sim::getMagnetization(double result[])
{
    std::copy(m_magnetization, m_magnetization + m_lNPos*3, result);
    return true;
}

void bloch_sim::print(std::complex<double> *b1, double *gr, double tp, double *b0, double *pr, double T1, double T2, std::complex<double> *sens, double *m0)
{
    std::cout<<"\nb1"<<std::endl;
    for(int i=0; i<m_lNTime; i++)
    {
        for(int j=0; j<m_lNCoils; j++)
            std::cout<<(b1+i*m_lNCoils+j)->real() << "+i" << (b1+i*m_lNCoils+j)->imag() << "   ";
         std::cout<<std::endl;
    }

    std::cout<<"\ngr"<<std::endl;
    for(int i=0; i<m_lNTime; i++)
        std::cout<<*(gr+3*i)<< " "<<*(gr+3*i+1)<<" "<<*(gr+3*i+2)<<std::endl;

    std::cout<<"\nb0"<<std::endl;
    for(int i=0; i<m_lNPos; i++)
        std::cout<<*(gr+i)<<std::endl;

    std::cout<<"\npr"<<std::endl;
    for(int i=0; i<m_lNPos; i++)
        std::cout<<*(pr+3*i)<< " "<<*(pr+3*i+1)<<" "<<*(pr+3*i+2)<<std::endl;

    std::cout<<"\nsens"<<std::endl;
    for(int i=0; i<m_lNCoils; i++)
    {
        for(int j=0; j<m_lNPos; j++)
            std::cout<<(sens+i*m_lNPos+j)->real() << "+i" << (sens+i*m_lNPos+j)->imag() << "   ";
         std::cout<<std::endl;
    }

    std::cout<<"\nm0"<<std::endl;
    for(int i=0; i<m_lNPos; i++)
        std::cout<<*(m0+3*i)<< " "<<*(m0+3*i+1)<<" "<<*(m0+3*i+2)<<std::endl;
}





