/* Bloch simulation
 *
 * Author : Ali Aghaeifar <ali.aghaeifar.mri [at] gmail [dot] com>
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
#include "bloch_sim.h"
#include "./gpu_matrix_mul/gpu_matrix_mul.h"


#define GAMMA_T 267522187.44
#define TWOPI	6.283185307179586

// Find the rotation matrix that rotates |n| radians about the vector given by nx,ny,nz

void apply_rot_CayleyKlein(double nx, double ny, double nz, double *m0, double *m1)
{
    double ar, ai, br, bi, hp, cp, sp;
    double arar, aiai, arai2, brbr, bibi, brbi2, arbi2, aibr2, arbr2, aibi2;
    double phi, rmat[9];

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

    m1[0] = std::inner_product(rmat+0, rmat+3, m0, 0.0);
    m1[1] = std::inner_product(rmat+3, rmat+6, m0, 0.0);
    m1[2] = std::inner_product(rmat+6, rmat+9, m0, 0.0);
}

void apply_rot_quaternion(double nx, double ny, double nz, double *m0, double *m1)
{
    // creating quaternion rotation vector
    double q[4], t[3];
    double phi = sqrt(nx*nx + ny*ny + nz*nz);
    double sp = sin(phi/2) / phi; // /phi because [nx, ny, nz] is unit length in defs.
    q[0] = nx * sp;
    q[1] = ny * sp;
    q[2] = nz * sp;
    q[3] = cos(phi/2);
    // see https://blog.molecular-matters.com/2013/05/24/a-faster-quaternion-vector-multiplication/
    t[0] = 2*(q[1]*m0[2] - q[2]*m0[1]);
    t[1] = 2*(q[2]*m0[0] - q[0]*m0[2]);
    t[2] = 2*(q[0]*m0[1] - q[1]*m0[0]);

    m1[0] = m0[0] + q[3]*t[0] + q[1]*t[2] - q[2]*t[1];
    m1[1] = m0[1] + q[3]*t[1] + q[2]*t[0] - q[0]*t[2];
    m1[2] = m0[2] + q[3]*t[2] + q[0]*t[1] - q[1]*t[0];
}

// Only gradients and off-resonance, no RF
void apply_rot_quaternion(double nz, double *m0, double *m1)
{
    // creating quaternion rotation vector
    double q[4], t[2];
    double phi = nz;
    double sp = sin(phi/2) / phi; // /phi because [nx, ny, nz] is unit length in defs.
    q[2] = nz * sp;
    q[3] = cos(phi/2);
    // see https://blog.molecular-matters.com/2013/05/24/a-faster-quaternion-vector-multiplication/
    t[0] =-2*q[2]*m0[1];
    t[1] = 2*q[2]*m0[0];

    m1[0] = m0[0] + q[3]*t[0] - q[2]*t[1];
    m1[1] = m0[1] + q[3]*t[1] + q[2]*t[0];
    m1[2] = m0[2];
}

void timekernel(std::complex<double> *b1, double *gr,
                double *pr, double b0, double dt_gamma, double *m0,
                double e1, double e2, long nNTime, double *output)
{
    double rotx, roty, rotz, m1[3];
    double e1_1 = e1 - 1;
    std::copy(m0, m0+3, output);
    for (int ct=0; ct<nNTime; ct++)
    {
        rotx = b1[ct].real() * dt_gamma * -1.0;
        roty = b1[ct].imag() * dt_gamma * -1.0;
        rotz = std::inner_product(gr, gr+3, pr, b0) * dt_gamma * -1.0; // -(gx*px + gy*py + gz*pz + b0) * dT * gamma
        gr += 3; // move to the next position

        //apply_rot_CayleyKlein(rotx, roty, rotz, output, m1);
        apply_rot_quaternion(-rotx, -roty, rotz, output, m1); // quaternion needs additional sign reverse because looking down the axis of rotation, positive rotations appears clockwise
        m1[0] *= e2;
        m1[1] *= e2;
        m1[2]  = m1[2] * e1 - e1_1;

        std::copy(m1, m1+3, output); // set magnetization for the next iteration
    }
}


bloch_sim::bloch_sim(long nPositions, long nTime, long nCoils)
{
    m_lNPos = nPositions;
    m_lNTime = nTime;
    m_lNCoils = nCoils;
    m_dMagnetization = new double[m_lNPos*3]();
    m_cdb1 = new std::complex<double>[nPositions*nTime](); // () init to zero
}

bloch_sim::~bloch_sim()
{
    delete[] m_dMagnetization;
    delete[] m_cdb1;
}


// ----------------------------------------------- //

bool bloch_sim::run(std::complex<double> *b1,   // m_lNTime x m_lNCoils [Volt] : {t0c0, t1c0, t2c0,...,t0c1, t1c1, t2c1,...}
                    double *gr,                 // m_lNTime x 3 [Tesla/m] : {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
                    double td,                  // m_lNTime x 1 [second]
                    double *b0,                 // m_lNPos x 1  [Radian]
                    double *pr,                 // m_lNPos x 3  [meter] : {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
                    double T1, double T2,       // [second]
                    std::complex<double> *sens, // m_lNCoils x m_lNPos [Tesla/Volt]: {c0p0, c1p0, c1p0,...,c0p1, c1p1, c1p1,...}
                    double *m0)                 // {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
{
    if (b1==NULL || sens==NULL)
    {
        std::cout << "Program found at least one input is NULL"<<std::endl;
        return false;
    }

    auto start = std::chrono::system_clock::now();

    // we can gain a lot in speed if we use float precision
    gpu_matrix_mul myGPU;
    myGPU.mul(b1, sens, m_cdb1, m_lNTime, m_lNCoils, m_lNPos);

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
    std::cout<< "preparation " << elapsed.count() << " millisecond" << std::endl;

    return run(m_cdb1, gr, td, b0, pr, T1, T2, m0);
}

bool bloch_sim::run(std::complex<double> *b1combined, double *gr, double td, double *b0, double *pr, double T1, double T2, double *m0)
{
    if (b1combined==NULL || gr==NULL || b0==NULL || pr==NULL || m0==NULL)
    {
        std::cout << "Program found at least one input is NULL"<<std::endl;
        return false;
    }
    auto start = std::chrono::system_clock::now();
    // Calculate the E1 and E2 values at each time step.
    double e1 = exp(-td/T1);
    double e2 = exp(-td/T2);
    double tp_gamma = td * GAMMA_T;
    // =================== Do The Simulation! ===================
    // b1combined : {t0p0, t1p0, t2p0,... , t0p1, t1p1, t2p1, ...}
    concurrency::parallel_for (int(0), (int)m_lNPos, [&](int cpos){
//    for (int cpos=0; cpos<(int)m_lNPos; cpos++){
        timekernel(b1combined+cpos*m_lNTime, gr, pr+cpos*3, *(b0+cpos), tp_gamma, m0+cpos*3, e1, e2, m_lNTime, m_dMagnetization+cpos*3);
    }
    );

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
    std::cout<< "Simulation " << elapsed.count() << " millisecond" << std::endl;

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
    std::copy(m_dMagnetization, m_dMagnetization + m_lNPos*3, result);
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





