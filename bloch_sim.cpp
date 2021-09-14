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
#include <mkl.h>
#include "tbb/parallel_for.h" // Intel Threading Building Blocks

#include "../eigen-3.4.0/Eigen/Dense"
#include "bloch_sim.h"

#define GAMMA_T 267522187.44

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

// q = {x, y, z, w}
void quat_rot(double q[4], double *m0, double *m1, bool noRF = false)
{
    double t[3];
    // see https://blog.molecular-matters.com/2013/05/24/a-faster-quaternion-vector-multiplication/
    if(noRF)
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

void apply_rot_quaternion(double nx, double ny, double nz, double *m0, double *m1)
{
    // creating quaternion rotation vector
    double q[4] = {0, 0, 0, 1}, phi;
    bool noRF = false;

    if(nx == 0 && ny == 0)
    {  // Only gradients and off-resonance, no RF
        noRF = true;
        phi  = abs(nz);
        q[2] = nz * sin(phi/2) / phi;
        q[3] = cos(phi/2);
    }
    else
    {
        phi = sqrt(nx*nx + ny*ny + nz*nz);
        double sp = sin(phi/2) / phi; // /phi because [nx, ny, nz] is unit length in defs.
        q[0] = nx * sp;
        q[1] = ny * sp;
        q[2] = nz * sp;
        q[3] = cos(phi/2);
    }
    quat_rot(q, m0, m1, noRF);
}

// ----------------------------------------------- //

void timekernel(std::complex<double> *b1, double *gr,
                double *pr, double b0, double dt_gamma, double *m0,
                double e1, double e2, long nNTime, double *output)
{
    double rotx, roty, rotz, m1[3], phi, sp;
    double e1_1 = e1 - 1;
    std::copy(m0, m0+3, output);

    if(e1 >= 0 && e2 >= 0) // including relaxations
    {
        for (int ct=0; ct<nNTime; ct++)
        {
            rotx = b1[ct].real() * dt_gamma * -1.0;
            roty = b1[ct].imag() * dt_gamma * -1.0;
            rotz = std::inner_product(gr, gr+3, pr, b0) * dt_gamma * -1.0; // -(gx*px + gy*py + gz*pz + b0) * dT * gamma
            gr += 3; // move to the next position

            //apply_rot_CayleyKlein(rotx, roty, rotz, output, m1);
            apply_rot_quaternion(-rotx, -roty, -rotz, output, m1); // quaternion needs additional sign reverse because looking down the axis of rotation, positive rotations appears clockwise

            m1[0] *= e2;
            m1[1] *= e2;
            m1[2]  = m1[2] * e1 - e1_1;

            std::copy(m1, m1+3, output); // set magnetization for the next iteration
        }
    }
    else // excluding relaxations, slightly faster
    {
        //cblas_zgemv
        Eigen::Quaterniond q0(1,0,0,0);
        Eigen::Quaterniond q1(1,0,0,0);

        for (int ct=0; ct<nNTime; ct++)
        {
            rotx = b1[ct].real() * dt_gamma * -1.0;
            roty = b1[ct].imag()  * dt_gamma * -1.0;
            rotz = std::inner_product(gr, gr+3, pr, b0) * dt_gamma * -1.0; // -(gx*px + gy*py + gz*pz + b0) * dT * gamma // ~40-50ms
            gr += 3; // move to the next position

            phi    = sqrt(rotx*rotx + roty*roty + rotz*rotz);
            sp     = sin(phi/2) / phi; // /phi because [nx, ny, nz] is unit length in defs.
            q1.x() = -rotx * sp;
            q1.y() = -roty * sp;
            q1.z() = -rotz * sp;
            q1.w() = cos(phi/2);
            q0     = q1 * q0; // ~20-30ms
        }
        double q[4] = {q0.x(), q0.y(), q0.z(), q0.w()};
        quat_rot(q, m0, output, false);
    }
}

// ----------------------------------------------- //

bloch_sim::bloch_sim(long nPositions, long nTime, long nCoils)
{
    m_lNPos   = nPositions;
    m_lNTime  = nTime;
    m_lNCoils = nCoils;
    m_dMagnetization = new double[m_lNPos*3];
    m_b1combined = NULL;
    if(m_lNCoils > 1)
        m_b1combined = new std::complex<double>[nPositions*nTime];
}

bloch_sim::~bloch_sim()
{
    delete[] m_dMagnetization;
    if(m_lNCoils > 1)
        delete[] m_b1combined;
}


// ----------------------------------------------- //

bool bloch_sim::run(std::complex<double> *b1,   // m_lNTime x m_lNCoils [Volt] : {t0c0, t1c0, t2c0,...,t0c1, t1c1, t2c1,...}
                    double *gr,                 // m_lNTime x 3 [Tesla/m] : {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
                    double td,                  // m_lNTime x 1 [second]
                    double *b0,                 // m_lNPos x 1  [Tesla]
                    double *pr,                 // m_lNPos x 3  [meter] : {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}                    
                    std::complex<double> *sens, // m_lNCoils x m_lNPos [Tesla/Volt]: {c0p0, c1p0, c1p0,...,c0p1, c1p1, c1p1,...}
                    double T1, double T2,       // [second]
                    double *m0)                 // {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
{
    if (b1==NULL || gr==NULL || b0==NULL || pr==NULL || m0==NULL)
    {
        std::cout << "At least one input is NULL. Terminate simulation..."<<std::endl;
        return false;
    }

    // Calculate the E1 and E2 values at each time step.
    double e1 = T1 < 0 ? -1.0 : exp(-td/T1);
    double e2 = T2 < 0 ? -1.0 : exp(-td/T2);
    double tp_gamma = td * GAMMA_T;

    if(sens != NULL && m_lNCoils>1)
    {
        MKL_Complex16 alpha, beta;
        alpha.real = 1.0; alpha.imag = 0.0;
        beta.real = 0.0; beta.imag = 0.0;
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m_lNTime, m_lNPos, m_lNCoils, &alpha, b1, m_lNTime, sens, m_lNCoils, &beta, m_b1combined, m_lNTime);
    }
    else
       m_b1combined = b1;

    // =================== Do The Simulation! ===================
    // b1combined : {t0p0, t1p0, t2p0,... , t0p1, t1p1, t2p1, ...}
    tbb::parallel_for(tbb::blocked_range<int>(0, m_lNPos), [&](tbb::blocked_range<int> r) {
        for (int cpos=r.begin(); cpos<r.end(); cpos++)
            timekernel(m_b1combined+cpos*m_lNTime, gr, pr+cpos*3, *(b0+cpos), tp_gamma, m0+cpos*3, e1, e2, m_lNTime, m_dMagnetization+cpos*3);
    });

    return true;
}


// ----------------------------------------------- //
bool bloch_sim::getMagnetization(double result[])
{
    std::copy(m_dMagnetization, m_dMagnetization + m_lNPos*3, result);
    return true;
}

void bloch_sim::print(std::complex<double> *b1, double *gr, double *b0, double *pr, std::complex<double> *sens, double *m0)
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
        std::cout<<*(b0+i)<<std::endl;

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





