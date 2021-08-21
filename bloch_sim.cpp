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
#include "bloch_sim.h"

#if defined(_MSC_VER) && defined(USE_PPL)
#include <windows.h>
#include <ppl.h>
#endif

#define GAMMA_T 267522187.44
#define TWOPI	6.283185307179586

// Find the rotation matrix that rotates |n| radians about the vector given by nx,ny,nz
void calcrotmat(double nx, double ny, double nz, Eigen::Matrix3d &rmat)
{
    double ar, ai, br, bi, hp, cp, sp;
    double arar, aiai, arai2, brbr, bibi, brbi2, arbi2, aibr2, arbr2, aibi2;
    double phi;

    phi = sqrt(nx*nx + ny*ny + nz*nz);

    if (phi == 0.0)
    {
        rmat.setZero();
        rmat.diagonal() << 1.0, 1.0, 1.0;
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
        rmat(0,0) =  arar  -aiai -brbr +bibi;
        rmat(0,1) = -arai2 -brbi2;
        rmat(0,2) = -arbr2 +aibi2;
        rmat(1,0) =  arai2 -brbi2;
        rmat(1,1) =  arar  -aiai +brbr -bibi;
        rmat(1,2) = -aibr2 -arbi2;
        rmat(2,0) =  arbr2 +aibi2;
        rmat(2,1) =  arbi2 -aibr2;
        rmat(2,2) =  arar  +aiai -brbr -bibi;
    }
}

void calcrotmat(double nx, double ny, double nz, double rmat[][3])
{
    double ar, ai, br, bi, hp, cp, sp;
    double arar, aiai, arai2, brbr, bibi, brbi2, arbi2, aibr2, arbr2, aibi2;
    double phi;

    phi = sqrt(nx*nx + ny*ny + nz*nz);

    if (phi == 0.0)
    {
        rmat[0][0] = 1; rmat[0][1] = 0; rmat[0][2] = 0;
        rmat[1][0] = 0; rmat[1][1] = 1; rmat[1][2] = 0;
        rmat[2][0] = 0; rmat[2][1] = 0; rmat[2][2] = 1;
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
        rmat[0][0] =  arar  -aiai -brbr +bibi;
        rmat[0][1] = -arai2 -brbi2;
        rmat[0][2] = -arbr2 +aibi2;
        rmat[1][0] =  arai2 -brbi2;
        rmat[1][1] =  arar  -aiai +brbr -bibi;
        rmat[1][2] = -aibr2 -arbi2;
        rmat[2][0] =  arbr2 +aibi2;
        rmat[2][1] =  arbi2 -aibr2;
        rmat[2][2] =  arar  +aiai -brbr -bibi;
    }
}


bloch_sim::bloch_sim()
{
    m_lNTime = m_lNPos = m_lNCoils = 1;
    m_result.resize(m_lNPos, vector<vector<double>>(m_lNTime, {0.0, 0.0, 0.0}));
}

bool bloch_sim::runkernel(MatrixXcd &e_b1, MatrixXd &e_gr, VectorXd &e_tp, VectorXd &e_b0, MatrixXcd &e_sens, MatrixXd &e_pr, MatrixXd &e_m0, double T1, double T2)
{
    Eigen::setNbThreads(6);
    // Resize m_result
    m_result_x.resize(m_lNPos, m_lNTime);
    m_result_y.resize(m_lNPos, m_lNTime);
    m_result_z.resize(m_lNPos, m_lNTime);

    // Calculate the E1 and E2 values at each time step.
    Eigen::VectorXd e_e1 = (e_tp/T1).array().exp();
    Eigen::VectorXd e_e2 = (e_tp/T2).array().exp();

    auto start = std::chrono::system_clock::now();
    // Combined B1    
    Eigen::MatrixXcd e_b1comb = e_b1 * e_sens.transpose(); // m_lNTime * m_lNPos
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
    std::cout<< "Preparation time1 " << elapsed.count() << " millisecond" << std::endl;

    e_b0 = e_b0 * TWOPI; // from Hz to rad;

    start = std::chrono::system_clock::now();
    // -(gx*px + gy*py + gz*pz + b0) * dT * gamma
    Eigen::MatrixXd rotz = ((e_gr * e_pr.transpose()).rowwise() + e_b0.transpose()).array().colwise() * e_tp.array() * -1.0 * GAMMA_T; // m_lNTime * m_lNPos
    Eigen::MatrixXd rotx = e_b1comb.real().array().colwise() * e_tp.array() * -1.0 * GAMMA_T;
    Eigen::MatrixXd roty = e_b1comb.imag().array().colwise() * e_tp.array() * GAMMA_T; // Hao Sun has changed this sign to '-', but I beleive the original '+' is correct.

    MatrixXd e_m0t = e_m0.transpose(); // m_lNPos x 3 --> 3 x m_lNPos
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
    std::cout<< "Preparation time2 " << elapsed.count() << " millisecond" << std::endl;
    start = std::chrono::system_clock::now();
#if defined(_MSC_VER) && defined(USE_PPL)
    std::cout << "Using PPL"<<std::endl;
    concurrency::parallel_for (int(0), (int)m_lNPos, [&](int cpos){
#else
#ifdef _OPENMP
    std::cout << "Using OpenMP"<<std::endl;
    #pragma omp parallel
    #pragma omp for
#endif
    for (int cpos=0; cpos<(int)m_lNPos; cpos++){
#endif
//        Eigen::Matrix3d rotmat;
//        Eigen::VectorXd m1;
        double rotmat[3][3], m1[3];
        for (int ct=0; ct<(int)m_lNTime; ct++)
        {
            // I don't use Eigen here because OpenMP and MKL don't work well toghether
            calcrotmat(rotx(ct,cpos), roty(ct,cpos), rotz(ct,cpos), rotmat);
            // m0: Current magnetization before rotation.
            //m1 = rotmat * e_m0t.col(cpos);

            m1[0] = rotmat[0][0]*e_m0t.col(cpos)(0) + rotmat[0][1]*e_m0t.col(cpos)(1) + rotmat[0][2]*e_m0t.col(cpos)(2);
            m1[1] = rotmat[1][0]*e_m0t.col(cpos)(0) + rotmat[1][1]*e_m0t.col(cpos)(1) + rotmat[1][2]*e_m0t.col(cpos)(2);
            m1[2] = rotmat[2][0]*e_m0t.col(cpos)(0) + rotmat[2][1]*e_m0t.col(cpos)(1) + rotmat[2][2]*e_m0t.col(cpos)(2);
            // Decay
            m_result_x(cpos, ct) = e_m0t(0,cpos) = m1[0] * e_e2(ct);
            m_result_y(cpos, ct) = e_m0t(1,cpos) = m1[1] * e_e2(ct);
            m_result_z(cpos, ct) = e_m0t(2,cpos) = m1[2] * e_e1(ct) + 1.0 - e_e1(ct);
        }
    }
#if defined(_MSC_VER) && defined(USE_PPL)
    );
#endif
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
    std::cout<< "Simulation time " << elapsed.count() << " millisecond" << std::endl;
    return true;
}


bool bloch_sim::prep(vector<vector<complex<double>>> &b1,
                     vector<vector<double>> &gr,
                     vector<double> &tp,
                     vector<double> &b0,
                     vector<vector<double>> &pr,
                     vector<vector<complex<double>>> &sens,
                     vector<vector<double>> &m0)
{
    m_lNTime  = b1.size();   // m_lNTime x m_lNCoils
    m_lNPos   = b0.size();   // m_lNPos x 1
    m_lNCoils = sens.empty()?1:sens[0].size(); // m_lNPos x m_lNCoils

    // =============== Initialize Unset Inputs =================
    if(m0.empty())
        m0.resize(m_lNPos, {0.0, 0.0, 1.0}); // Set magnetization to Equilibrium

    if(sens.empty())
    {   // create only one coil
        complex<double> cTemp(1.0, 0.0);
        sens.resize(m_lNPos, vector<complex<double>>(1, cTemp));
    }

    // m_b1: m_lNTime x 1 --> m_lNTime x m_lNCoils
    if (b1[0].size() == 1 && m_lNCoils>1)
        for(size_t i=0; i<m_lNTime; i++)
            b1[i].resize(m_lNCoils, b1[i][0]);

    // ====================== Input Size =========================
    //cout<<"Time points = "<<m_lNTime<<endl;
    if(m_lNTime != gr.size())
    {
        cout<<"Mismatch in number of time points"<<endl;
        return false;
    }

    //cout<<"Number of positions = "<<m_lNPos<<endl;
    if(m_lNPos != pr.size() ||
            m_lNPos != sens.size() ||
            m_lNPos != m0.size() )
    {
        cout<<"Mismatch in number of positions"<<endl;
        return false;
    }

    //cout<<"Number of coils = "<<m_lNCoils<<endl;
    if(m_lNCoils != b1[0].size())
    {
        cout<<"Mismatch in number of coils"<<endl;
        return false;
    }

    if(gr[0].size() != 3 ||  m0[0].size() != 3 || pr[0].size() != 3 )
    {
        cout<<"Unacceptable input size."<<endl;
        return false;
    }

    // ====================== Time points ======================

    // Single value given --> this is the interval length for all.
    if (tp.size() == 1)
        tp.resize(m_lNTime, tp[0]);
    // Monotonically INCREASING list of end times given --> calculate intervals
    else if (tp.size() == m_lNTime)
    {
        vector<double> tp_temp = tp;
        std::adjacent_difference(tp_temp.begin(), tp_temp.end(), tp.begin());
        for (size_t i=0; i<tp.size(); i++)
            if(tp[i] < 0)
            {
                cout<<"Time points are not monotonically increasing!"<<endl;
                return false;
            }
    }
    else
    {
        cout<<"Time-point length differs from B1 length."<<endl;
        return false;
    }

    return true;
}


bool bloch_sim::run( vector<vector<complex<double>>> b1,// m_lNTime x m_lNCoils
                     vector<vector<double>> gr,         // m_lNTime x 3
                     vector<double> tp,                 // m_lNTime x 1
                     vector<double> b0,                 // m_lNPos x 1
                     vector<vector<double>> pr,         // m_lNPos x 3
                     double T1, double T2,
                     vector<vector<complex<double>>> sens, // m_lNPos x m_lNCoils
                     vector<vector<double>> m0)         // m_lNPos x 3
{
    if(prep(b1, gr, tp, b0, pr, sens, m0) == false)
        return false;

    // =================== Create Eigen Vars ===================
    Eigen::MatrixXcd e_b1(m_lNTime, m_lNCoils);
    Eigen::MatrixXd e_gr(m_lNTime, 3);
    for (size_t i = 0; i < m_lNTime; i++)
    {
        e_b1.row(i) = Eigen::Map<Eigen::VectorXcd>(&b1[i][0], m_lNCoils);
        e_gr.row(i) = Eigen::Map<Eigen::VectorXd>(&gr[i][0], 3);
    }

    Eigen::VectorXd e_tp = Eigen::Map<Eigen::VectorXd>(&tp[0], m_lNTime);
    Eigen::VectorXd e_b0 = Eigen::Map<Eigen::VectorXd>(&b0[0], m_lNPos);

    Eigen::MatrixXcd e_sens(m_lNPos, m_lNCoils);
    Eigen::MatrixXd  e_pr(m_lNPos, 3);
    Eigen::MatrixXd  e_m0(m_lNPos, 3);
    for (size_t i = 0; i < m_lNPos; i++)
    {
        e_sens.row(i) = Eigen::Map<Eigen::VectorXcd>(&sens[i][0], m_lNCoils);
        e_pr.row(i)   = Eigen::Map<Eigen::VectorXd>(&pr[i][0], 3);
        e_m0.row(i)   = Eigen::Map<Eigen::VectorXd>(&m0[i][0], 3);
    }

    // =================== Do The Simulation! ===================
    return(runkernel(e_b1, e_gr, e_tp, e_b0, e_sens, e_pr, e_m0, T1, T2));
}

// ----------------------------------------------- //

bool bloch_sim::run(Eigen::MatrixXcd b1,   // m_lNTime x m_lNCoils
                    Eigen::MatrixXd gr,    // m_lNTime x 3
                    Eigen::VectorXd tp,    // second  m_lNTime x 1
                    Eigen::VectorXd b0,    // m_lNPos x 1
                    Eigen::MatrixXd pr,    // m_lNPos x 3
                    double T1, double T2,  // second
                    Eigen::MatrixXcd sens, // m_lNPos x m_lNCoils
                    Eigen::MatrixXd m0)
{
    m_lNCoils = b1.cols();
    m_lNPos   = pr.rows();
    m_lNTime  = b1.rows();

    //std::cout<< "Pos Time Coil = " << m_lNPos << " " << m_lNTime << " " << m_lNCoils << std::endl;
    // =================== Do The Simulation! ===================
    runkernel(b1, gr, tp, b0, sens, pr, m0, T1, T2);
    return true;
}

/*
bool bloch_sim::fastrun(MatrixXd rotx, VectorXd roty, VectorXd rotz, VectorXd e_e1, VectorXd e_e2, MatrixXd e_m0t)
{
    clock_t t = clock();
    concurrency::parallel_for (size_t(0), m_lNPos, [&](size_t cpos){
    //for (size_t cpos=0; cpos<m_lNPos; cpos++){
        Eigen::Matrix3d rotmat;
        Eigen::VectorXd m1;
        for (size_t ct=0; ct<m_lNTime; ct++)
        {
            //std::cout<< "cpos:\n" << cpos << " " << ct << " " << m_lNTime<<std::endl;
            calcrotmat(rotx(ct,cpos), roty(ct,cpos), rotz(ct,cpos), rotmat);
            // m0: Current magnetization before rotation.
            m1 = rotmat * e_m0t.col(cpos);
            // Decay
            m_result_x(cpos, ct) = e_m0t(0,cpos) = m1(0) * e_e2(ct);
            m_result_y(cpos, ct) = e_m0t(1,cpos) = m1(1) * e_e2(ct);
            m_result_z(cpos, ct) = e_m0t(2,cpos) = m1(2) * e_e1(ct) + 1.0 - e_e1(ct);
            //Eigen::VectorXd::Map(&m_result[cpos][ct][0], 3) = e_m0t.col(cpos);
        }
    });
    t = clock() - t;
    std::cout<< "Simulation time " << (float)t/CLOCKS_PER_SEC << " second" << std::endl;
}
*/

// ----------------------------------------------- //
bool bloch_sim::getMagnetization(vector<vector<double>> &result)
{

    const double* begin_x = &m_result_x.col(m_result_x.cols()-1).data()[0];
    result.push_back(std::vector<double>(begin_x, begin_x + m_result_x.rows()));
    const double* begin_y = &m_result_y.col(m_result_y.cols()-1).data()[0];
    result.push_back(std::vector<double>(begin_y, begin_y + m_result_y.rows()));
    const double* begin_z = &m_result_z.col(m_result_z.cols()-1).data()[0];
    result.push_back(std::vector<double>(begin_z, begin_z + m_result_z.rows()));

    return true;
}

bool bloch_sim::getMagnetizationAll(vector<vector<double>> &result_x, vector<vector<double>> &result_y, vector<vector<double>> &result_z)
{
    for (int i=0; i<m_result_x.rows(); i++)
    {
        const double* begin_x = &m_result_x.row(i).data()[0];
        result_x.push_back(std::vector<double>(begin_x, begin_x + m_result_x.cols()));
        const double* begin_y = &m_result_y.row(i).data()[0];
        result_y.push_back(std::vector<double>(begin_y, begin_y + m_result_y.cols()));
        const double* begin_z = &m_result_z.row(i).data()[0];
        result_z.push_back(std::vector<double>(begin_z, begin_z + m_result_z.cols()));
    }
    return true;
}

bool bloch_sim::getMagnetization(Eigen::MatrixXd &result)
{
    result.resize(m_lNPos, 3);
    result.col(0) = m_result_x.col(m_result_x.cols()-1);
    result.col(1) = m_result_y.col(m_result_y.cols()-1);
    result.col(2) = m_result_z.col(m_result_z.cols()-1);
    return true;
}





