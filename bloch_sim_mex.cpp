#include "mex.h"
#include "matrix.h"
#include <iostream>
#include "bloch_sim.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    Eigen::MatrixXcd b1;   // m_lNTime x m_lNCoils
    Eigen::MatrixXd gr;    // m_lNTime x 3
    Eigen::VectorXd tp;    // second  m_lNTime x 1
    Eigen::VectorXd b0;    // m_lNPos x 1
    Eigen::MatrixXd pr;    // m_lNPos x 3
    double T1, T2;         // second
    Eigen::MatrixXcd sens; // m_lNPos x m_lNCoils
    Eigen::MatrixXd m0;

    long nRow, nCol;

    //std::cout << "Parsing data" << std::endl;
    if (nrhs >= 5)
    {
        // map B1
        if(!mxIsComplex(prhs[0]))
            mexErrMsgTxt("B1 must be complex");
        nRow = mxGetM(prhs[0]); // m_lNTime
        nCol = mxGetN(prhs[0]); // m_lNCoils
        std::complex<double>* pb1 = reinterpret_cast<std::complex<double>*>(mxGetComplexDoubles(prhs[0]));
        b1 = Eigen::Map<Eigen::MatrixXcd>(pb1, nRow, nCol);

        // map gr
        nRow = mxGetM(prhs[1]); // m_lNTime
        nCol = mxGetN(prhs[1]); // m_lNCoils
        if (nCol != 3)
            mexErrMsgTxt("Expected input for gr is Nx3");
        double* pgr = mxGetPr(prhs[1]);
        gr = Eigen::Map<Eigen::MatrixXd>(pgr, nRow, nCol);

        // map tp
        nRow = mxGetM(prhs[2]); // m_lNTime
        nCol = mxGetN(prhs[2]); // m_lNCoils
        if (nCol != 1)
            mexErrMsgTxt("Expected input for tp is Nx1");
        double* ptp = mxGetPr(prhs[2]);
        tp = Eigen::Map<Eigen::VectorXd>(ptp, nRow);

        // map b0
        nRow = mxGetM(prhs[3]); // m_lNTime
        nCol = mxGetN(prhs[3]); // m_lNCoils
        if (nCol != 1)
            mexErrMsgTxt("Expected input for b0 is Nx1");
        double* pb0 = mxGetPr(prhs[3]);
        b0 = Eigen::Map<Eigen::VectorXd>(pb0, nRow);

        // map pr
        nRow = mxGetM(prhs[4]); // m_lNTime
        nCol = mxGetN(prhs[4]); // m_lNCoils
        if (nCol != 3)
            mexErrMsgTxt("Expected input for pr is Nx3");
        double* ppr = mxGetPr(prhs[4]);
        pr = Eigen::Map<Eigen::MatrixXd>(ppr, nRow, nCol);

    }
    else
        mexErrMsgTxt("Not enough inputs.");


    if (nrhs >= 7)
    {
        T1 = *mxGetPr(prhs[5]);
        T2 = *mxGetPr(prhs[6]);
    } else
    {
        T1 = 999999;
        T2 = 999999;
    }

    if (nrhs >= 8)
    {
        // map sens
        if(!mxIsComplex(prhs[7]))
            mexErrMsgTxt("B1 must be complex");
        nRow = mxGetM(prhs[7]); // m_lNTime
        nCol = mxGetN(prhs[7]); // m_lNCoils
        std::complex<double>* psens = reinterpret_cast<std::complex<double>*>(mxGetComplexDoubles(prhs[7]));
        sens = Eigen::Map<Eigen::MatrixXcd>(psens, nRow, nCol);

    } else
    {
        sens = Eigen::MatrixXcd::Ones(pr.rows(), b1.cols()); // Only real is 1, imag is zero
    }

    if (nrhs >= 9)
    {
        // map m0
        nRow = mxGetM(prhs[8]); // m_lNTime
        nCol = mxGetN(prhs[8]); // m_lNCoils
        if (nCol != 3)
            mexErrMsgTxt("Expected input for m0 is Nx3");
        double* pm0 = mxGetPr(prhs[8]);
        m0 = Eigen::Map<Eigen::MatrixXd>(pm0, nRow, nCol);

    } else
    {
        m0 = Eigen::MatrixXd::Zero(pr.rows(), 3);
        m0.col(2).setOnes();
    }

    
    
    // start maltab using -nojvm then stdout and stderr are visible
//    std::cout << "b1:\n" << b1<< std::endl;
//    std::cout << "gr:\n" << gr<< std::endl;
//    std::cout << "tp:\n" << tp<< std::endl;
//    std::cout << "b0:\n" << b0<< std::endl;
//    std::cout << "pr:\n" << pr<< std::endl;
//    std::cout << "sens:\n" << sens<< std::endl;
//    std::cout << "m0:\n" << m0<< std::endl;
    
    bloch_sim sim;
    if(sim.run(b1, gr, tp, b0, pr, T1, T2, sens, m0))
    {
        //mexPrintf("Finished successfully\n");
        Eigen::MatrixXd result(pr.rows(), 3);
        sim.getMagnetization(result);
        //std::cout << "Results:\n"<< result << std::endl;
    }
    else
        mexErrMsgTxt("Failed.");
}



