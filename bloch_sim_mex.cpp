#include "mex.h"
#include "matrix.h"
#include <iostream>
#include "bloch_sim.h"
#include <chrono>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double tp = 0;       // second  m_lNTime x 1
    double T1, T2;   // second
    long nRow, nCol, nPos, nTime, nCoil;
    std::complex<double> *pb1=NULL, *psens=NULL;
    double* pgr=NULL, *pb0=NULL, *ppr=NULL, *pm0=NULL;

    //std::cout << "Parsing data" << std::endl;
    if (nrhs >= 5)
    {
        // map B1
        if(!mxIsComplex(prhs[0]))
            mexErrMsgTxt("B1 must be complex");
        nTime = mxGetM(prhs[0]); // m_lNTime
        nCoil = mxGetN(prhs[0]); // m_lNCoils
        pb1 = reinterpret_cast<std::complex<double>*>(mxGetComplexDoubles(prhs[0]));

        // map gr
        nRow = mxGetM(prhs[1]); // pos
        nCol = mxGetN(prhs[1]); // 3
        if (nRow != 3)
            mexErrMsgTxt("Expected input for gr is 3xN");
        pgr = mxGetPr(prhs[1]);

        // map tp
        nRow = mxGetM(prhs[2]); // 1
        nCol = mxGetN(prhs[2]); // 1
        if (nCol != 1)
            mexErrMsgTxt("Expected input for tp is Nx1");
        tp= *mxGetPr(prhs[2]);

        // map b0
        nRow = mxGetM(prhs[3]); // pos
        nCol = mxGetN(prhs[3]); // 1
        if (nCol != 1)
            mexErrMsgTxt("Expected input for b0 is Nx1");
        pb0 = mxGetPr(prhs[3]);

        // map pr
        nRow = mxGetM(prhs[4]); // pos
        nPos = mxGetN(prhs[4]); // 3
        if (nRow != 3)
            mexErrMsgTxt("Expected input for pr is 3xN");
        ppr = mxGetPr(prhs[4]);
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
        psens = reinterpret_cast<std::complex<double>*>(mxGetComplexDoubles(prhs[7]));
    }

    if (nrhs >= 9)
    {
        // map m0
        nRow = mxGetM(prhs[8]); // m_lNTime
        nCol = mxGetN(prhs[8]); // m_lNCoils
        if (nRow != 3)
            mexErrMsgTxt("Expected input for m0 is 3xN");
        pm0 = mxGetPr(prhs[8]);
    }


    plhs[0] = mxCreateDoubleMatrix(3,nPos, mxREAL);
    double *result = mxGetPr(plhs[0]);

    bloch_sim sim(nPos, nTime, nCoil);
     auto start = std::chrono::system_clock::now();
    //for(int i=0; i<10; i++)
    if(sim.run(pb1, pgr, tp, pb0, ppr, T1, T2, psens, pm0))
    {
        //mexPrintf("Finished successfully\n");
        sim.getMagnetization(result);
        //std::cout << "Results:\n"<< result << std::endl;
    }
    else
        mexErrMsgTxt("Failed.");
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
    std::cout<< "Simulation Matlab " << elapsed.count() << " millisecond" << std::endl;

}



