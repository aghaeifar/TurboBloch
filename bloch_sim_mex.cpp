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
    std::stringstream buffer;

    if (nrhs != 9)
        mexErrMsgTxt("Wrong number of inputs.");

    // map B1
    if(!mxIsComplex(prhs[0]))
        mexErrMsgTxt("B1 must be complex");
    nTime = mxGetM(prhs[0]); // m_lNTime
    nCoil = mxGetN(prhs[0]); // m_lNCoils
    pb1 = reinterpret_cast<std::complex<double>*>(mxGetComplexDoubles(prhs[0]));

    // map gr
    if (mxGetM(prhs[1]) != 3 || mxGetN(prhs[1]) != nTime)
    {
        buffer << "Expected input for gr is 3x" <<nTime<<std::endl;
        mexErrMsgTxt(buffer.str().c_str());
    }
    pgr = mxGetPr(prhs[1]);

    // map tp
    if (mxGetM(prhs[2]) * mxGetN(prhs[2]) != 1)
        mexErrMsgTxt("Expected input for tp is 1x1");
    tp= *mxGetPr(prhs[2]);

    // map b0
    if (mxGetM(prhs[3]) * mxGetN(prhs[3]) != mxGetN(prhs[3])+mxGetM(prhs[3])-1)
        mexErrMsgTxt("B0 must be vector.");
    nPos = mxGetM(prhs[3]) * mxGetN(prhs[3]);
    pb0 = mxGetPr(prhs[3]);

    // map pr
    if (mxGetM(prhs[4]) != 3 || mxGetN(prhs[4]) != nPos)
    {
        buffer << "Expected input for gr is 3x" <<nPos<<std::endl;
        mexErrMsgTxt(buffer.str().c_str());
    }
    ppr = mxGetPr(prhs[4]);

    if (mxGetM(prhs[5]) == 0) // empty input == []
        T1 = -1.0;
    else
        T1 = *mxGetPr(prhs[5]);

    if (mxGetM(prhs[6]) == 0)
        T2 = -1.0;
    else
        T2 = *mxGetPr(prhs[6]);

    // map sens
    if (mxGetM(prhs[7]) == 0)
    {
        nCoil = 1;
        psens = NULL;
    }
    else
    {
        if(!mxIsComplex(prhs[7]))
            mexErrMsgTxt("B1 must be complex");
        if(mxGetM(prhs[7]) != nCoil || mxGetN(prhs[7]) != nPos)
        {
            buffer << "Expected input for sensitivity is "<<nCoil<<"x"<<nPos<<std::endl;
            mexErrMsgTxt(buffer.str().c_str());
        }
        psens = reinterpret_cast<std::complex<double>*>(mxGetComplexDoubles(prhs[7]));
    }

    // m0
    if (mxGetM(prhs[8]) != 3 || mxGetN(prhs[8]) != nPos)
    {
        buffer << "Expected input for m0 is 3x" <<nPos<<std::endl;
        mexErrMsgTxt(buffer.str().c_str());
    }
    pm0 = mxGetPr(prhs[8]);


    plhs[0] = mxCreateDoubleMatrix(3,nPos, mxREAL);
    double *result = mxGetPr(plhs[0]);

    bloch_sim sim(nPos, nTime, nCoil);
    //auto start = std::chrono::system_clock::now();
    //for(int i=0; i<10; i++)
    if(sim.run(pb1, pgr, tp, pb0, ppr, psens, T1, T2, pm0))
    {
        //mexPrintf("Finished successfully\n");
        sim.getMagnetization(result);
        //std::cout << "Results:\n"<< result << std::endl;
    }
    else
        mexErrMsgTxt("Failed.");
    //auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
    //std::cout<< "Simulation Matlab " << elapsed.count() << " millisecond" << std::endl;

}



