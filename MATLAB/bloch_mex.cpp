
#include <iostream>
#include <chrono>
#include "../bloch.h"

#include "mex.h"
#include "matrix.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double tp = 0;   // second  m_lNTime x 1
    double T1, T2;   // second
    long nRow, nCol, nPos, nTime, nCoil;
    std::complex<double> *pb1=NULL, *psens=NULL;
    double* pgr=NULL, *pb0=NULL, *ppr=NULL, *pm0=NULL;
    std::stringstream buffer;
    
    if (nrhs < 9)
        mexErrMsgTxt("Wrong number of inputs.");

    // --------- map B1 ---------
    if(!mxIsComplex(prhs[0]))
        mexErrMsgTxt("B1 must be complex");
    nTime = mxGetM(prhs[0]); // m_lNTime
    nCoil = mxGetN(prhs[0]); // m_lNCoils
    pb1 = reinterpret_cast<std::complex<double>*>(mxGetComplexDoubles (prhs[0]));

    // --------- map gr ---------
    if (mxGetM(prhs[1]) != 3 || mxGetN(prhs[1]) != nTime)
    {
        buffer << "Expected input for gr is 3x" <<nTime<<std::endl;
        mexErrMsgTxt(buffer.str().c_str());
    }
    pgr = mxGetPr(prhs[1]);

    // --------- map tp ---------
    if (mxGetM(prhs[2]) * mxGetN(prhs[2]) != 1)
        mexErrMsgTxt("Expected input for tp is 1x1");
    tp= *mxGetPr(prhs[2]);

    // --------- map b0 ---------
    if (mxGetM(prhs[3]) * mxGetN(prhs[3]) != mxGetN(prhs[3])+mxGetM(prhs[3])-1)
        mexErrMsgTxt("B0 must be vector.");
    nPos = mxGetM(prhs[3]) * mxGetN(prhs[3]);
    pb0 = mxGetPr(prhs[3]);

    // --------- map pr ---------
    if (mxGetM(prhs[4]) != 3 || mxGetN(prhs[4]) != nPos)
    {
        buffer << "Expected input for pr is 3x" <<nPos<<std::endl;
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

    // --------- map sens ---------
    if (mxGetM(prhs[7]) == 0)
    {
        nCoil = 1;
        psens = NULL;
    }
    else 
    {
        if(!mxIsComplex(prhs[7]))
            mexErrMsgTxt("Sensitivity map must be complex");
        if(mxGetM(prhs[7]) != nCoil || mxGetN(prhs[7]) != nPos)
        {
            buffer << "Expected input for sensitivity is "<<nCoil<<"x"<<nPos<<std::endl;
            mexErrMsgTxt(buffer.str().c_str());
        }
        psens = reinterpret_cast<std::complex<double>*>(mxGetComplexDoubles(prhs[7]));
    }

    // --------- m0 ---------
    if (mxGetM(prhs[8]) != 3 || mxGetN(prhs[8]) != nPos)
    {
        buffer << "Expected input for m0 is 3x" <<nPos<<std::endl;
        mexErrMsgTxt(buffer.str().c_str());
    }
    pm0 = mxGetPr(prhs[8]);

    // --------- save all -----------
    bool saveAll = false;
    if(nrhs > 9)
        if (mxIsLogicalScalar(prhs[9]) && mxGetM(prhs[9]) != 0)
            saveAll = *mxGetLogicals(prhs[9]);

    mwSize info, dims[3];
    dims[0] = 3;
    dims[1] = saveAll ? nTime+1 : 1;
    dims[2] = nPos;
    plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
    double *presult = mxGetPr(plhs[0]);
    
    // auto start = std::chrono::system_clock::now();
    auto start = std::chrono::system_clock::now();

    bloch bloch_obj(nPos, nTime, nCoil, saveAll, true);
    
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
    std::cout<< "Memory allocation took " << elapsed.count() << " millisecond" << std::endl;
    start = std::chrono::system_clock::now();

    if(bloch_obj.run(pb1, pgr, tp, pb0, ppr, psens, T1, T2, pm0, presult) == false)
        mexErrMsgTxt("Failed.");  

    //mxSetDimensions(plhs[0], dims, 3);
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
    std::cout<< "Bloch simulation took " << elapsed.count() << " millisecond" << std::endl;

}

