
#include <iostream>
#include <chrono>
#include "../bloch.h"

#include "mex.h"
#include "matrix.h"

// B1 (complex), gr, td, b0, pr, t1, T2, m0, save_all
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    _T td = 0;   // second  m_lNTime x 1
    _T T1a = 10000., T2a=10000. ;
    _T *T1=&T1a, *T2=&T2a;   // second

    long nRow, nCol, nPos, nTime;
    std::complex<_T> *pb1=NULL;
    _T* pgr=NULL, *pb0=NULL, *ppr=NULL, *pm0=NULL;
    std::stringstream buffer;
    
    if (nrhs < 8)
        mexErrMsgTxt("Wrong number of inputs.");
    
    for (int i=1; i<8; i++)
    #ifdef __SINGLE_PRECISION__
        if (mxIsDouble(prhs[i]))
            mexErrMsgTxt("all inputs must be single!");
    #else
        if (mxIsSingle(prhs[i]))
            mexErrMsgTxt("all inputs must be double!");
    #endif

    // --------- map B1 ---------
    if(!mxIsComplex(prhs[0]))
        mexErrMsgTxt("B1 must be complex");

    nTime = mxGetM(prhs[0]); // m_lNTime
    nPos  = mxGetN(prhs[0]); // m_lNPos

    #ifdef __SINGLE_PRECISION__
    pb1 = reinterpret_cast<std::complex<_T>*>(mxGetComplexSingles (prhs[0]));
    #else    
    pb1 = reinterpret_cast<std::complex<_T>*>(mxGetComplexDoubles (prhs[0]));
    #endif
    
    
    // --------- map gr ---------
    if (mxGetM(prhs[1]) != 3 || mxGetN(prhs[1]) != nTime)
    {
        buffer << "Expected input for gr is 3x" <<nTime<<std::endl;
        mexErrMsgTxt(buffer.str().c_str());
    }
    #ifdef __SINGLE_PRECISION__
    pgr = mxGetSingles(prhs[1]);
    #else
    pgr = mxGetDoubles(prhs[1]);
    #endif
    

    // --------- map td ---------
    if (mxGetM(prhs[2]) * mxGetN(prhs[2]) != 1)
        mexErrMsgTxt("Expected input for td is 1x1");
    #ifdef __SINGLE_PRECISION__
    td = *mxGetSingles(prhs[2]);
    #else
    td = *mxGetDoubles(prhs[2]);
    #endif

    // --------- map b0 ---------
    if (mxGetM(prhs[3]) * mxGetN(prhs[3]) != mxGetN(prhs[3])+mxGetM(prhs[3])-1)
        mexErrMsgTxt("B0 must be vector.");
    if(nPos != mxGetM(prhs[3]) * mxGetN(prhs[3]))
    {
        buffer << "B0 must be a vector or a column of size " <<nPos<<" but it is " << mxGetM(prhs[3]) * mxGetN(prhs[3]) <<std::endl;
        mexErrMsgTxt(buffer.str().c_str());
    }
    #ifdef __SINGLE_PRECISION__
    pb0 = mxGetSingles(prhs[3]);
    #else
    pb0 = mxGetDoubles(prhs[3]);
    #endif

    // --------- map pr ---------
    if (mxGetM(prhs[4]) != 3 || mxGetN(prhs[4]) != nPos)
    {
        buffer << "Expected input for pr is 3x" <<nPos<<std::endl;
        mexErrMsgTxt(buffer.str().c_str());
    }
    #ifdef __SINGLE_PRECISION__
    ppr = mxGetSingles(prhs[4]);
    #else
    ppr = mxGetDoubles(prhs[4]);
    #endif
    
    // --------- T1 & T2 ---------
    if (mxGetM(prhs[5]) * mxGetN(prhs[5]) != 1)
        mexErrMsgTxt("Expected input for T1 is 1x1");
    #ifdef __SINGLE_PRECISION__
    T1 = mxGetSingles(prhs[5]);
    #else
    T1 = mxGetDoubles(prhs[5]);
    #endif

    if (mxGetM(prhs[6]) * mxGetN(prhs[6]) != 1)
        mexErrMsgTxt("Expected input for T2 is 1x1");
    #ifdef __SINGLE_PRECISION__
    T2 = mxGetSingles(prhs[6]);
    #else
    T2 = mxGetDoubles(prhs[6]);
    #endif


    bool isT1T2Constant = true;
    // if (mxGetM(prhs[5]) * mxGetN(prhs[5]) == nPos && mxGetM(prhs[6]) * mxGetN(prhs[6]) == nPos)
    //     isT1T2Constant = false;

    // --------- m0 ---------
    if (mxGetM(prhs[7]) != 3 || mxGetN(prhs[7]) != nPos)
    {
        buffer << "Expected input for m0 is 3x" <<nPos<<std::endl;
        mexErrMsgTxt(buffer.str().c_str());
    }
    #ifdef __SINGLE_PRECISION__
    pm0 = mxGetSingles(prhs[7]);
    #else
    pm0 = mxGetDoubles(prhs[7]);
    #endif

    // --------- save all -----------
    bool saveAll = false;
    if(nrhs > 7)
        if (mxIsLogicalScalar(prhs[8]) && mxGetM(prhs[8]) != 0)
            saveAll = *mxGetLogicals(prhs[8]);

    mwSize info, dims[3];
    dims[0] = 3;
    dims[1] = saveAll ? nTime+1 : 1;
    dims[2] = nPos;

    
    plhs[0] = mxCreateNumericArray(3, dims, 
    #ifdef __SINGLE_PRECISION__
            mxSINGLE_CLASS, 
    #else
            mxDOUBLE_CLASS,
    #endif
            mxREAL);
    
    
    _T *presult = 
    #ifdef __SINGLE_PRECISION__
        mxGetSingles(plhs[0]);
    #else
        mxGetDoubles(plhs[0]);
    #endif

    
    // auto start = std::chrono::system_clock::now();
    auto start = std::chrono::system_clock::now();

    bloch bloch_obj(nPos, nTime, saveAll, isT1T2Constant);
    if(bloch_obj.run(pb1, pgr, td, pb0, ppr, T1, T2, pm0, presult) == false)
        mexErrMsgTxt("Failed.");  

    //mxSetDimensions(plhs[0], dims, 3);
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
    std::cout<< "Bloch simulation took " << elapsed.count() << " millisecond" << std::endl;
}

