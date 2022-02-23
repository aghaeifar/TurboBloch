
#include <iostream>
#include <chrono>
#include "./CPU/bloch.h"

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#include "matrix.h"
#endif

int main()
{
    std::cout << "call sim(...) to execute a simulation" << std::endl;
    return 0;
}

bool sim(std::complex<float> *pb1,   // m_lNTime x m_lNCoils [Volt] : column-major order {t0c0, t1c0, t2c0,...,t0c1, t1c1, t2c1,...}
         float *pgr,                 // 3 x m_lNTime [Tesla/m] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
         float td,                   // m_lNTime x 1 [second]
         float *pb0,                 // m_lNPos x 1  [Tesla]
         float *ppr,                 // 3 x m_lNPos  [meter] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}                    
         std::complex<float> *psens, // m_lNCoils x m_lNPos [Tesla/Volt]: column-major order {c0p0, c1p0, c1p0,...,c0p1, c1p1, c1p1,...}
         float T1, float T2,        // [second]
         float *pm0,                 // 3 x m_lNPos : column-maj
         long nPosition, 
         long nTime, 
         long nCoil,
         float *presult)            // 3 x m_lNPos  [meter] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}  
{
    bloch bloch_obj(nPosition, nTime, nCoil);
    if(bloch_obj.run(pb1, pgr, td, pb0, ppr, psens, T1, T2, pm0) == false)
        return false;

    bloch_obj.getMagnetization(presult);
    return true;
}


#ifdef MATLAB_MEX_FILE

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    float tp = 0;   // second  m_lNTime x 1
    float T1, T2;   // second
    long nRow, nCol, nPos, nTime, nCoil;
    std::complex<float> *pb1=NULL, *psens=NULL;
    float* pgr=NULL, *pb0=NULL, *ppr=NULL, *pm0=NULL;
    std::stringstream buffer;

    if (nrhs != 9)
        mexErrMsgTxt("Wrong number of inputs.");

    // --------- map B1 ---------
    if(!mxIsComplex(prhs[0]))
        mexErrMsgTxt("B1 must be complex");
    nTime = mxGetM(prhs[0]); // m_lNTime
    nCoil = mxGetN(prhs[0]); // m_lNCoils
    pb1 = reinterpret_cast<std::complex<float>*>(mxGetComplexSingles (prhs[0]));

    // --------- map gr ---------
    if (mxGetM(prhs[1]) != 3 || mxGetN(prhs[1]) != nTime)
    {
        buffer << "Expected input for gr is 3x" <<nTime<<std::endl;
        mexErrMsgTxt(buffer.str().c_str());
    }
    pgr = mxGetSingles(prhs[1]);

    // --------- map tp ---------
    if (mxGetM(prhs[2]) * mxGetN(prhs[2]) != 1)
        mexErrMsgTxt("Expected input for tp is 1x1");
    tp= *mxGetSingles(prhs[2]);

    // --------- map b0 ---------
    if (mxGetM(prhs[3]) * mxGetN(prhs[3]) != mxGetN(prhs[3])+mxGetM(prhs[3])-1)
        mexErrMsgTxt("B0 must be vector.");
    nPos = mxGetM(prhs[3]) * mxGetN(prhs[3]);
    pb0 = mxGetSingles(prhs[3]);

    // --------- map pr ---------
    if (mxGetM(prhs[4]) != 3 || mxGetN(prhs[4]) != nPos)
    {
        buffer << "Expected input for pr is 3x" <<nPos<<std::endl;
        mexErrMsgTxt(buffer.str().c_str());
    }
    ppr = mxGetSingles(prhs[4]);

    if (mxGetM(prhs[5]) == 0) // empty input == []
        T1 = -1.0;
    else
        T1 = *mxGetSingles(prhs[5]);

    if (mxGetM(prhs[6]) == 0)
        T2 = -1.0;
    else
        T2 = *mxGetSingles(prhs[6]);

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
        psens = reinterpret_cast<std::complex<float>*>(mxGetComplexSingles(prhs[7]));
    }

    // --------- m0 ---------
    if (mxGetM(prhs[8]) != 3 || mxGetN(prhs[8]) != nPos)
    {
        buffer << "Expected input for m0 is 3x" <<nPos<<std::endl;
        mexErrMsgTxt(buffer.str().c_str());
    }
    pm0 = mxGetSingles(prhs[8]);

    mwSize info, dims[2];
    dims[0] = 3;
    dims[1] = nPos;
    plhs[0] = mxCreateNumericArray (2, dims, mxSINGLE_CLASS, mxREAL);
    float *presult = mxGetSingles(plhs[0]);

    // auto start = std::chrono::system_clock::now();
    if(sim(pb1, pgr, tp, pb0, ppr, psens, T1, T2, pm0, nPos, nTime, nCoil, presult) == false)
        mexErrMsgTxt("Failed.");

    // auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - start);
    // std::cout<< "Simulation Matlab " << elapsed.count() << " millisecond" << std::endl;

}

#endif // MATLAB_MEX_FILE

