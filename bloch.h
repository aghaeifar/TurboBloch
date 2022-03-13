#ifndef _BLOCH_
#define _BLOCH_

#include <complex>

extern "C"
{
    bool bloch_sim(std::complex<float> *pB1,   // m_lNTime x m_lNCoils [Volt] : column-major order {t0c0, t1c0, t2c0,...,t0c1, t1c1, t2c1,...}
                   float *pGr,                 // 3 x m_lNTime [Tesla/m] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
                   float td,                   // [second]
                   float *pB0,                 // m_lNPos x 1  [Tesla]
                   float *pPos,                // 3 x m_lNPos  [meter] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
                   std::complex<float> *pSens, // m_lNCoils x m_lNPos [Tesla/Volt]: column-major order {c0p0, c1p0, c1p0,...,c0p1, c1p1, c1p1,...}
                   float T1, float T2,         // [second]
                   float *pM0,                 // 3 x m_lNPos : column-maj
                   long nPosition,
                   long nTime,
                   long nCoil,
                   float *pResult,             // 3 x (m_lNTime+1) x m_lNPos : column-major order {x1t0,y1t0,z1t0,...,x1tn,y1tn,z1tn,x2t0,y2t0,z2t0,...}, result equals m0 at t0
                   bool saveAll);              // return all time-points or only the final magnetization
} // extern

class bloch
{
public:
    bloch(long nPosition, long nTime, long nCoil, bool saveAll);
    ~bloch();

    bool run(std::complex<float> *pB1,   // m_lNTime x m_lNCoil [Volt]: column-major order {t0c0, t1c0, t2c0,...,t0c1, t1c1, t2c1,...}
             float *pGr,                 // 3 x m_lNTime [Tesla/m] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
             float td,                   // [second]
             float *pB0,                 // m_lNPos x 1  [Tesla]
             float *pPos,                // 3 x m_lNPos  [meter] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
             std::complex<float> *pSens, // m_lNCoil x m_lNPos [Tesla/Volt]: column-major order {c0p0, c1p0, c2p0,...,c0p1, c1p1, c2p1,...}
             float T1, float T2,         // [second]
             float *pM0,                 // 3 x m_lNPos : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
             float *pResult);            // 3 x (m_lNTime+1) x m_lNPos : column-major order {x1t0,y1t0,z1t0,...,x1tn,y1tn,z1tn,x2t0,y2t0,z2t0,...}, result equals m0 at t0

protected:
    void timekernel(std::complex<float> *b1xy, 
                    float *gr,
                    float *pr, 
                    float b0, 
                    float td_gamma, 
                    float *m0,
                    float e1, float e2, 
                    float *pResult);    

private:
    std::complex<float> *m_fB1combined; // combined b1
    int m_lNTime;	// Number of time points
    int m_lNPos;    // Number of positions
    int m_lNCoil;   // Number of Coil
    int m_lStepPos, m_lStepTime;
};

#endif // _BLOCH_