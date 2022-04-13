#ifndef _BLOCH_
#define _BLOCH_

#include <complex>

extern "C"
{
    bool bloch_sim(std::complex<double> *pB1,   // m_lNTime x m_lNCoils [Volt] : column-major order {t0c0, t1c0, t2c0,...,t0c1, t1c1, t2c1,...}
                   double *pGr,                 // 3 x m_lNTime [Tesla/m] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
                   double td,                   // [second]
                   double *pB0,                 // m_lNPos x 1  [Tesla]
                   double *pPos,                // 3 x m_lNPos  [meter] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
                   std::complex<double> *pSens, // m_lNCoils x m_lNPos [Tesla/Volt]: column-major order {c0p0, c1p0, c1p0,...,c0p1, c1p1, c1p1,...}
                   double T1, double T2,        // [second]
                   double *pM0,                 // 3 x m_lNPos : column-maj
                   long nPosition,
                   long nTime,
                   long nCoil,
                   double *pResult,             // 3 x (m_lNTime+1) x m_lNPos : column-major order {x1t0,y1t0,z1t0,...,x1tn,y1tn,z1tn,x2t0,y2t0,z2t0,...}, result equals m0 at t0
                   bool saveAll);               // return all time-points or only the final magnetization
} // extern

class bloch
{
public:
    bloch(long nPosition, long nTime, long nCoil, bool saveAll);
    ~bloch();

    bool run(std::complex<double> *pB1,   // m_lNTime x m_lNCoil [Volt]: column-major order {t0c0, t1c0, t2c0,...,t0c1, t1c1, t2c1,...}
             double *pGr,                 // 3 x m_lNTime [Tesla/m] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
             double td,                   // [second]
             double *pB0,                 // m_lNPos x 1  [Tesla]
             double *pPos,                // 3 x m_lNPos  [meter] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
             std::complex<double> *pSens, // m_lNCoil x m_lNPos [Tesla/Volt]: column-major order {c0p0, c1p0, c2p0,...,c0p1, c1p1, c2p1,...}
             double T1, double T2,        // [second]
             double *pM0,                 // 3 x m_lNPos : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
             double *pResult);            // 3 x (m_lNTime+1) x m_lNPos : column-major order {x1t0,y1t0,z1t0,...,x1tn,y1tn,z1tn,x2t0,y2t0,z2t0,...}, result equals m0 at t0

protected:
    void timekernel(std::complex<double> *b1xy, 
                    double *gr,
                    double *pr, 
                    double b0, 
                    double td_gamma, 
                    double *m0,
                    double e1, double e2, 
                    double *pResult);    

private:
    std::complex<double> *m_dB1combined; // combined b1
    int m_lNTime;	// Number of time points
    int m_lNPos;    // Number of positions
    int m_lNCoil;   // Number of Coil
    int m_lStepPos, m_lStepTime;
};

#endif // _BLOCH_