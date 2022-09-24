#ifndef _BLOCH_
#define _BLOCH_

#include <complex>

#ifdef WIN32
#ifdef __EXPORT_CLASS_BLOCH__
#define DllExport __declspec(dllexport)
#else
#define DllExport __declspec(dllimport)
#endif
#else
#define DllExport
#endif


#ifdef __SINGLE_PRECISION__
typedef float _T;
#else
typedef double _T;
#endif

extern "C"
{
    bool bloch_sim(std::complex<_T> *pB1,   // RF; m_lNTime x m_lNCoils [Volt] : column-major order {t0c0, t1c0, t2c0,...,t0c1, t1c1, t2c1,...}
                   _T *pGr,                 // gradients; 3 x m_lNTime [Tesla/m] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
                   _T td,                   // dwell-time; [second]
                   _T *pB0,                 // off-resonance; m_lNPos x 1  [Tesla]
                   _T *pPos,                // spatial positions; 3 x m_lNPos  [meter] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
                   std::complex<_T> *pSens, // TX coils sensitivity; m_lNCoils x m_lNPos [Tesla/Volt]: column-major order {c0p0, c1p0, c1p0,...,c0p1, c1p1, c1p1,...}
                   _T T1, _T T2,            // initial magnetization; relaxations; [second]
                   _T *pM0,                 // output; 3 x m_lNPos : column-maj
                   long nPosition,
                   long nTime,
                   long nCoil,
                   _T *pResult,             // 3 x (m_lNTime+1) x m_lNPos : column-major order {x1t0,y1t0,z1t0,...,x1tn,y1tn,z1tn,x2t0,y2t0,z2t0,...}, result equals m0 at t0
                   bool saveAll,            // return all time-points or only the final magnetization
                   bool isPTx);
} // extern

class DllExport bloch
{
public:
    bloch(long nPosition, long nTime, long nCoil, bool saveAll, bool isPTx = true);
    ~bloch();

    bool run(std::complex<_T> *pB1,   // RF; m_lNTime x m_lNCoil [Volt]: column-major order {t0c0, t1c0, t2c0,...,t0c1, t1c1, t2c1,...}
             _T *pGr,                 // gradients; 3 x m_lNTime [Tesla/m] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
             _T td,                   // dwell-time; [second]
             _T *pB0,                 // off-resonance; m_lNPos x 1  [Tesla]
             _T *pPos,                // spatial positions; 3 x m_lNPos  [meter] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
             std::complex<_T> *pSens, // TX coils sensitivity; m_lNCoil x m_lNPos [Tesla/Volt]: column-major order {c0p0, c1p0, c2p0,...,c0p1, c1p1, c2p1,...}
             _T T1, _T T2,            // relaxations; [second]
             _T *pM0,                 // initial magnetization; 3 x m_lNPos : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
             _T *pResult);            // output; 3 x (m_lNTime+1) x m_lNPos : column-major order {x1t0,y1t0,z1t0,...,x1tn,y1tn,z1tn,x2t0,y2t0,z2t0,...}, result equals m0 at t0

protected:
    void timekernel(std::complex<_T> *b1xy, 
                    _T *gr,
                    _T *pr, 
                    _T b0, 
                    _T td_gamma, 
                    _T *m0,
                    _T e1, _T e2, 
                    _T *pResult);   

private:
    std::complex<_T> *m_pB1combined; // combined b1
    int m_lNTime;	// Number of time points
    int m_lNPos;    // Number of positions
    int m_lNCoil;   // Number of Coil
    int m_lStepPos, m_lStepTime, m_lStepB1;
};

#endif // _BLOCH_