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


class DllExport bloch
{
public:
    bloch(long nPosition, long nTime, bool hasDiffusion=false, bool saveAll=false);
    ~bloch();

    bool run(const std::complex<_T> *pB1,   // RF; m_lNTime x 1 [Volt]: 
             const _T *pGr,                 // gradients; 3 x m_lNTime [Tesla/m] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
             const _T td,                   // dwell-time; [second]
             const _T *pB0,                 // off-resonance; m_lNPos x 1  [Tesla]
             const _T *pPos,                // spatial positions; 3 x m_lNPos  [meter] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
             const _T T1, const _T T2,      // relaxations; [second]
             const _T *pM0,                 // initial magnetization; 3 x m_lNPos : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
             _T *pResult);                  // output; 3 x (m_lNTime+1) x m_lNPos : column-major order {x1t0,y1t0,z1t0,...,x1tn,y1tn,z1tn,x2t0,y2t0,z2t0,...}, result equals m0 at t0

protected:
    void timekernel(const std::complex<_T> *b1xy, 
                    const _T *gr,
                    const _T *pr, 
                    const _T *b0, 
                    const _T td_gamma, 
                    const _T *m0,
                    const _T e1, const _T e2, 
                    _T *pResult);   

private:
    int m_lNTime;	// Number of time points
    int m_lNPos;    // Number of positions
    int m_lStepPos, m_lStepTime;
    bool m_bHasDiffusion;
};

#endif // _BLOCH_