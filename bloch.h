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
    bloch(size_t n_spatial_position,          // number of spatial positions or spins 
          size_t n_timepoints,                // number of time points 
          bool save_all_timepoints = false, // save final magnetization or whole evolution of spins
          bool is_relaxation_constant = true,
          bool is_dwelltime_constant = true
          );
    ~bloch();

    bool run(const std::complex<_T> *pB1,// RF pulse [T]          ; n_timepoints x n_spatial_position
             const _T *pGr,              // gradients [T/m]       ; 3 x n_timepoints: column-major order {gx1,gy1,gz1,gx2,gy2,gz2,...,gxm,gym,gzm}
             const _T *dt,               // time step [Sec]       ; 1 x n_timepoints
             const _T *pB0,              // off-resonance [T]     ; 1 x n_spatial_position
             const _T *pPos,             // spatial positions [m] ; 3 x n_spatial_position: column-major order {x1,y1,z1,...,xm,ym,zm}
             const _T *T1,               // relaxations T1 [Sec]  ; 1 or 1 x n_spatial_position depends on "is_relaxation_constant"
             const _T *T2,               // relaxations T2 [Sec]  ; 1 or 1 x n_spatial_position depends on "is_relaxation_constant"
             const _T *pM0,              // initial magnetization ; 3 x n_spatial_position : column-major order {x1,y1,z1,x2,y2,z2,...,xm,ym,zm}
             _T *pResult                 // output                ; 3 x n_spatial_position or 3 x (n_timepoints+1) x n_spatial_position, depends on "save_all_timepoints": column-major order {x1t1,y1t1,z1t1,...,x1tn,y1tn,z1tn,x2t1,y2t1,z2t1,...,x2tn,y2tn,z2tn,...,xmt1,ymt1,zmt1,...,xmtn,ymtn,zmtn}, result equals m0 at t0
            );
protected:
    void timekernel(const std::complex<_T> *b1, 
                    const _T *gr,
                    const _T *pr, 
                    const _T *b0, 
                    const _T *m0,
                    const _T T1, 
                    const _T T2, 
                    _T *pResult);   

private:
    _T *m_dt;       // Time step
    size_t m_lNTime;	// Number of time points
    size_t m_lNPos;    // Number of positions
    size_t m_lStepPos, m_lStepTime;
    bool m_bConstantT2T2, m_bConstantDwellTime;
};

#endif // _BLOCH_