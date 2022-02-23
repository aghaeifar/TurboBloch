#include <complex>
#include <iostream>

//#define EIGEN_USE_MKL_ALL

using namespace std;

class bloch
{
public:
    bloch(long nPosition, long nTime, long nCoil = 1);
    ~bloch();

    bool run(std::complex<float> *b1=NULL,   // m_lNTime x m_lNCoil [Volt]: column-major order {t0c0, t1c0, t2c0,...,t0c1, t1c1, t2c1,...}
             float *gr=NULL,                 // 3 x m_lNTime [Tesla/m] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
             float td=1e-3,                  // m_lNTime x 1 [second]
             float *b0=NULL,                 // m_lNPos x 1  [Tesla]
             float *pr=NULL,                 // 3 x m_lNPos  [meter] : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
             std::complex<float> *sens=NULL, // m_lNCoil x m_lNPos [Tesla/Volt]: column-major order {c0p0, c1p0, c2p0,...,c0p1, c1p1, c2p1,...}
             float T1=-1, float T2=-1,      // [second]
             float *m0=NULL);                // 3 x m_lNPos : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}

    bool getMagnetization(float result[]);   // 3 x m_lNPos : column-major order {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}

protected:
    void print(std::complex<float> *b1,   // m_lNTime x m_lNCoil : {t0c0, t1c0, t2c0,...,t0c1, t1c1, t2c1,...}
               float *gr,                 // m_lNTime x 3 [Tesla] : {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}                  // m_lNTime x 1 [second]
               float *b0,                 // m_lNPos x 1  [Tesla]
               float *pr,                 // m_lNPos x 3  [meter] : {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
               std::complex<float> *sens, // m_lNCoil x m_lNPos : {c0p0, c1p0, c1p0,...,c0p1, c1p1, c1p1,...}
               float *m0);

    float *m_dMagnetization; // m_lNTime * m_lNPos
    std::complex<float> *m_b1combined; // combined b1
    int m_lNTime;	/* Number of time points. 	 */
    int m_lNPos;    /* Number of positions.  Calculated from nposN and nposM, depends on them. */
    int m_lNCoil;  /* Number of Coil */

};
