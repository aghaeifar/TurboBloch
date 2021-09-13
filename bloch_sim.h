#include <complex>
#include <iostream>

//#define EIGEN_USE_MKL_ALL

using namespace std;

class bloch_sim
{
public:
    bloch_sim(long nPositions, long nTime, long nCoils = 1);
    ~bloch_sim();

    bool run(std::complex<double> *b1=NULL,   // m_lNTime x m_lNCoils [Volt]: {t0c0, t1c0, t2c0,...,t0c1, t1c1, t2c1,...}
             double *gr=NULL,                 // m_lNTime x 3 [Tesla/m] : {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
             double td=1e-3,                  // m_lNTime x 1 [second]
             double *b0=NULL,                 // m_lNPos x 1  [Tesla]
             double *pr=NULL,                 // 3 x m_lNPos  [meter] : {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
             std::complex<double> *sens=NULL, // m_lNCoils x m_lNPos [Tesla/Volt]: {c0p0, c1p0, c2p0,...,c0p1, c1p1, c2p1,...}
             double T1=NULL, double T2=NULL,  // [second]
             double *m0=NULL);                // {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}

    bool getMagnetization(double result[]);   // m_lNPos x 3 : {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}

protected:
    void print(std::complex<double> *b1,   // m_lNTime x m_lNCoils : {t0c0, t1c0, t2c0,...,t0c1, t1c1, t2c1,...}
               double *gr,                 // m_lNTime x 3 [Tesla] : {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}                  // m_lNTime x 1 [second]
               double *b0,                 // m_lNPos x 1  [Tesla]
               double *pr,                 // m_lNPos x 3  [meter] : {x1,y1,z1,x2,y2,z2,x3,y3,z3,...}
               std::complex<double> *sens, // m_lNCoils x m_lNPos : {c0p0, c1p0, c1p0,...,c0p1, c1p1, c1p1,...}
               double *m0);

    double *m_dMagnetization; // m_lNTime * m_lNPos
    std::complex<double> *m_b1combined; // combined b1
    int m_lNTime;	/* Number of time points. 	 */
    int m_lNPos;    /* Number of positions.  Calculated from nposN and nposM, depends on them. */
    int m_lNCoils;  /* Number of coils */

};
