#include <vector>
#include <complex>
#include <iostream>
#include "./eigen-3.3.9/Eigen/Dense"

using namespace std;
using namespace Eigen;


class bloch_sim
{
public:
    bloch_sim();
    bool run(vector<vector<complex<double>>> b1,// m_lNTime x m_lNCoils
             vector<vector<double>> gr,         // m_lNTime x 3
             vector<double> tp,                 // second  m_lNTime x 1
             vector<double> b0,                 // m_lNPos x 1
             vector<vector<double>> pr,         // m_lNPos x 3
             double T1, double T2,              // second
             vector<vector<complex<double>>> sens = vector<vector<complex<double>>>(), // m_lNPos x m_lNCoils
             vector<vector<double>> m0 = vector<vector<double>>()); // m_lNPos x 3

    bool run(Eigen::MatrixXcd b1,   // m_lNTime x m_lNCoils
             Eigen::MatrixXd gr,    // m_lNTime x 3
             Eigen::VectorXd tp,    // second  m_lNTime x 1
             Eigen::VectorXd b0,    // m_lNPos x 1
             Eigen::MatrixXd pr,    // m_lNPos x 3
             double T1, double T2,  // second
             Eigen::MatrixXcd sens, // m_lNPos x m_lNCoils
             Eigen::MatrixXd m0);

    bool getMagnetization(vector<vector<double>> &result); // 3 x m_lNPos ; is not [m_lNPos x 3] for copying faster
    bool getMagnetizationAll(vector<vector<double>> &result_x, vector<vector<double>> &result_y, vector<vector<double>> &result_z); // m_lNPos x m_lNTime
    bool getMagnetization(Eigen::MatrixXd &result); // m_lNPos x 3;  use this method where possible


    vector<vector<vector<double>>> getMagnetizationAll();

protected:
    bool prep(vector<vector<complex<double>>> &b1,
              vector<vector<double>> &gr,
              vector<double> &tp,
              vector<double> &b0,
              vector<vector<double>> &pr,
              vector<vector<complex<double>>> &sens,
              vector<vector<double>> &m0);

    bool runkernel(Eigen::MatrixXcd &e_b1,
                   Eigen::MatrixXd  &e_gr,
                   Eigen::VectorXd  &e_tp,
                   Eigen::VectorXd  &e_b0,
                   Eigen::MatrixXcd &e_sens,
                   Eigen::MatrixXd  &e_pr,
                   Eigen::MatrixXd  &e_m0,
                   double T1, double T2);

    // m_lNPos is in the first dim as loop is over it
    vector<vector<vector<double>>> m_result;
    Eigen::MatrixXd m_result_x, m_result_y, m_result_z; // m_lNPos x m_lNTime
    size_t m_lNTime;	/* Number of time points. 	 */
    size_t m_lNPos;     /* Number of positions.  Calculated from nposN and nposM, depends on them. */
    size_t m_lNCoils;   /* Number of coils */

};
