#ifndef __APPROXSTAT_INCLUDED__
#define __APPROXSTAT_INCLUDED__

#include <Eigen/Dense>

/*! \brief ApproxStat class
 *
 * The ApproxStat class implements functions that return different. Statistical
 * information about the state of an approximation. These include:
 * - PRD
 * - Sum of Fourier factors
 * - Covariance Matrix
 * - Correlation Matrix
*/

namespace APPRSDK
{
    template <typename T>
    class ApproxStat
    {
        protected:

        public:
            double GetPrd(Eigen::Matrix<T, 1, Eigen::Dynamic>& approximation, 
                          Eigen::Matrix<T, 1, Eigen::Dynamic>& signal);
            
            Eigen::MatrixXd GetCovM(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& weights,
                                    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& dPhi,
                                    unsigned int numberOfDataPoints,
                                    unsigned int numberOfParameters);

            Eigen::MatrixXd GetCorM(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& weights,
                                    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& dPhi,
                                    unsigned int numberOfDataPoints,
                                    unsigned int numberOfParameters,
                                    unsigned int functionSystemDegrees);
    };

    /*! \brief double GetPrd
    *
    * Public method GetPrd returns the calculated PRD of the approximation.
    * This is calculated by dividing the square norm of the difference between the
    * approximation and the signal by the square norm of the difference between the signal
    * and its avarege value. 
    */
    template<typename T>
    double ApproxStat::GetPrd(Eigen::Matrix<T, 1, Eigen::Dynamic>& approximation, 
            Eigen::Matrix<T, 1, Eigen::Dynamic>& signal)
    {
        return ((signal - approximation).norm() / (signal.array() - signal.mean()).matrix().norm());
    }
}

#endif