#ifndef __IFUNCTION_SYSTEM_INCLUDED__
#define __IFUNCTION_SYSTEM_INCLUDED__

#include <Eigen/Dense>
#include <vector>

namespace APPRSDK
{
    /*! \brief Interface for function system classes.
    *
    *  The IApproximator interface has to be implemented by all concrete FunctionSystem
    *  classes. The interface has a template parameter:
    *  T : type of the signal and function system used for the approximation.
    */
    template<typename T>
    class IFunctionSystem
    {
        protected:

        public:
            virtual Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> GetFunctionSystem() = 0;
            virtual void ApplyParameters(const std::vector<T> parameters) = 0;
    };
}

#endif