#ifndef __IAPPROXIMATOR_INCLUDED__
#define __IAPPROXIMATOR_INCLUDED__

#include <vector>

namespace APPRSDK
{

/*! \brief Interface for concrete Approximator classes
 *         Contains the methods to be implemented by children Approximator classes.
 *
 *  The IApproximator interface has to be implemented by all concrete Approximator
 *  classes. The interface has a template parameter:
 *  T : type of the signal and function system used for the approximation.
 */
template<typename T>
class IApproximator
{
	protected:

	public:
	IApproximator() {}
	virtual ~IApproximator() {}

	virtual std::vector<T> getApproximation() = 0;
	virtual std::vector<T> getOptimalNonLinParamSet() = 0;
	virtual std::vector<T> getOptimalLinParamSet() = 0;
	virtual ApproxStat<T> getStatistics() = 0;
	virtual T costFun() = 0;
	virtual void applyParameters() = 0;
};

}
#endif
