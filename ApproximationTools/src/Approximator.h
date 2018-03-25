#ifndef __APPROXIMATOR_INCLUDED__
#define __APPROXIMATOR_INCLUDED__

#include "IApproximator.h"

namespace APPRSDK
{

/*! \brief Approximator class
 * 		   The basic take on the approximator class
 *
 * The Approximator class implements the IApproximator interface.
 * Approximator type objects are capable of conducting an approximation
 * with different parameters, and sharing statistical information about it.
 */
template<typename T>
class Approximator : public IApproximator<T>
{
	protected:
		//Declare function pointer type that can point to costfun and Jacobian functions
		typedef T (Approximator::*methodPointer)();
	
		std::vector<T> _signal;
		std::vector<T> _nonLinParams;
		std::vector<T> _linParams;
		std::vector<methodPointer> _costFunctionPointers;

		IApproximationStrategy* _approximationStrategy;
		IFunctionSystem* _functionSystem;
		ApproxStat<T> _statistics;
	public:
		Approximator();
		virtual ~Approximator();

		std::vector<T> getApproximation();
		std::vector<T> getOptimalNonLinParamSet();
		std::vector<T> getOptimalLinParamSet();
		ApproxStat<T> getStatistics();
		T getCostFun();

};
#endif
