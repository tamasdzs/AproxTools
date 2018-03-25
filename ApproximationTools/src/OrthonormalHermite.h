#ifndef __ORTHONORMAL_HERMITE_INCLUDED__
#define __ORTHONORMAL_HERMITE_INCLUDED__

#include <math.h>
#include <iostream>
#include "OrthogonalPolynomialBase.h"

namespace APPRSDK
{
    /*! \brief OrthonormalHermite function system

    The OrthonormalHermite class is derived
    from OrthogonalPolynomialBase and implements the orthonormal
    variant of the classical Hermite orthogonal polynomials. 
    The weight function e^(-x^2/2) is included to each function during
    the generation process, therefore the actual functions will be 
    orthogonal to the weight 1.

    A template parameter T is added to the class so that the user can
    decide the precision of the floating point representation. Please
    note however that these functions are defined over R.
    */
    template <typename T>  
    class OrthonormalHermite: public OrthogonalPolynomialBase<T>
    {
        protected:
            unsigned _degrees;
            T _dilatation;
            T _translation;

            void setDomain();
            void setCristoffelDarboux();
            void setOrtPolynomials();
            void setDFunctionSystem();
            void setPartialDerivativesFunctionSystem();    
        public:
            OrthonormalHermite(unsigned int numberOfValues, unsigned int degrees);
            ~OrthonormalHermite();

            void ApplyParameters(std::vector<T> parameters);
            void GenerateWithCostumDomain(Eigen::Array<T, 1, Eigen::Dynamic> domain, unsigned int deg);
    };

    /*! \brief Constructor

    After calling the constructor of the parent class,
    the constructor of OrthonormalHermite proceeds to call
    appropriate methods to calculate the domain, the Cristoffel-Darboux
    numbers and the orthonormal Hermite functions themselves.
    */
    template <typename T>
    OrthonormalHermite<T>::OrthonormalHermite(unsigned int numberOfValues, unsigned int degrees):
        OrthogonalPolynomialBase<T>(numberOfValues, degrees)
    {
        _degrees = degrees;
        _dilatation = 1;
        _translation = 0;
        setDomain();
        setOrtPolynomials();
        setCristoffelDarboux();
        setDFunctionSystem();
        setPartialDerivativesFunctionSystem();
    }

    /*! \brief Destructor
    */
    template <typename T>
    OrthonormalHermite<T>::~OrthonormalHermite()
    {

    }

    /*! \brief setDomain()
    
    Private method setDomain() calculates the domain over
    which the orthonormal Hermite system is considered.

    The domain's points consist of the roots of the (n+1)th
    Hermite polynomial where n is the number of data points
    specified in the constructor's numberOfValues parameter.

    The algorithm to calculate the roots first creates a tridiagonal
    matrix with the help of the alpha and beta values from the three-term
    recurrence formula. We acquire the roots by finding the eigenvalues 
    of this matrix. For further detail please refer to the documentation.
    */
    template <typename T>
    void OrthonormalHermite<T>::setDomain()
    {
        const int n = this->_domain.cols();
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> recMat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(n, n);

        for (int i = 0; i < n; ++i)
        {
            if (i == 0)
            {
                recMat(0, i+1) = sqrt((T)0.5 * (T)(i+1));
            }
            else if (i > 0 && i < (n-1))
            {
                recMat(i, i+1) = sqrt((T)0.5 * (T)(i+1));
                recMat(i, i-1) = sqrt((T)0.5 * (T)(i));
            }
            else
            {
                recMat(i, n-2) = sqrt((T)0.5*(T)i);
            }
        }

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > eigensolver(recMat);
        this->_domain = eigensolver.eigenvalues();
    }

    /*! \brief void setOrtPolynomails()

    The private method setOrtPolynomials() generates the discrete
    orthonormal Hermite function system over the points contained 
    in this->_domain. To achieve this, the method implements the orthonormal
    form of the recursion forumula for Hermite polynomials.
    */
    template <typename T>
    void OrthonormalHermite<T>::setOrtPolynomials()
    {
        const int dataPoints = this->_domain.cols();

        T firstPartialValue, secondPartialValue;

        Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> phi = Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(dataPoints, dataPoints);
        Eigen::Array<T, 1, Eigen::Dynamic> domainArray = this->_domain.array();
        Eigen::Array<T, 1, Eigen::Dynamic> weight = ((T)(-1)*(domainArray*domainArray)/(T)(2.0)).exp();

        phi.col(0) = weight;
        phi.col(1) = (T)(2)*(domainArray*weight)/(T)(sqrt(2));

        for (int i = 2; i < _degrees; ++i)
        {
            firstPartialValue = (T)1.0/(T)sqrt(2.0*double(i));
            secondPartialValue = (T)1.0/(T)sqrt(2.0*(double(i)-1.0)) * firstPartialValue;
            phi.col(i) = (T)2.0*(domainArray*phi.col(i-1).transpose()*firstPartialValue - ((T)i - (T)1)*phi.col(i-2).transpose()*secondPartialValue);
        }

        phi *= (T)pow(4.0*atan(1.0), -1.0/4.0);
        this->_functionSystem = phi;
    }

    /*! \brief void setCristoffelDarboux()

    The private method setCristoffelDarboux() sets the CD
    numbers for the domain the the orthonormal Hermite system.
    */
    template <typename T>
    void OrthonormalHermite<T>::setCristoffelDarboux()
    {
        this->_lambda = (this->_functionSystem)*(this->_functionSystem).transpose();
    }

    /*! \brief void GenerateWithCostumDomain()

    The public method void GenerateWithCostumDomain() generates
    a discrete orthonormal function system given over the data points
    defined in input parameter domain
    */
    template <typename T>
    void OrthonormalHermite<T>::GenerateWithCostumDomain(Eigen::Array<T, 1, Eigen::Dynamic> domain, unsigned int deg)
    {
        this->_domain = domain.matrix();
        _degrees = deg;
        setOrtPolynomials();
    }

    /*! \brief void applyParameters(std::vector<T> parameters)

    The public method applyParameters() applies the parameters given as the input
    to the function system (and if it is considered, also its derivative). In the case
    of the OrthonormalHermite class this will mean the following. First the domain will be replace
    by an equidistant interval (Gauss-Newton formulas cannot be applied
    with costum parameters), then the affine transforms of the Hermite functions will
    be determined over this interval as given in [1]. Parameters:
    -parameters[0] : the dilatation of the function system (the choice of the length of a unit)
    -parameters[1] : the translation of the function system (the choice of the place of origin)

    */
    template<typename T>
    void OrthonormalHermite<T>::ApplyParameters(const std::vector<T> parameters)
    {
        //TODO: Check parameters sanity (i.e length, type) throw exception in case of mismatch
        _dilatation = parameters[0];
        _translation = parameters[1];
        
        Eigen::Array<T, 1, Eigen::Dynamic> domain = Eigen::Array<T, 1, Eigen::Dynamic>::LinSpaced(this->_domain.cols(), this->_domain(0,0), this->_domain(0, this->_domain.cols()-1));
        this->GenerateWithCostumDomain(_dilatation*(this->_domain.array() + _translation), _degrees);
        this->_functionSystem *= (T)sqrt(_dilatation);
    }

    /*! \brief private setDFunctionSystem calculates the 
    matrix of partial derivatives of the orthonormal function
    system with respect to dilatation and translation
    */
    template <typename T>
    void OrthonormalHermite<T>::setDFunctionSystem()
    {
        this->_dFunctionSystem.col(0) = (-1*(_translation + _dilatation*this->_domain.transpose().array())*this->_functionSystem.col(0).array()).matrix();
        for (uint i = 1; i < _degrees; ++i)
        {
            this->_dFunctionSystem.col(i) = ((T)sqrt(2*_degrees)*this->_functionSystem.col(i-1).array() - ((_translation + _dilatation*this->_domain.transpose().array())*this->_functionSystem.col(0).array())).matrix();
        }
    }

    /*! \brief private setPartialDerivativesFunctionSystem() calculates
    the partial derivatives of the orthonormal Hermite system with regards
    to dilatation and translation.
    */
    template <typename T>
    void OrthonormalHermite<T>::setPartialDerivativesFunctionSystem()
    {
        this->_partialDerivativesFunctionSystem.resize(this->_functionSystem.rows(), _degrees*2);
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dDilat;
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dTrans;

        dDilat.resize(this->_functionSystem.rows(), _degrees);
        dTrans.resize(this->_functionSystem.rows(), _degrees);

        for (uint i = 0; i < _degrees; ++i)
        {
            dDilat.col(i) = ((this->_domain.transpose().array() + _translation)*this->_dFunctionSystem.col(i).array()).matrix();
            //DOT: TODO Ask Peter why the -1 is needed
            dTrans.col(i) = -1*(this->_dFunctionSystem.col(i).array()*_dilatation);
        }
        
        uint dilatInd = 0;
        uint transInd = 0;
        for (uint i = 0; i < 2*_degrees; ++i)
        {
            if (i%2 == 0)
            {
                this->_partialDerivativesFunctionSystem.col(i) = dDilat.col(dilatInd);
                this->_index.insert(std::pair<uint, uint>(i, 0));
                dilatInd++;
            }
            else
            {
                this->_partialDerivativesFunctionSystem.col(i) = dTrans.col(transInd);
                this->_index.insert(std::pair<uint, uint>(i, 1));
                transInd++;
            }
        }
    }
}

#endif