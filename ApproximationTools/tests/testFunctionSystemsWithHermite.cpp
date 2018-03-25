#include <iostream>
#include <unistd.h>
#include <vector>
#include "OrthonormalHermite.h"

void printDetailsOfObject(unsigned int longSleep, unsigned int shortSleep, APPRSDK::OrthonormalHermite<double>& H1)
{
    std::cout<<"Retrieving calculated roots"<<std::endl;
    usleep(shortSleep);
    std::cout<<H1.GetDomain()<<std::endl;
    std::cout<<"END OF ROOTS OUTPUT"<<std::endl;
    usleep(longSleep);

    std::cout<<"Retrieving calculated third degree Hermite function"<<std::endl;
    usleep(shortSleep);
    std::cout<<H1.GetFunctionSystem().col(4)<<std::endl;
    std::cout<<"END OF HERMITE FUNCTION OUTPUT"<<std::endl;
    usleep(longSleep);

    std::cout<<"Retrieving calculated derivative (x)"<<std::endl;
    usleep(shortSleep);
    std::cout<<H1.GetDFunctionSystem()<<std::endl;
    std::cout<<"END OF DERIVATIVE OUTPUT"<<std::endl;
    usleep(longSleep);

    std::cout<<"Retrieving partial derivatives (lambda, t)"<<std::endl;
    usleep(shortSleep);
    std::cout<<H1.GetPartialDerivativesFunctionSystem()<<std::endl;
    std::cout<<"END OF PARTIAL DERIVATIVE OUTPUT"<<std::endl;
    usleep(longSleep);
}

int main()
{
    const unsigned int longSleep = 3000;
    const unsigned int shortSleep = 1000;
    std::vector<double> params;

    std::cout<<"Function system test begun..."<<std::endl;
    APPRSDK::OrthonormalHermite<double> H1(100, 10);
    std::cout<<"OrthonormalHermite object succesfully created"<<std::endl;

    printDetailsOfObject(longSleep, shortSleep, H1);

    std::cout<<"Applying parameters: lambda = 0.5, t = 1"<<std::endl;
    params.push_back(0.5);
    params.push_back(1.0);
    H1.ApplyParameters(params);
    std::cout<<"Paramtere application succesful. Printing results...";
    usleep(shortSleep);

    printDetailsOfObject(longSleep, shortSleep, H1);

    return 0;
}