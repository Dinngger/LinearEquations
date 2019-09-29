/*
 * LinearEquationsTest.cpp
 * 
 * Created Date: Saturday, September 28th 2019, 12:28:06 am
 * by Dinger
 */

#include <iostream>
#include <string>
#include <sstream>
#include "Matrix.hpp"
#include "GaussEliminationWithPivoting.hpp"
#include "JacobiIteration.hpp"

int main()
{
    Matrix<double> Ab(3, 4);
    std::cout << "The Matrix A|b is \n";
    std::string sMatrix = "10 -1 -2 72\n"
                          "-1 10 -2 83\n"
                          "-1 -1  5 42\n";
    std::cout << sMatrix;
    std::stringstream ss;
    ss << sMatrix;
    for (int i=0; i<3; i++) {
        for (int j=0; j<4; j++) {
            ss >> Ab.access(i, j);
        }
    }
    GaussEliminationWithPivoting<double> GEP(Ab);
    JacobiIteration<double> JI(Ab);
    std::vector<double> result(3);
    JI.solve(result);
    std::cout << "Jacobi result: \n";
    for (int i=0; i<3; i++) {
        std::cout << result[i] << std::endl;
    }
    GEP.solve(result);
    std::cout << "Gauss result: \n";
    for (int i=0; i<3; i++) {
        std::cout << result[i] << std::endl;
    }
    return 0;
}