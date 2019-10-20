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
    uint32_t n = 4;
    Matrix<double> Ab(n, n + 1);
    std::cout << "The Matrix A|b is \n";
    std::string sMatrix = "2 1 0 0 3\n"
                          "1 2 -3 0 -3\n"
                          "0 3 -7 4 -10\n"
                          "0 0 2 5 2\n";
    std::cout << sMatrix;
    std::stringstream ss;
    ss << sMatrix;
    for (int i=0; i<n; i++) {
        for (int j=0; j<n+1; j++) {
            ss >> Ab.access(i, j);
        }
    }
    GaussEliminationWithPivoting<double> GEP(Ab);
    JacobiIteration<double> JI(Ab);
    std::vector<double> result(n);
    if (!JI.solve(result)) {
        std::cout << "Jacobi error!\n";
    } else {
        std::cout << "Jacobi result: \n";
        for (int i=0; i<n; i++) {
            std::cout << result[i] << std::endl;
        }
    }
    if (!GEP.solve(result)) {
        std::cout << "Gauss error!\n";
    } else {
        std::cout << "Gauss result: \n";
        for (int i=0; i<n; i++) {
            std::cout << result[i] << std::endl;
        }
    }
    return 0;
}
