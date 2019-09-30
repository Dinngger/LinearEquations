/*
 * JacobiIteration.hpp
 * 
 * Created Date: Sunday, September 29th 2019, 4:20:34 pm
 * by Dinger
 */

#ifndef __JACOBI_ITERATION_HPP
#define __JACOBI_ITERATION_HPP

#include <vector>
#include "Matrix.hpp"

template <typename T>
class JacobiIteration
{
private:
    Matrix<T> _Ab;
public:
    JacobiIteration() {}
    JacobiIteration(Matrix<T> Ab) : _Ab(Ab) {}
    void setMatrixAb(Matrix<T> Ab) {_Ab = Ab;}
    bool solve(std::vector<T>& result) const;
};

template <typename T>
bool JacobiIteration<T>::solve(std::vector<T>& result) const
{
    if (_Ab.rows() != _Ab.cols() - 1) {
        std::cout << "Wrong input!\n";
        return false;
    }
    int n = _Ab.rows();
    Matrix<T> Bf = _Ab;
    for (int i=0; i<n; i++) {
        if (abs(_Ab.read(i, i)) <= 1e-10) {
            std::cout << "Wrong matrix! Zero element(s) on the principal diagonal.\n";
            return false;
        }
        for (int j=0; j<n+1; j++) {
            Bf.access(i, j) /= _Ab.read(i, i);
            if (j < n) {
                Bf.access(i, j) *= -1;
            }
        }
        Bf.access(i, i) = 0;
    }
    T error, last_error=0;
    do {
        std::vector<T> last_res = result;
        error = 0;
        for (int i=0; i<n; i++) {
            result[i] = 0;
            for (int j=0; j<n; j++) {
                result[i] += Bf.read(i, j) * last_res[j];
            }
            result[i] += Bf.read(i, n);
            error = error < abs(result[i] - last_res[i]) ? abs(result[i] - last_res[i]) : error;
        }
        if (last_error != 0 && error > last_error * 2) {
            std::cout << "Wrong! Non convergence!\n";
            return false;
        }
        last_error = error;
    } while (error > 1e-10);
    return true;
}

#endif // __JACOBI_ITERATION_HPP
