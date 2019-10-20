/*
 * TridiagonalMatrix.hpp
 * 
 * Created Date: Sunday, October 20th 2019, 11:37:58 am
 * by Dinger
 */

#ifndef __TRIDIAGONAL_MATRIX_HPP
#define __TRIDIAGONAL_MATRIX_HPP

#include "Matrix.hpp"

template <typename T>
class TridiagonalMatrix : Matrix<T> {
    TridiagonalMatrix() : Matrix<T>::Matrix() {}
    TridiagonalMatrix(uint rows);
    TridiagonalMatrix(uint rows, uint cols);
    TridiagonalMatrix(const TridiagonalMatrix<T> &a) : _rows(a._rows), _cols(a._cols), _size(a._size) {
        _data.resize(a._data.size());
        std::copy(a._data.begin(), a._data.end(), _data.begin());
    };
    Matrix<T> toMatrix() const;
    ~TridiagonalMatrix() {}
    inline T& access(uint row, uint col) override {
        if (abs((int)col - row) < 2) {
            return _data[col * 2 + row];
        } else {
            std::cout << "Tridiagonal Matrix Wroung: Can't access none-diagonal item.\n";
            return nullptr;
        }
    }
    inline T read(uint row, uint col) const override {
        return abs((int)col - row) < 2 ? _data[col * 2 + row] : 0;
    }
    void swapRow(uint i, uint j) override {
        std::cout << "Tridiagonal Matrix Woung: Can't swap row!\n";
    }
    void swapCol(uint i, uint j) override {
        std::cout << "Tridiagonal Matrix Woung: Can't swap col!\n";
    }
};

template <typename T>
TridiagonalMatrix<T>::TridiagonalMatrix(uint rows) :
    _rows(rows),
    _cols(rows) {
    _size = rows * rows;
    _data.resize(3 * rows - 2);
}

template <typename T>
TridiagonalMatrix<T>::TridiagonalMatrix(uint rows, uint cols) {
    if (rows == cols) {
        TridiagonalMatrix(rows);
    } else {
        std::cout << "Tridiagonal Matrix Woung: rows != cols\n";
        TridiagonalMatrix();
    }
}

template<typename T>
Matrix<T> TridiagonalMatrix<T>::toMatrix() const{
    Matrix<T> tempMatrix(_rows, _cols);
    for (uint i=0; i<_rows; i++) {
        for (uint j=0; j<_cols; j++) {
            tempMatrix.access(i, j) = read(i, j);
        }
    }
    return tempMatrix;
};

#endif // __TRIDIAGONAL_MATRIX_HPP
