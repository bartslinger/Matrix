/**
 * @file QRDecomposition.hpp
 *
 * Description.
 *
 * @author Bart Slinger <bartslinger@gmail.com>
 */

#pragma once

#include "math.hpp"

namespace matrix {

template<typename Type, size_t M, size_t N>
class QRDecomposition
{
public:

    QRDecomposition() = default;

    // Constructor from A matrix
    QRDecomposition(Matrix<Type, M, N> A)
    {
        // Copy contentents of matrix A
        memcpy(_data, A._data, sizeof(_data));

    }

    Matrix<Type, N, N> getR() {
        Type r_data[N*N];

        Matrix<Type, N, N> R(r_data);
        return R;
    }

private:
    Type _data[M][N] {};

};


} // namespace matrix

/* vim: set et fenc=utf-8 ff=unix sts=0 sw=4 ts=4 : */
