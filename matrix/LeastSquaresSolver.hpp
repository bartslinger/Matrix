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
class LeastSquaresSolver
{
public:

    LeastSquaresSolver() = default;

    /**
     * @brief Class calculates QR decomposition which can be used for linear
     * least squares
     * @param A Matrix of size MxN
     *
     * Initialize the class with a MxN matrix. The constructor starts the
     * QR decomposition. This class does not check the rank of the matrix.
     * The user needs to make sure that rank(A) = N.
     */
    LeastSquaresSolver(Matrix<Type, M, N> A)
    {
        // Copy contentents of matrix A
        memcpy(_data, A._data, sizeof(_data));

        for (size_t j = 0; j < N; j++) {
            float normx = 0.0f;
            for (size_t i = j; i < M; i++) {
                normx += _data[i][j] * _data[i][j];
            }
            normx = sqrt(normx);
            float s = _data[j][j] > 0 ? -1.0f : 1.0f;
            float u1 = _data[j][j] - s*normx;
            // prevent divide by zero
            // also covers u1. normx is never negative
            if (normx < 1e-8f) {
                return;
            }
            float w[M] = {};
            w[0] = 1.0f;
            for (size_t i = j+1; i < M; i++) {
                w[i-j] = _data[i][j] / u1;
                _data[i][j] = w[i-j];
            }
            _data[j][j] = s*normx;
            _tau[j] = -s*u1/normx;

            for (size_t k = j+1; k < N; k++) {
                float tmp = 0.0f;
                for (size_t i = j; i < M; i++) {
                    tmp += w[i-j] * _data[i][k];
                }
                for (size_t i = j; i < M; i++) {
                    _data[i][k] -= _tau[j] * w[i-j] * tmp;
                }
            }

        }
    }

    /**
     * @brief qtb Calculate Q^T * b
     * @param b
     * @return Q^T*b
     *
     * This function calculates Q^T * b. This is useful for the solver
     * because R*x = Q^T*b.
     */
    Vector<Type, M> qtb(Vector<Type, M> b) {
        Vector<Type, M> qtbv = b;

        for (size_t j = 0; j < N; j++) {
            float w[M];
            w[0] = 1.0f;
            // fill vector w
            for (size_t i = j+1; i < M; i++) {
                w[i-j] = _data[i][j];
            }
            float tmp = 0.0f;
            for (size_t i = j; i < M; i++) {
                tmp += w[i-j] * qtbv(i);
            }

            for (size_t i = j; i < M; i++) {
                qtbv(i) -= _tau[j] * w[i-j] * tmp;
            }
        }
        return qtbv;
    }

    /**
     * @brief Solve Ax=b for x
     * @param b
     * @return Vector x
     *
     * Find x in the equation Ax = b.
     * A is provided in the initializer of the class.
     */
    Vector<Type, N> solve(Vector<Type, M> b) {
        Vector<Type, M> qtbv = qtb(b);
        Vector<Type, N> x;

        for (size_t l = N; l > 0 ; l--) {
            size_t i = l - 1;
            x(i) = qtbv(i);
            for (size_t r = i+1; r < N; r++) {
                x(i) -= _data[i][r] * x(r);
            }
            // divide by zero, return vector of zeros
            if (fabs(_data[i][i]) < 1e-8f) {
                for (size_t z = 0; z < N; z++) {
                    x(z) = 0.0f;
                }
            }
            x(i) = x(i) / _data[i][i];
        }
        return x;
    }

private:
    Type _data[M][N] {};
    Type _tau[N] {};

};

} // namespace matrix

/* vim: set et fenc=utf-8 ff=unix sts=0 sw=4 ts=4 : */
