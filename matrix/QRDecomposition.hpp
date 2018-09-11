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

        for (size_t j = 0; j < N; j++) {
            float normx = 0.0f;
            for (size_t i = j; i < M; i++) {
                normx += _data[i][j] * _data[i][j];
            }
            normx = sqrt(normx);
            float s = _data[j][j] > 0 ? -1.0f : 1.0f;
            float u1 = _data[j][j] - s*normx;
            float w[M] = {};
            w[0] = 1.0f;
            for (size_t i = j+1; i < M; i++) {
                w[i-j] = _data[i][j] / u1;
                _data[i][j] = w[i-j];
            }
            _data[j][j] = s*normx;
            _tau[j] = -s*u1/normx;
            printf("tau: %1.7f\n", _tau[j]);

            //printf("w: %1.3f %1.3f %1.3f %1.3f\n", w[0], w[1], w[2], w[3]);
            //printf("diag: %1.7f\n", _data[j][j]);
            for (size_t k = j+1; k < N; k++) {
                float tmp = 0.0f;
                for (size_t i = j; i < M; i++) {
                    tmp += w[i-j] * _data[i][k];
                }
                printf("tmp: %1.7f\n", tmp);
                for (size_t i = j; i < M; i++) {
                    _data[i][k] -= _tau[j] * w[i-j] * tmp;
                    //if (i == 1 && j == 1) {
                    //    printf("JAAA %1.7f\n", _tau[j] * w[i-j] * tmp);
                    //}
                }
            }
            //printf("diag+1: %1.7f\n", _data[j+1][j+1]);

        }
        /*
        return;

        for (size_t j = 0; j < N; j++) {
            float s = 0.0f;
            for (size_t i = j; i < M; i++) {
                s += _data[i][j] * _data[i][j];
            }
            s = sqrt(s);
            _tau[j] = _data[j][j] > 0 ? -s : s;
            float fak = sqrt(s * (s + abs(_data[j][j])));
            _data[j][j] -= _tau[j];
            //printf("d[%zu] = %1.7f\n", j, _d[j]);
            for (size_t k = j; k < M; k++) {
                _data[k][j] /= fak;
            }
            for (size_t i = j+1; i < N; i++) {
                s = 0.0f;
                for (size_t k = j; k < M; k++) {
                    s += _data[k][j] * _data[k][i];
                }
                for (size_t k = j; k < M; k++) {
                    _data[k][i] -= _data[k][j] * s;
                }
            }
        }
        */
    }

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

    /* Find x for Ax = b */
    Vector<Type, N> solve(Vector<Type, M> b) {
        Vector<Type, M> qtbv = qtb(b);
        Vector<Type, N> x;

//        printf("N: %zu;\n", size_t(N));
        for (size_t l = N; l > 0 ; l--) {
            size_t i = l - 1;
            x(i) = qtbv(i);
            for (size_t r = i+1; r < N; r++) {
                x(i) -= _data[i][r] * x(r);
            }
            x(i) = x(i) / _data[i][i];
        }
        return x;
    }

    Matrix<Type, N, N> getR() {
        // super hacky right now, need to add another test first to cover other dimensions
        Type r_data[N][N] = {};
        for (size_t i = 0; i < N; i++) {
            for (size_t j = i; j < N; j++) {
                r_data[i][j] = _data[i][j];
            }
            r_data[i][i] = _tau[i];
        }

        Matrix<Type, N, N> R(_data);
        return R;
    }

private:
    Type _data[M][N] {};
    Type _tau[N] {};

};


} // namespace matrix

/* vim: set et fenc=utf-8 ff=unix sts=0 sw=4 ts=4 : */
