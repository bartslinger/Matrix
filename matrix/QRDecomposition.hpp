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
            float s = 0.0f;
            for (size_t i = j; i < M; i++) {
                s += _data[i][j] * _data[i][j];
            }
            s = sqrt(s);
            _d[j] = _data[j][j] > 0 ? -s : s;
            float fak = sqrt(s * (s + abs(_data[j][j])));
            _data[j][j] -= _d[j];
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
    }

    Matrix<Type, N, N> getR() {
        // super hacky right now, need to add another test first to cover other dimensions
        Type r_data[N*N];
        r_data[0] = _d[0];
        r_data[1] = _data[0][1];
        r_data[2] = _data[0][2];
        r_data[3] = 0.f;
        r_data[4] = _d[1];
        r_data[5] = _data[1][2];
        r_data[6] = 0.f;
        r_data[7] = 0.f;
        r_data[8] = _d[2];

        Matrix<Type, N, N> R(r_data);
        return R;
    }

private:
    Type _data[M][N] {};
    Type _d[N] {};

};


} // namespace matrix

/* vim: set et fenc=utf-8 ff=unix sts=0 sw=4 ts=4 : */
