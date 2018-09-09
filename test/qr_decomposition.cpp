#include "test_macros.hpp"
#include <matrix/math.hpp>

using namespace matrix;

int main()
{

    // Start with an (m x n) A matrix
    float data[12] = {20.f , -10.f , -13.f ,
                      17.f ,  16.f , -18.f ,
                       0.7f,  -0.8f,   0.9f,
                      -1.f ,  -1.1f,  -1.2f};
    Matrix<float, 4, 3> A(data);

    float data_check[9] = {-26.27717641f,  -2.76057058f,  21.4699628f ,
                             0.f        , -18.7144129f ,   5.24358701f,
                             0.f        ,   0.f        ,  -2.60681656f};

    Matrix<float, 3, 3> R_check(data_check);

    QRDecomposition<float, 4, 3> qrd = QRDecomposition<float, 4, 3>(A);

    Matrix<float, 3, 3> R = qrd.getR();

    TEST(isEqual(R, R_check));

    return 0;
}

/* vim: set et fenc=utf-8 ff=unix sts=0 sw=4 ts=4 : */
