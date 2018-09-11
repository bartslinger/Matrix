#include "test_macros.hpp"
#include <matrix/math.hpp>

using namespace matrix;

int test_4x3(void);
int test_4x4(void);

int main()
{
    int ret;

    //ret = test_4x3();
    //if (ret != 0) return ret;

    ret = test_4x4();
    if (ret != 0) return ret;

    return 0;
}

int test_4x3() {
    // Start with an (m x n) A matrix
    float data[12] = {20.f , -10.f , -13.f ,
                      17.f ,  16.f , -18.f ,
                       0.7f,  -0.8f,   0.9f,
                      -1.f ,  -1.1f,  -1.2f};
    Matrix<float, 4, 3> A(data);

    float data_check[9] = {-26.27717641f,  -2.76057058f,  21.4699628f ,
                             0.f        , -18.7144129f ,   5.24358701f,
                             0.f        ,   0.f        ,   -2.60681656f};
    Matrix<float, 3, 3> R_check(data_check);

    QRDecomposition<float, 4, 3> qrd = QRDecomposition<float, 4, 3>(A);

    Matrix<float, 3, 3> R = qrd.getR();
    TEST(isEqual(R, R_check));

    float qtb_check_data[4] = {-3.37935852f,
                       -0.53280017f,
                        2.70502894f,
                        5.91429441f};
    Vector<float, 4> qtb_check(qtb_check_data);

    float b_data[4] = {2.0, 3.0, 4.0, 5.0};
    Vector<float, 4> b(b_data);

    Vector<float, 4> qtbv = qrd.qtb(b);
    TEST(isEqual(qtbv, qtb_check));

    float x_check_data[3] = {-0.69168233f,
                             -0.26227593f,
                             -1.03767522f};
    Vector<float, 3> x_check(x_check_data);

    Vector<float, 3> x = qrd.solve(b);

//    TEST(isEqual(x, x_check));

    return 0;
}

int test_4x4() {
    // Start with an (m x n) A matrix
    float data[16] = { 20.f , -10.f , -13.f ,  21.f ,
                       17.f ,  16.f , -18.f , -14.f ,
                        0.7f,  -0.8f,   0.9f,  -0.5f,
                       -1.f ,  -1.1f,  -1.2f,  -1.3f};
    Matrix<float, 4, 4> A(data);

    float data_check[16] = {-26.27717641f,  -2.76057058f,  21.4699628f ,  -6.96231578f,
                              0.f        , -18.7144129f ,   5.24358701f,  24.1199105f ,
                              0.f        ,   0.f        ,  -2.60681656f,  -1.19525533f,
                              0.f        ,   0.f        ,   0.f        ,  -2.69581922f};

    Matrix<float, 4, 4> R_check(data_check);

    QRDecomposition<float, 4, 4> qrd = QRDecomposition<float, 4, 4>(A);

    Matrix<float, 4, 4> R = qrd.getR();

    //TEST(isEqual(R, R_check));
    float b_data[4] = {2.0, 3.0, 4.0, 5.0};
    Vector<float, 4> b(b_data);
    Vector<float, 4> qtbv = qrd.solve(b);
    printf("qtbv: %1.7f %1.7f %1.7f %1.7f", qtbv(0), qtbv(1), qtbv(2), qtbv(3));
    return 0;
}

/* vim: set et fenc=utf-8 ff=unix sts=0 sw=4 ts=4 : */
