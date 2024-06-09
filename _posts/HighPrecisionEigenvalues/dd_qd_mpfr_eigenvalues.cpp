#include <stdio.h>
#include <iostream>
#include <string>
#include <sstream>
#include <Eigen/Eigenvalues>
#include <mpreal.h>
#include <algorithm>
#include <complex>
#include <qd/qd_real.h>
#include <qd/dd_math.h>

using namespace Eigen;
using namespace std;

// define template structure for precision
template<typename T>
struct Precision;

template<>
struct Precision<dd_real> {
    static constexpr int value = 32; // DoubleDouble-precision
};

template<>
struct Precision<qd_real> {
    static constexpr int value = 64; // QuadDouble-precision
};

template<>
struct Precision<mpfr::mpreal> {
    static constexpr int value = 32; // mpreal at 32 digit precision
};

template <typename Real>
void CalcEigenvals(int mat_size) {
    srand(42);
    // Define matrix and vector types used in program. Xmp ending corresponds to mpreal, Xz endings correspond to dd_Real
    typedef Eigen::Matrix<Real,Dynamic,Dynamic,RowMajor>  MatrixXz;
    typedef Eigen::Matrix<Real,Dynamic,1>  VectorXz;

    // generate random square matrices
    MatrixXz lhs_mat = MatrixXz::Random(mat_size, mat_size);

    srand(26);
    MatrixXz rhs_mat = MatrixXz::Random(mat_size, mat_size);

    // calculate eigenvalues using Eigen
    GeneralizedSelfAdjointEigenSolver<MatrixXz> es(lhs_mat, rhs_mat, EigenvaluesOnly);
    VectorXz EigenvalueVec = es.eigenvalues();

    // extract first eigenvalue from vector
    Real ExtractedEigenValue = EigenvalueVec[0];

    cout << "\n" << setprecision(Precision<Real>::value) << ExtractedEigenValue;
}
int main(int argc, char* argv[])
{	
    int mat_size = atoi(argv[1]);

	CalcEigenvals<dd_real>(mat_size);

    CalcEigenvals<qd_real>(mat_size);

    // Specify the global precision to be used (in bits): 
    // 16bit = half, 
    // 32bit = single, 
    // 64bit = double, 
    // 128bit = quadruple
    mpfr::mpreal::set_default_prec(128); 

    CalcEigenvals<mpfr::mpreal>(mat_size);

	return 0;
}