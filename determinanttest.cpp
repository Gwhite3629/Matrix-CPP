#include "matrix.cpp"
#include "vector.cpp"

int main(void)
{
    Matrix<double> M(3,3);

    M(0,0) = 1;M(0,1) = 2;M(0,2) = 3;
    M(1,0) = 3;M(1,1) = 2;M(1,2) = 1;
    M(2,0) = 2;M(2,1) = 1;M(2,2) = 3;

    double D = M.laplace_determinant();

    std::cout << "Laplace Determinant:" << D << std::endl;

    D = M.echelon_determinant();

    std::cout << "Echelon Determinant:" << D << std::endl;

    return 0;
}