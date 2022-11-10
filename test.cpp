#include "matrix.cpp"
#include "vector.cpp"

#include <iostream>
#include <complex>

int main(int argc, char *argv[])
{
    if (argc < 5) {
        printf("Correct syntax:\n\t./%s <rows> <cols> <start> <end>\n", argv[0]);
        return -1;
    }

    unsigned int rows = atoi(argv[1]);
    unsigned int cols = atoi(argv[2]);
    double start = strtod(argv[3], NULL);
    double end = strtod(argv[4], NULL);

    Matrix<double> M(rows, cols, start, end);

    M.print(5);

    double D = M.laplace_determinant();
    std::cout << "Det: " << D << std::endl;

    Matrix<double> Q(rows, cols);
    Matrix<double> R(rows, cols);

    M.QR(&Q, &R);

    //Q.print(5);
    //R.print(5);
    M = Q.outer(R);
    M.print(5);

    Vector<std::complex<double>> E(M.rows, 1);

    E = M.eigenvalues(10000);

    E.print(5);

    std::complex<double> Cinf = M.condition(0);
    std::complex<double> Cone = M.condition(1);
    std::complex<double> Ctwo = M.condition(2);
    
    std::cout << "Condition Number: \n";
    std::cout << "One Norm: " << Cone << "\n";
    std::cout << "Two Norm: " << Ctwo << "\n";
    std::cout << "Inf Norm: " << Cinf << std::endl;

    return 0;
}