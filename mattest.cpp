#include "matrix.cpp"
#include "vector.cpp"

int main(int argc, char *argv[])
{
    unsigned int r = 4;
    unsigned int c = 3;
    double s = -2;
    double e = 2;

    Matrix<double> M(r, c, s, e);
    Vector<double> rV(3, 1, 2.0, 2.0);
    Vector<double> cV(4, 0, 1.0, 1.0);

    M.print(3);
    std::cout << "\n";
    rV.print(3);
    std::cout << "\n";
    cV.print(3);
    std::cout << "\n";

    M.AddVecTo(rV, 1);
    M.print(3);
    std::cout << "\n";
    M.AddVecTo(cV, 2);
    M.print(3);
    std::cout << "\n";

    Matrix<double> m = M.del_row(3);
    m.print(3);
    std::cout << "\n";

    double D = m.laplace_determinant();
    std::cout << "Det: " << D << std::endl;

    Vector<std::complex<double>> E(m.rows, 1);

    E = m.eigenvalues(10000);

    E.print(3);

    std::complex<double> Cinf = m.condition(0);
    std::complex<double> Cone = m.condition(1);
    std::complex<double> Ctwo = m.condition(2);
    
    std::cout << "Condition Number: \n";
    std::cout << "One Norm: " << Cone << "\n";
    std::cout << "Two Norm: " << Ctwo << "\n";
    std::cout << "Inf Norm: " << Cinf << std::endl;

    Matrix<double> L1(100, 100, -1, 1);
    Matrix<double> L2(100, 100, -1, 1);

    Matrix<double> R = L1.outer(L2);

    std::cout << "Sum: " << R.sum() << std::endl;

    D = R.laplace_determinant();
    std::cout << "Det: " << D << std::endl;

    Cinf = R.condition(0);
    Cone = R.condition(1);
    Ctwo = R.condition(2);
    
    std::cout << "Condition Number: \n";
    std::cout << "One Norm: " << Cone << "\n";
    std::cout << "Two Norm: " << Ctwo << "\n";
    std::cout << "Inf Norm: " << Cinf << std::endl;

    return 0;
}