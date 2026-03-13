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

    float D = m.echelon_determinant();
    std::cout << "Det: " << D << std::endl;
    D = m.laplace_determinant();
    std::cout << "Det: " << D << std::endl;

    Vector<std::complex<double>> E(m.nrow(), 1);

    E = m.eigenvalues(10000);

    E.print(3);

    std::complex<double> Cinf = m.condition(0);
    std::complex<double> Cone = m.condition(1);
    std::complex<double> Ctwo = m.condition(2);
    
    std::cout << "Condition Number: \n";
    std::cout << "One Norm: " << Cone << "\n";
    std::cout << "Two Norm: " << Ctwo << "\n";
    std::cout << "Inf Norm: " << Cinf << std::endl;

    Matrix<float> L1(100, 100, -10, 10);
    Matrix<float> L2(100, 100, -10, 10);

    Matrix<float> R = L1.outer(L2);

    std::cout << "Sum: " << R.sum() << std::endl;
    float FD;
    FD = R.echelon_determinant();
    std::cout << "Det: " << D << std::endl;
    FD = R.laplace_determinant();
    std::cout << "Det: " << D << std::endl;

    Matrix<float> L(100,100);
    Matrix<float> U(100,100);

    R.LU_fast(&L,&U);

    float Ldet = L.echelon_determinant();
    float Udet = U.echelon_determinant();

    FD = Ldet * Udet;
    std::cout << "Det with LU method: " << D << std::endl;

    Matrix<int> Ident(100, 100);
    Ident.I();
    D = Ident.echelon_determinant();
    std::cout << "Identity Det Echelon: " << D << std::endl;
    D = Ident.laplace_determinant();
    std::cout << "Identity Det Laplace: " << D << std::endl;

    //Cinf = R.condition(0);
    //Cone = R.condition(1);
    //Ctwo = R.condition(2);
    
    //std::cout << "Condition Number: \n";
    //std::cout << "One Norm: " << Cone << "\n";
    //std::cout << "Two Norm: " << Ctwo << "\n";
    //std::cout << "Inf Norm: " << Cinf << std::endl;

    Matrix<float> Re(3,4);

    Re(0,0) =  2;Re(0,1) =  1;Re(0,2) = -1;Re(0,3) =   8;
    Re(1,0) = -3;Re(1,1) = -1;Re(1,2) =  2;Re(1,3) = -11;
    Re(2,0) = -2;Re(2,1) =  1;Re(2,2) =  2;Re(2,3) =  -3;

    Matrix<float> In = Re.copy().slice(0,2,0,2);
    Matrix<float> RR = Re.copy();
    Matrix<float> RE = Re.copy();
    Matrix<float> RRE = Re.copy();

    In.invert();
    RR.reduce();
    RE.echelon();
    RRE.echelon();
    RRE.reduce();

    Re.print(4);
    std::cout << std::endl;
    In.print(4);
    std::cout << std::endl;
    RR.print(4);
    std::cout << std::endl;
    RE.print(4);
    std::cout << std::endl;
    RRE.print(4);

    return 0;
}