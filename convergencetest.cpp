#include "matrix.cpp"
#include "vector.cpp"

#include <iostream>
#include <complex>
#include <fstream>

int main(void)
{
    std::ofstream file;

    file.open("data.txt");

    Matrix<double> M(3, 3);
    M.set(0,0,1);
    M.set(0,1,-5);
    M.set(0,2,-6);

    M.set(1,0,3);
    M.set(1,1,1);
    M.set(1,2,2);

    M.set(2,0,1);
    M.set(2,1,7);
    M.set(2,2,4);

    M.print(5);

    double D = M.laplace_determinant();
    std::cout << "Det: " << D << std::endl;

    std::complex<double> Cinf = M.condition(0);
    std::complex<double> Cone = M.condition(1);
    std::complex<double> Ctwo = M.condition(2);
    
    std::cout << "Condition Number: \n";
    std::cout << "One Norm: " << Cone << "\n";
    std::cout << "Two Norm: " << Ctwo << "\n";
    std::cout << "Inf Norm: " << Cinf << std::endl;

    Vector<std::complex<double>> E(M.rows, 1);
    unsigned int k = 0;
    unsigned int b = 1;
    while(k < 100001) {
        for (unsigned int i = 1; i < 101; i++)
        {
            k = b*i;
            E = M.eigenvalues(k);
            file << k << ",";
            for (unsigned int j = 0; j < E.length; j++) {
                file << E.get(j).real();
                if (E.get(j).imag() > 0)
                    file << "+";
                file << E.get(j).imag();
                if (abs(E.get(j).imag()) > 0)
                file << "j";
                if (j < E.length-1)
                    file << ",";
            }
            file << "\n";
        }
        b = b*10;
    }
    file.close();
    return 0;
}