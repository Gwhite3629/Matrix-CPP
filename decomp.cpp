#include "matrix.cpp"
#include "vector.cpp"

int main(void)
{
    Matrix<double> M(3,3);

    M(0,0) = 25;M(0,1) = 15;M(0,2) =  5;
    M(1,0) = 15;M(1,1) = 13;M(1,2) = 11;
    M(2,0) =  5;M(2,1) = 11;M(2,2) = 21;

    M.print(2);

    {
        Matrix<double> Q(3, 3);
        Matrix<double> R(3, 3);
        M.QR(&Q, &R);
        std::cout << "Beginning of Gram-Schmidt QR Algorithm" << std::endl;
        Q.print(2, 2);
        std::cout << std::endl;
        R.print(2, 2);
        std::cout << std::endl;
        Q.outer(R).print(2, 2);
        std::cout << std::endl;
    }

    {
        Matrix<double> Q(3, 3);
        Matrix<double> R(3, 3);
        M.QR_fast(&Q, &R);
        std::cout << "Beginning of Givens QR Algorithm" << std::endl;
        Q.print(2, 2);
        std::cout << std::endl;
        R.print(2, 2);
        std::cout << std::endl;
        Q.outer(R).print(2, 2);
        std::cout << std::endl;
    }

    {
        Matrix<double> L(3, 3);
        Matrix<double> U(3, 3);
        M.LU_fast(&L, &U);
        std::cout << "Beginning of LU Algorithm" << std::endl;
        L.print(2, 2);
        std::cout << std::endl;
        U.print(2, 2);
        std::cout << std::endl;
        L.outer(U).print(2, 2);
        std::cout << std::endl;
    }

    return 0;
}