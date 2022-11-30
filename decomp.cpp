#include "matrix.cpp"
#include "vector.cpp"

int main(void)
{
    Matrix<double> M(3,3);

    M(0,0) = 25;M(0,1) = 15;M(0,2) =  5;
    M(1,0) = 15;M(1,1) = 13;M(1,2) = 11;
    M(2,0) =  5;M(2,1) = 11;M(2,2) = 21;

    {
        Matrix<double> Q(3, 3);
        Matrix<double> R(3, 3);
        M.QR(&Q, &R);

        Q.print(2, 2);
        R.print(2, 2);
        Q.outer(R).print(2, 2);
        std::cout << std::endl;
    }

    {
        Matrix<double> Q(3, 3);
        Matrix<double> R(3, 3);
        M.QR_fast(&Q, &R);

        Q.print(2, 2);
        R.print(2, 2);
        Q.outer(R).print(2, 2);
        std::cout << std::endl;
    }

    return 0;
}