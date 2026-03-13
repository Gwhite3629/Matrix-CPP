#include "matrix.cpp"
#include "vector.cpp"

int main(int argc, char *argv[])
{
        Matrix<double> A(3,3);
        A(0, 0) = 6; A(0, 1) = 5; A(0, 2) = 0;
        A(1, 0) = 5; A(1, 1) = 1; A(1, 2) = 4;
        A(2, 0) = 0; A(2, 1) = 4; A(2, 2) = 3;

        Matrix<double> givens1 = A.Givens(1,0);
        std::cout << "First givens:" << std::endl;
        givens1.print(4);

        A = givens1.outer(A);

        std::cout << "First rotation:" << std::endl;
        A.print(4);

        Matrix<double> givens2 = A.Givens(2,1);
        std::cout << "Second givens:" << std::endl;
        givens2.print(4);

        A = givens2.outer(A);

        std::cout << "Second rotation:" << std::endl;
        A.print(4);

        return 0;
}
