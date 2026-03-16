#include "matrix.cpp"
#include "vector.cpp"

int main(int argc, char *argv[])
{
        Matrix<double> A(3,3);
        A(0, 0) = 12; A(0, 1) = -51; A(0, 2) = 4;
        A(1, 0) = 6; A(1, 1) = 167; A(1, 2) = -68;
        A(2, 0) = -4; A(2, 1) = 24; A(2, 2) = -41;

        Matrix<double> givens1 = A.Givens(2,0);
        std::cout << "First givens:" << std::endl;
        givens1.print(4);

        A = givens1.outer(A);

        std::cout << "First rotation:" << std::endl;
        A.print(4);

        return 0;
}
