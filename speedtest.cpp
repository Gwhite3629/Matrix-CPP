#include "matrix.cpp"
#include "vector.cpp"

#include <chrono>

int main(int argc, char *argv[])
{
        char *endptr;
        unsigned int matrixSize = 0;

        if (argc < 2) {
            std::cout << "Correct Usage: \n" << argv[0] << " (matrixSize)" << std::endl;
            return -1;
        }

        matrixSize = std::strtod(argv[1], &endptr);

        Matrix<double> A(matrixSize,matrixSize,-100,100);

        Matrix<double> QGramS(matrixSize,matrixSize);
        Matrix<double> RGramS(matrixSize,matrixSize);
        Matrix<double> ErrGramS(matrixSize,matrixSize);
        double cumErrGramS = 0;

        Matrix<double> QGivens(matrixSize,matrixSize);
        Matrix<double> RGivens(matrixSize,matrixSize);
        Matrix<double> ErrGivens(matrixSize,matrixSize);
        double cumErrGivens = 0;

        auto GramSstart = std::chrono::high_resolution_clock::now();
        A.QR(&QGramS,&RGramS);
        auto GramSstop = std::chrono::high_resolution_clock::now();

        auto GramSduration = std::chrono::duration_cast<std::chrono::microseconds>(GramSstop - GramSstart);
        ErrGramS = A - QGramS.outer(RGramS);
        cumErrGramS = ErrGramS.sum();

        

        auto Givensstart = std::chrono::high_resolution_clock::now();
        A.QR_fast(&QGivens,&RGivens);
        auto Givensstop = std::chrono::high_resolution_clock::now();

        auto Givensduration = std::chrono::duration_cast<std::chrono::microseconds>(Givensstop - Givensstart);
        ErrGivens = A - QGivens.outer(RGivens);
        cumErrGivens = ErrGivens.sum();

        std::cout << "Gram-Schmidt QR Time: " << GramSduration.count() << std::endl;
        std::cout << "Gram-Schmidt QR Error: " << cumErrGramS << std::endl;

        std::cout << "Givens QR Time: " << Givensduration.count() << std::endl;
        std::cout << "Givens QR Error: " << cumErrGivens << std::endl;

        return 0;
}
