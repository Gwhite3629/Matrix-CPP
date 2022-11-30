#include "../matrix.cpp"
#include "../vector.cpp"

#include "benchmark.h"

int main(int argc, char *argv[])
{
    timespec start = {0, 0};
    timespec stop = {0, 0};

    int test_size = TEST_SIZE;
    double r_start = R_START;
    double r_end = R_END;

    srand(time(NULL));

    if (argc == 4) {
        test_size = atoi(argv[1]);
        r_start = strtod(argv[2], NULL);
        r_end = strtod(argv[3], NULL);
    }

    Matrix<double> M(test_size, test_size, r_start, r_end);

    begin(start);

    double det = M.laplace_determinant();

    end(stop);

    std::cout << "Total time: " << get_time(start, stop) << std::endl;
    std::cout << "Determinant: " << det << std::endl;

    return 0;
}