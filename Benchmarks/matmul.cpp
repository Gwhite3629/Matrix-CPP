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
    int flag = 0; // 0 is float

    srand(time(NULL));

    if (argc == 5) {
        test_size = atoi(argv[1]);
        r_start = strtod(argv[2], NULL);
        r_end = strtod(argv[3], NULL);
        flag = atoi(argv[4]);
    }

    if (flag) {
        Matrix<double> M1(test_size, test_size, r_start, r_end);
        Matrix<double> M2(test_size, test_size, r_start, r_end);

        begin(start);

        M1.outer(M2);

        end(stop);
    } else {
        Matrix<float> M1(test_size, test_size, r_start, r_end);
        Matrix<float> M2(test_size, test_size, r_start, r_end);

        begin(start);

        M1.outer(M2);

        end(stop);
    }

    std::cout << "Total time: " << get_time(start, stop) << std::endl;

    return 0;
}