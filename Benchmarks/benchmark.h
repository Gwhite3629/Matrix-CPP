#include <ctime>
#include <iostream>
#include <cstdlib>

#define TEST_SIZE 500

#define R_START 100

#define R_END -100

#define GETRAND(T, start, end) \
    ((T)rand()/(T)RAND_MAX) * (end - start) + start;

#define begin(start) \
    clock_gettime(CLOCK_REALTIME, &start)

#define end(stop) \
    clock_gettime(CLOCK_REALTIME, &stop)

#define get_time(start, end) \
    (double)(stop.tv_sec - start.tv_sec) + (double)((stop.tv_nsec - start.tv_nsec)) * 1.e-9
