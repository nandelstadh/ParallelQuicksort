#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "quicksort.h"

#define NUM_THREADS 8

static double get_wall_seconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
    return seconds;
}

static int count_values(const double* list, int n, int x) {
    int count = 0;
    int i;
    for (i = 0; i < n; i++) {
        if (list[i] == x) count++;
    }
    return count;
}

double rand_normal() {
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached) {
        double x, y, r;
        do {
            x = (rand() % 2) - 1;
            y = (rand() % 2) - 1;
            r = x * x + y * y;
        } while (r == 0.0 || r > 1.0);
        double d = sqrt(-2.0 * log(r) / r);
        double n1 = x * d, n2 = y * d;
        n2_cached = 1;
        return n1;
    } else {
        n2_cached = 0;
        return n2;
    }
}

double rand_expo(double lambda) {
    double u;
    u = rand() / (RAND_MAX + 1.0);
    return -log(1 - u) / lambda;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        printf("Please give 1 argument: N (number of elements to sort).\n");
        return -1;
    }
    int N = atoi(argv[1]);
    printf("N = %d\n", N);
    if (N < 1) {
        printf("Error: (N < 1).\n");
        return -1;
    }
    double* list_to_sort = (double*)malloc(N * sizeof(double));

    // Creating uniformly distributed list
    /* for (int i = 0; i < N; i++) list_to_sort[i] = rand() % 100; */

    // Creating normally distributed list
    /* for (int i = 0; i < N; i++) list_to_sort[i] = rand_normal(); */

    // Creating exponentially distributed list
    for (int i = 0; i < N; i++) list_to_sort[i] = rand_expo(1);

    int count7 = count_values(list_to_sort, N, 7);
    printf("Before sort: the number 7 occurs %d times in the list.\n", count7);

    // Sort list
    double time1 = get_wall_seconds();
    gsort(list_to_sort, N, NUM_THREADS);
    printf("Sorting list with length %d took %7.3f wall seconds.\n", N,
           get_wall_seconds() - time1);

    int count7_again = count_values(list_to_sort, N, 7);
    printf("After sort : the number 7 occurs %d times in the list.\n",
           count7_again);

    // Check that list is really sorted
    for (int i = 0; i < N - 1; i++) {
        if (list_to_sort[i] > list_to_sort[i + 1]) {
            printf("Error! List not sorted!\n");
            return -1;
        }
    }
    printf("OK, list is sorted!\n");

    free(list_to_sort);

    return 0;
}
