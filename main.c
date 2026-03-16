#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "quicksort.h"

int T;

static double get_wall_seconds() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
    return seconds;
}

double rand_normal() {
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached) {
        double x, y, r;
        do {
            x = 2.0 * (rand() / (RAND_MAX + 1.0)) - 1.0;
            y = 2.0 * (rand() / (RAND_MAX + 1.0)) - 1.0;
            r = x * x + y * y;
        } while (r == 0.0 || r > 1.0);
        double d = sqrt(-2.0 * log(r) / r);
        double n1 = x * d;
        n2 = y * d;
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
    // Verifying correctness of macros and arguments
    if (argc != 4) {
        printf("Please give 2 arguments: N (number of elements to sort), D (u, n, e for unif, normal and exp distributions resp.), T (number of threads).\n");
        return -1;
    }
    int N = atoi(argv[1]);
    printf("N = %d\n", N);
    if (N < 1) {
        printf("Error: (N < 1).\n");
        return -1;
    }

    // Allocating list on size N
    double* list_to_sort = (double*)malloc(N * sizeof(double));

    char D = argv[2][0];
    printf("D = %c\n", D);
    if (D == 'u') {
        // Creating uniformly distributed list
        for (int i = 0; i < N; i++) list_to_sort[i] = rand() % 100;
    }

    else if (D == 'n') {
        // Creating normally distributed list
        for (int i = 0; i < N; i++) list_to_sort[i] = rand_normal();
    }

    else if (D == 'e') {
        // Creating exponentially distributed list
        for (int i = 0; i < N; i++) list_to_sort[i] = rand_expo(1);
    }

    else {
        printf("Invalid distribution. Valid types are u, n, e.\n");
        return -1;
    }

    int T = atoi(argv[3]);
    if ((T & (T - 1)) != 0) {
        printf("T must be a power of 2\n");
        return -1;
    }
    if (T <= 0) {
        printf("T must be nonnegative\n");
        return -1;
    }

    printf("List sorting starts\n");

    // Sort list
    double time1 = get_wall_seconds();
    double* sorted = gsort(list_to_sort, N, T);
    if (sorted != list_to_sort) {
        free(list_to_sort);
        list_to_sort = sorted;
    }
    /* seqSort(list_to_sort, N); */
    printf("Sorting list with length %d took %7.3f wall seconds.\n", N,
           get_wall_seconds() - time1);

    printf("List sorting finishes\n");

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
