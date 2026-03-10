#include "quicksort.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void swap(int* a, int* b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

int compare(const void* a, const void* b) { return (*(int*)a - *(int*)b); }

void seqSort(int arr[], int N) { qsort(arr, N, sizeof(int), compare); }

void gsortHelper(int arr[], int N, int n_threads, int depth) {
    printf("We got into the helper function");
    printf("I know what my ID is, it is: %d \n", omp_get_thread_num());
    // a. Select pivot
#pragma omp barrier
    // b. Divide into smaller and larger
    // c. Split processors into two groups and exchange data pairwise
    // d. Merge data into sorted list
}

void gsort(int arr[], int N, int n_threads) {
    int block_size = N / n_threads;
    int remainder = N % n_threads;
#pragma omp parallel num_threads(n_threads)
    {
        int id = omp_get_thread_num();
        int start, length;
        if (id < remainder) {
            start = id * block_size + id;
            length = block_size + 1;
        } else {
            start = id * block_size + remainder;
            length = block_size;
        }
        int* localList = malloc(length * sizeof(int));
        memcpy(localList, arr + start, length * sizeof(int));
        seqSort(localList, length);

        for (int i = 0; i < log2(n_threads); i++) {
            gsortHelper(localList, length, n_threads, i);
        }
        free(localList);
    }
}
