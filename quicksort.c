#include "quicksort.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Normal < comparison for qsort
int compare(const void* a, const void* b) {
    double da = *(const double*)a;
    double db = *(const double*)b;
    return (da > db) - (da < db);
}

// Wrapper for sequentially sorting lists
void seqSort(double arr[], int N) { qsort(arr, N, sizeof(double), compare); }

// Sets each pivot to the median of the data in the first processor of each group.
void pivotFinderA(int group_id, double* localLists[], int localLen[], double pivots[], int group_size) {
    int id = group_size * group_id;
    pivots[group_id] = localLists[id][localLen[id] / 2];
}

// Mean of medians.
void pivotFinderB(int group_id, double* localLists[], int localLen[], double pivots[], int group_size) {
    double median = 0;
    for (int i = 0; i < group_size; i++) {
        int id = group_size * group_id + i;
        median += localLists[id][localLen[id] / 2];
    }
    pivots[group_id] = median / group_size;
}

// Mean of middlemost medians.
double pivotFinderC(int group_id, double localLists[], int localLen[], int group_size) {
    double medians[group_size];
    for (int i = 0; i < group_size; i++) {
        int id = group_size * group_id + i;
        // The choice of index here is clearly wrong, we need to think about how to do this when the recursive logic works
        medians[i] = localLists[localLen[id] / 2];
    }
    qsort(medians, group_size, sizeof(double), compare);
    return (medians[group_size / 2] + medians[group_size / 2 - 1]) / 2;
}

// Binary search to find pivot positions
int binarySearch(double arr[], int start, int length, double key) {
    int left = start;
    int right = start + length;
    while (left < right) {
        int mid = left + ((right - left) / 2);

        if (arr[mid] < key) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }

    return left;
}

// Merges two sorted lists into a sorted output buffer
void mergeData(double m_buf[], double left_list[], int left_len, double right_list[], int right_len) {
    int a = 0;
    int b = 0;
    int out_len = left_len + right_len;
    for (int i = 0; i < out_len; i++) {
        if (a >= left_len)
            m_buf[i] = right_list[b++];

        else if (b >= right_len)
            m_buf[i] = left_list[a++];

        else if (left_list[a] <= right_list[b])
            m_buf[i] = left_list[a++];

        else
            m_buf[i] = right_list[b++];
    }
}

// Helper function for one iteration inside the active gsort parallel region.
void gsortHelper(double* arr, double* localLists[], int localLen[], int n_threads, int partitions[], int id) {
    int local_id = n_threads % id;
    int group_id = id / n_threads;
    int group_start = group_id * n_threads;

    if (n_threads == 1) return;

    // a. Select pivot
    double pivot;
    if (local_id == 0) {
        pivot = pivotFinderC(group_id, localLists[id], localLen, n_threads);
    }

#pragma omp barrier
    partitions[id] = binarySearch(localLists[id], 0, localLen[id], pivot);
#pragma omp barrier

    int right_len, left_len;
    double* left_list;
    double* right_list;

    // Exchanging information about lengths of lower and upper lists to be exchanged
    if (local_id < n_threads / 2) {
        // We receive lower list
        left_list = localLists[id];
        left_len = partitions[id];
        right_list = localLists[id + n_threads / 2];
        right_len = partitions[id + n_threads / 2];
    } else {
        // We receive higher list
        left_list = localLists[id] + partitions[id];
        left_len = localLen[id] - partitions[id];
        right_list = localLists[id - n_threads / 2] + partitions[id - n_threads / 2];
        right_len = localLen[id - n_threads / 2] - partitions[id - n_threads / 2];
    }

    int new_len = left_len + right_len;

    // Recompute where the data lies in the buffer

    // Merging the data in the buffer
    double* buffer = malloc(new_len * sizeof(double));
    mergeData(buffer, left_list, left_len, right_list, right_len);
    localLists[id] = buffer;

    gsortHelper();
    free(buffer);

#pragma omp barrier
}

double* gsort(double arr[], int N, int n_threads) {
    // If we only have one thread, we sort sequentially
    if (N <= n_threads || n_threads <= 1) {
        seqSort(arr, N);
        return arr;
    }

    int block_size = N / n_threads;
    int remainder = N % n_threads;
    int localOffsetA[n_threads], localOffsetB[n_threads], localLenA[n_threads], localLenB[n_threads];
    int* localOffset = localOffsetA;
    int* nextOffset = localOffsetB;
    int* localLen = localLenA;
    int* nextLen = localLenB;
    int iterations = log2(n_threads);
    int partitions[n_threads];
    double* localLists[n_threads];
    double pivots[n_threads];

    double* temp = malloc(N * sizeof(double));
    double* buffer = temp;

#pragma omp parallel num_threads(n_threads)
    {
        // Splitting lists evenly between threads
        int id = omp_get_thread_num();
        int start, length;
        if (id < remainder) {
            start = id * block_size + id;
            length = block_size + 1;
        } else {
            start = id * block_size + remainder;
            length = block_size;
        }

        // Sort each thread's local segment in-place.
        localOffset[id] = start;
        localLen[id] = length;
        seqSort(arr + start, length);

        // Perform one iteration of group split/merge inside this parallel region.

        // Set the src array to equal the dst array
        /* double* tmpBuf = arr; */
        /* arr = buffer; */
        /* buffer = tmpBuf; */
        /**/
        /* // Change offsets */
        /* int* tmpOffset = localOffset; */
        /* localOffset = nextOffset; */
        /* nextOffset = tmpOffset; */
        /**/
        /* // Change lengths */
        /* int* tmpLen = localLen; */
        /* localLen = nextLen; */
        /* nextLen = tmpLen; */
    }

    if (arr != temp) {
        free(temp);
    }

    return arr;
}
