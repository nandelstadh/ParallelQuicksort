#include "quicksort.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Normal < comparison for qsort
int compare(const void* a, const void* b) { return (*(double*)a - *(double*)b); }

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
void pivotFinderC(int group_id, double* localLists[], int localLen[], double pivots[], int group_size) {
    double medians[group_size];
    for (int i = 0; i < group_size; i++) {
        int id = group_size * group_id + i;
        medians[i] = localLists[id][localLen[id] / 2];
    }
    qsort(medians, group_size, sizeof(double), compare);
    pivots[group_id] = (medians[group_size / 2] + medians[group_size / 2 - 1]) / 2;
}

// Binary search to find pivot positions
int binarySearch(double arr[], int start, int length, int key) {
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

// Merges sorted lists into a bigger sorted list in a buffer
void mergeData(double* localLists[], int localLen[], double m_buf[], double left_list[], int left_len, double right_list[], int right_len, int id) {
    localLen[id] = left_len + right_len;

    int a = 0;
    int b = 0;

    for (int i = 0; i < localLen[id]; i++) {
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

void gsortHelper(double* localLists[], int localLen[], double* m_bufs[], int N, int n_threads, int iter) {
    int group_size = n_threads >> iter;
    int group_num = n_threads / group_size;
    double pivots[group_num];
    long partitions[n_threads];

#pragma omp parallel num_threads(n_threads)
    {
        int id = omp_get_thread_num();
        int group_id = id / group_size;
        int group_start = group_id * group_size;
        // a. Select pivot
        if (id % group_size == 0)
            pivotFinderC(group_id, localLists, localLen, pivots, group_size);

#pragma omp barrier
        // b. Divide into smaller and larger
        partitions[id] = binarySearch(localLists[id], 0, localLen[id], pivots[group_id]);

        // c. Split processors into two groups and exchange data pairwise
        int pair_id = group_start + ((id - group_start + (group_size >> 1)) % group_size);

        // Part d. Merge data from pair into sorted list

#pragma omp barrier

        int right_len;
        int left_len;
        double* left_list;
        double* right_list;

        // Exchanging information about lengths of lower and upper lists to be exchanged
        if (id < pair_id) {
            // We receive lower list
            left_list = localLists[id];
            left_len = partitions[id];
            right_list = localLists[pair_id];
            right_len = partitions[pair_id];
        } else {
            // We receive higher list
            left_list = localLists[id] + partitions[id];
            left_len = localLen[id] - partitions[id];
            right_list = localLists[pair_id] + partitions[pair_id];
            right_len = localLen[pair_id] - partitions[pair_id];
        }

#pragma omp barrier
        // Merging the data into our merger buffer
        mergeData(localLists, localLen, m_bufs[id], left_list, left_len, right_list, right_len, id);

        double* temp = localLists[id];
        localLists[id] = m_bufs[id];
        m_bufs[id] = temp;
    }
}

void gsort(double arr[], int N, int n_threads) {
    int block_size = N / n_threads;
    int remainder = N % n_threads;
    double* localLists[n_threads];
    int localLen[n_threads];
    double* m_bufs[n_threads];

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

        // Splitting and sorting the lists locally
        double* localList = malloc(N * sizeof(double));
        memcpy(localList, arr + start, length * sizeof(double));
        localLists[id] = localList;
        localLen[id] = length;
        seqSort(localList, length);

        double* m_buf = malloc(N * sizeof(double));
        m_bufs[id] = m_buf;
    }

    // Main recursive loop
    for (int i = 0; i < log2(n_threads); i++) {
        gsortHelper(localLists, localLen, m_bufs, N, n_threads, i);
    }

    // Finally copying the local lists into the original array
    int location = 0;
    for (int i = 0; i < n_threads; i++) {
        memcpy(arr + location, localLists[i], localLen[i] * sizeof(double));
        location += localLen[i];
    }

//  Freeing the lists and buffers
#pragma omp parallel num_threads(n_threads)
    {
        int id = omp_get_thread_num();
        free(localLists[id]);
        free(m_bufs[id]);
    }
}
