#include "quicksort.h"

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

static double selectPivot(const double arr[], int length) {
    if (length <= 0) return 0.0;
    return arr[length / 2];
}

static double groupMedianPivot(int group_start, int group_size, double** localLists, int* localLen) {
    double medians[group_size];
    for (int i = 0; i < group_size; i++) {
        int id = group_start + i;
        medians[i] = selectPivot(localLists[id], localLen[id]);
    }
    qsort(medians, (size_t)group_size, sizeof(double), compare);
    return (medians[group_size / 2 - 1] + medians[group_size / 2]) / 2.0;
}

// Lower bound search to find the split point around a pivot.
static int findSplit(const double arr[], int length, double key) {
    int left = 0, right = length;
    while (left < right) {
        int mid = left + ((right - left) / 2);
        if (arr[mid] < key)
            left = mid + 1;
        else
            right = mid;
    }
    return left;
}

// Merges two sorted lists into a sorted output buffer.
static void mergeData(double out[], const double left[], int left_len, const double right[], int right_len) {
    int a = 0, b = 0, out_len = left_len + right_len;
    for (int i = 0; i < out_len; i++) {
        if (a >= left_len)
            out[i] = right[b++];
        else if (b >= right_len)
            out[i] = left[a++];
        else if (left[a] <= right[b])
            out[i] = left[a++];
        else
            out[i] = right[b++];
    }
}

/*
 * Per-group recursive step:
 * - runs in its own parallel region with exactly n_threads workers
 * - barriers are group-local (not global)
 */
static void gsortGroup(int n_threads, int group_start,
                       double** localLists, double** nextLists,
                       int* localLen, int* nextLen,
                       int* partitions, double* pivots) {
    if (n_threads == 1) return;

#pragma omp parallel num_threads(n_threads)
    {
        int local_id = omp_get_thread_num();
        int id = group_start + local_id;

        if (local_id == 0) {
            pivots[group_start] = groupMedianPivot(group_start, n_threads, localLists, localLen);
        }

#pragma omp barrier
        partitions[id] = findSplit(localLists[id], localLen[id], pivots[group_start]);

#pragma omp barrier

        int half = n_threads / 2;
        int partner = (local_id < half) ? (id + half) : (id - half);

        int my_split = partitions[id];
        int partner_split = partitions[partner];

        const double* left_part;
        int left_len;
        const double* right_part;
        int right_len;

        if (local_id < half) {
            left_part = localLists[id];
            left_len = my_split;
            right_part = localLists[partner];
            right_len = partner_split;
        } else {
            left_part = localLists[partner] + partner_split;
            left_len = localLen[partner] - partner_split;
            right_part = localLists[id] + my_split;
            right_len = localLen[id] - my_split;
        }

        int merged_len = left_len + right_len;
        double* merged = NULL;
        if (merged_len > 0) {
            merged = (double*)malloc((size_t)merged_len * sizeof(double));
            mergeData(merged, left_part, left_len, right_part, right_len);
        }

        nextLists[id] = merged;
        nextLen[id] = merged_len;

#pragma omp barrier

        free(localLists[id]);
        localLists[id] = nextLists[id];
        localLen[id] = nextLen[id];
    }

    int half = n_threads / 2;

    // Recurse on two subgroups in parallel sections for speed
#pragma omp parallel sections
    {
#pragma omp section
        gsortGroup(half, group_start, localLists, nextLists, localLen, nextLen, partitions, pivots);

#pragma omp section
        gsortGroup(half, group_start + half, localLists, nextLists, localLen, nextLen, partitions, pivots);
    }
}

void gsort(double arr[], int N, int n_threads) {
    omp_set_nested(1);
    // If we only have one thread (or tiny input), sort sequentially.
    if (N <= n_threads || n_threads <= 1) {
        seqSort(arr, N);
        return;
    }

    // Best with power-of-two thread count for partner logic.
    // Minimal fallback: if not power of two, use sequential sort.
    if ((n_threads & (n_threads - 1)) != 0) {
        seqSort(arr, N);
        return;
    }

    int block_size = N / n_threads;
    int remainder = N % n_threads;

    double** localLists = (double**)calloc((size_t)n_threads, sizeof(double*));
    double** nextLists = (double**)calloc((size_t)n_threads, sizeof(double*));
    int* localLen = (int*)calloc((size_t)n_threads, sizeof(int));
    int* nextLen = (int*)calloc((size_t)n_threads, sizeof(int));
    int* partitions = (int*)calloc((size_t)n_threads, sizeof(int));
    double* pivots = (double*)calloc((size_t)n_threads, sizeof(double));

    // Initial local sorts (embarrassingly parallel)
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

        double* chunk = (double*)malloc((size_t)length * sizeof(double));
        memcpy(chunk, arr + start, (size_t)length * sizeof(double));
        seqSort(chunk, length);

        localLists[id] = chunk;
        localLen[id] = length;
    }

    // Group-recursive global sort
    gsortGroup(n_threads, 0, localLists, nextLists, localLen, nextLen, partitions, pivots);

    // Gather back
    int write_pos = 0;
    for (int i = 0; i < n_threads; i++) {
        if (localLen[i] > 0) {
            memcpy(arr + write_pos, localLists[i], (size_t)localLen[i] * sizeof(double));
        }
        write_pos += localLen[i];
        free(localLists[i]);
    }

    free(localLists);
    free(nextLists);
    free(localLen);
    free(nextLen);
    free(partitions);
    free(pivots);
}
