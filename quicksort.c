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

// Lower bound search to find the split point around a pivot.
static int findSplit(const double arr[], int length, double key) {
    int left = 0;
    int right = length;
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

// Merges two sorted lists into a sorted output buffer.
static void mergeData(double out[], const double left[], int left_len, const double right[], int right_len) {
    int a = 0;
    int b = 0;
    int out_len = left_len + right_len;
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

double* gsort(double arr[], int N, int n_threads) {
    // If we only have one thread (or tiny input), sort sequentially.
    if (N <= n_threads || n_threads <= 1) {
        seqSort(arr, N);
        return arr;
    }

    int block_size = N / n_threads;
    int remainder = N % n_threads;
    int allocation_failed = 0;

    double** localLists = calloc((size_t)n_threads, sizeof(double*));
    double** nextLists = calloc((size_t)n_threads, sizeof(double*));
    int* localLen = calloc((size_t)n_threads, sizeof(int));
    int* nextLen = calloc((size_t)n_threads, sizeof(int));
    int* partitions = calloc((size_t)n_threads, sizeof(int));
    double* pivots = calloc((size_t)n_threads, sizeof(double));

    if (!localLists || !nextLists || !localLen || !nextLen || !partitions || !pivots) {
        free(localLists);
        free(nextLists);
        free(localLen);
        free(nextLen);
        free(partitions);
        free(pivots);
        seqSort(arr, N);
        return arr;
    }

#pragma omp parallel num_threads(n_threads)
    {
        // Split list into near-equal contiguous chunks and sort each chunk locally.
        int id = omp_get_thread_num();
        int start, length;
        if (id < remainder) {
            start = id * block_size + id;
            length = block_size + 1;
        } else {
            start = id * block_size + remainder;
            length = block_size;
        }

        double* chunk = NULL;
        if (length > 0) {
            chunk = malloc((size_t)length * sizeof(double));
            if (!chunk) {
#pragma omp atomic write
                allocation_failed = 1;
            } else {
                memcpy(chunk, arr + start, (size_t)length * sizeof(double));
                seqSort(chunk, length);
            }
        }

        localLists[id] = chunk;
        localLen[id] = length;

#pragma omp barrier

        for (int size = n_threads; size > 1; size /= 2) {
            int failed_now;
#pragma omp atomic read
            failed_now = allocation_failed;

            int local_id = id % size;
            int group_id = id / size;

            if (!failed_now && local_id == 0) {
                pivots[group_id] = selectPivot(localLists[id], localLen[id]);
            }

#pragma omp barrier

            if (!failed_now) {
                partitions[id] = findSplit(localLists[id], localLen[id], pivots[group_id]);
            }

#pragma omp barrier

            if (!failed_now) {
                int partner = (local_id < size / 2) ? (id + size / 2) : (id - size / 2);
                int my_split = partitions[id];
                int partner_split = partitions[partner];

                const double* left_part;
                int left_len;
                const double* right_part;
                int right_len;

                if (local_id < size / 2) {
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
                    merged = malloc((size_t)merged_len * sizeof(double));
                    if (!merged) {
#pragma omp atomic write
                        allocation_failed = 1;
                    } else {
                        mergeData(merged, left_part, left_len, right_part, right_len);
                    }
                }

                nextLists[id] = merged;
                nextLen[id] = merged_len;
            }

#pragma omp barrier

#pragma omp atomic read
            failed_now = allocation_failed;

            if (!failed_now) {
                free(localLists[id]);
                localLists[id] = nextLists[id];
                localLen[id] = nextLen[id];
            } else {
                free(nextLists[id]);
                nextLists[id] = NULL;
                nextLen[id] = 0;
            }

#pragma omp barrier

            if (failed_now) break;
        }
    }

    if (allocation_failed) {
        for (int i = 0; i < n_threads; i++) free(localLists[i]);
        free(localLists);
        free(nextLists);
        free(localLen);
        free(nextLen);
        free(partitions);
        free(pivots);
        seqSort(arr, N);
        return arr;
    }

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

    return arr;
}
