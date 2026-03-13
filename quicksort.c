#include "quicksort.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Normal < comparison for qsort
int compare(const void* a, const void* b) { return (*(int*)a - *(int*)b); }

// Wrapper for sequentially sorting lists
void seqSort(int arr[], int N) { qsort(arr, N, sizeof(int), compare); }

// Sets each pivot to the median of the data in the first processor of each group.
void pivotFinderA(int group_id, int arr[], int localLen, int pivots[]) {
    pivots[group_id] = arr[localLen / 2];
}

// Binary search to find pivot positions
int binarySearch(int arr[], int start, int length, int key) {
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
void mergeData(int* localLists[], int localLen[], int m_buf[], int left_list[], int left_len, int right_list[], int right_len, int id) {
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

void gsortHelper(int* localLists[], int localLen[], int* m_bufs[], int N, int n_threads, int iter) {
    int group_size = n_threads >> iter;
    int group_num = n_threads / group_size;
    int pivots[group_num];
    long partitions[n_threads];

#pragma omp parallel num_threads(n_threads)
    {
        int id = omp_get_thread_num();
        int group_id = id / group_size;
        int group_start = group_id * group_size;
        // a. Select pivot
        if (id % group_size == 0)
            pivotFinderA(group_id, localLists[id], localLen[id], pivots);

#pragma omp barrier
        // b. Divide into smaller and larger
        partitions[id] = binarySearch(localLists[id], 0, localLen[id], pivots[group_id]);

        // c. Split processors into two groups and exchange data pairwise
        int pair_id = group_start + ((id - group_start + (group_size >> 1)) % group_size);

        // Part d. Merge data from pair into sorted list

#pragma omp barrier

        int right_len;
        int left_len;
        int* left_list;
        int* right_list;

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

        int* temp = localLists[id];
        localLists[id] = m_bufs[id];
        m_bufs[id] = temp;
    }
}

void gsort(int arr[], int N, int n_threads) {
    int block_size = N / n_threads;
    int remainder = N % n_threads;
    int* localLists[n_threads];
    int localLen[n_threads];
    int* m_bufs[n_threads];

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
        int* localList = malloc(N * sizeof(int));
        memcpy(localList, arr + start, length * sizeof(int));
        localLists[id] = localList;
        localLen[id] = length;
        seqSort(localList, length);

        int* m_buf = malloc(N * sizeof(int));
        m_bufs[id] = m_buf;
    }

    // Main recursive loop
    for (int i = 0; i < log2(n_threads); i++) {
        gsortHelper(localLists, localLen, m_bufs, N, n_threads, i);
    }

    // Finally copying the local lists into the original array
    int location = 0;
    for (int i = 0; i < n_threads; i++) {
        memcpy(arr + location, localLists[i], localLen[i] * sizeof(int));
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
