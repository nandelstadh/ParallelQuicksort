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

// Sets each pivot to the median of the data in the first processor of each
// group. Should be optimized later with a better algorithm maybe
void pivotFinderA(int group_id, int arr[], int localLen, int pivots[]) {
    pivots[group_id] = arr[localLen / 2];
}

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

void mergeData(int* localLists[], int localLen[], int m_buf[], int left_list[], int left_len, int r_buf[], int r_buf_len, int id) {
    localLen[id] = left_len + r_buf_len;

    int a = 0;
    int b = 0;
    int comp = 0;

    for (int i = 0; i < localLen[id]; i++) {
        if (a >= left_len)
            m_buf[i] = r_buf[b++];

        else if (b >= r_buf_len)
            m_buf[i] = left_list[a++];

        else if (left_list[a] <= r_buf[b])
            m_buf[i] = left_list[a++];

        else
            m_buf[i] = r_buf[b++];
    }
}

void gsortHelper(int* localLists[], int localLen[], int* r_bufs[], int* m_bufs[], int N, int n_threads, int iter) {
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

        /* printf("Group_num: %d; Group_size: %d; Group_id: %d; Group_start: %d; Id: %d; Pair_Id: %d\n", group_num, group_size, group_id, group_start, id, pair_id); */
        /* printf("The pivot is %d: \n", pivots[group_id]); */
        /* printf("Here we have partition of id %d and pair_id %d: %ld and %ld\n ", id, pair_id, partitions[id], partitions[pair_id]); */

// This barrier is only for debugging
#pragma omp barrier

        int r_buf_len;
        int left_len;
        int* left_list;
        if (id < pair_id) {
            // We receive lower list
            left_list = localLists[id];
            r_buf_len = partitions[pair_id];
            left_len = partitions[id];
            memcpy(r_bufs[id], (int*)localLists[pair_id], r_buf_len * sizeof(int));
        } else {
            // We receive higher list
            left_list = localLists[id] + partitions[id];
            r_buf_len = localLen[pair_id] - partitions[pair_id];
            left_len = localLen[id] - partitions[id];
            memcpy(r_bufs[id], (int*)localLists[pair_id] + partitions[pair_id], r_buf_len * sizeof(int));
        }

#pragma omp barrier
        // Finally, we merge the data into our merger buffer
        mergeData(localLists, localLen, m_bufs[id], left_list, left_len, r_bufs[id], r_buf_len, id);
        memcpy(localLists[id], m_bufs[id], localLen[id] * sizeof(int));
    }
}

/* void gsort(int arr[], int N, int n_threads) { */
/*     seqSort(arr, N); */
/* } */
void gsort(int arr[], int N, int n_threads) {
    int block_size = N / n_threads;
    int remainder = N % n_threads;
    int* localLists[n_threads];
    int localLen[n_threads];
    int* r_bufs[n_threads];
    int* m_bufs[n_threads];

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

        int* localList = malloc(N * sizeof(int));
        memcpy(localList, arr + start, length * sizeof(int));
        localLists[id] = localList;
        localLen[id] = length;
        seqSort(localList, length);

        int* r_buf = malloc(N * sizeof(int));
        int* m_buf = malloc(N * sizeof(int));
        r_bufs[id] = r_buf;
        m_bufs[id] = m_buf;
    }

    for (int i = 0; i < log2(n_threads); i++) {
        gsortHelper(localLists, localLen, r_bufs, m_bufs, N, n_threads, i);
    }

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
        free(r_bufs[id]);
        free(m_bufs[id]);
    }
}
