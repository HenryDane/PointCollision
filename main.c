#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "files.h"
#include "collisions.h"

void fast_test(size_t n_pts, float* xs, float* vs, float* rs) {
    printf("==================== FAST TEST ====================\n");
    clock_t start = clock();

    size_t n_pairs;
    pair_t *pairs;
    make_collision_pairs(n_pts, xs, vs, rs, &n_pairs, &pairs);

    // count collisions
    size_t n_colls = count_collisions(n_pts, xs, vs, rs, n_pairs, pairs);

    clock_t end = clock();
    float seconds = (float)(end - start) / CLOCKS_PER_SEC;

    printf("Generated %d pairs.\n", n_pairs);
    printf("Found %d collisions.\n", n_colls);
    printf("Operation took %f ms\n\n\n", seconds*1000.0f);

    // clean up
    free(pairs);
}

void slow_test(size_t n_pts, float* xs, float* vs, float* rs) {
    printf("==================== SLOW TEST ====================\n");
    clock_t start = clock();

    size_t n_pairs;
    pair_t *pairs;
    make_collision_pairs_naiive(n_pts, xs, vs, rs, &n_pairs, &pairs);

    // count collisions
    size_t n_colls = count_collisions(n_pts, xs, vs, rs, n_pairs, pairs);

    clock_t end = clock();
    float seconds = (float)(end - start) / CLOCKS_PER_SEC;

    printf("Generated %d pairs.\n", n_pairs);
    printf("Found %d collisions.\n", n_colls);
    printf("Operation took %f ms\n\n\n", seconds*1000.0f);

    // clean up
    free(pairs);
}

int main() {
//    size_t arr_shape[4] = {6, 7, 8, 9};
//    size_t idxs[4] = {0};
//    unravel_index(1621, &arr_shape[0], &idxs[0]);
//    printf("idxs: %lu %lu %lu %lu\n", idxs[0], idxs[1], idxs[2], idxs[3]);
//    return 0;

    // read inputs
    float *xs, *vs, *rs;
    size_t n_rows, n_cols;
    if (!read_file("testing_coll/xs.txt", &xs, &n_rows, &n_cols)) {
        printf("Failed to read positions\n");
        return -1;
    }
    if (!read_file("testing_coll/vs.txt", &vs, &n_rows, &n_cols)) {
        printf("Failed to read positions\n");
        return -1;
    }
    if (!read_file("testing_coll/rs.txt", &rs, &n_rows, &n_cols)) {
        printf("Failed to read positions\n");
        return -1;
    }

    // first test
    fast_test(n_rows, xs, vs, rs);

    // second test
    slow_test(n_rows, xs, vs, rs);

    free(rs);
    free(vs);
    free(xs);

    return 0;
}
