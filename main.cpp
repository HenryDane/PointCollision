#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "files.hpp"
#include "collisions.hpp"

void fast_test(size_t n_pts, float* xs, float* vs, float* rs) {
    printf("==================== FAST TEST ====================\n");
    clock_t start = clock();

    /*
    size_t n_pairs;
    pair_t *pairs;
    make_collision_pairs(n_pts, xs, vs, rs, &n_pairs, &pairs);

    // count collisions
    size_t n_colls = count_collisions(n_pts, xs, vs, rs, n_pairs, pairs);
    */

    size_t n_pairs = 0;
    size_t n_colls = count_collisions_fast(n_pts, xs, vs, rs);

    clock_t end = clock();
    float seconds = (float)(end - start) / CLOCKS_PER_SEC;

    printf("Generated %lu pairs.\n", n_pairs);
    printf("Found %lu collisions.\n", n_colls);
    printf("Operation took %f ms\n\n\n", seconds*1000.0f);

    // clean up
//    free(pairs);
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

    printf("Generated %lu pairs.\n", n_pairs);
    printf("Found %lu collisions.\n", n_colls);
    printf("Operation took %f ms\n\n\n", seconds*1000.0f);

    // clean up
    free(pairs);
}

float rand_float(float min, float max) {
    return min + ((rand() / ((float) RAND_MAX)) * (max - min));
}

//float time_test() {
//    size_t n_pts = 50000;
//    float* xs = (float*) malloc(n_pts * 3 * sizeof(float));
//    float* vs = (float*) malloc(n_pts * 3 * sizeof(float));
//    float* rs = (float*) malloc(n_pts * sizeof(float));
//    for (size_t i = 0; i < n_pts; i++) {
//        size_t j = i * 3;
//        // positions [-10, 10]
//        xs[j + 0] = rand_float(-10, 10);
//        xs[j + 1] = rand_float(-10, 10);
//        xs[j + 2] = rand_float(-10, 10);
//        // velocities [-1, 1] (v_z=0)
//        vs[j + 0] = rand_float(-0.1, 0.01);
//        vs[j + 1] = rand_float(-0.1, 0.01);
//        vs[j + 2] = 0.0f;
//        // radii
//        rs[i] = 1.0f;
//    }
//
//    clock_t start = clock();
//
//    size_t n_pairs;
//    pair_t *pairs;
//    make_collision_pairs(n_pts, xs, vs, rs, &n_pairs, &pairs);
//
//    clock_t end = clock();
//    float seconds = (float)(end - start) / CLOCKS_PER_SEC;
//
//    free(rs);
//    free(vs);
//    free(xs);
//
//    return seconds;
//}

int main() {
//    pair_set_t pset;
//    init_set(&pset);
//
//    for (size_t i = 0; i < 100; i++) {
//        size_t a = rand();
//        size_t b = rand();
//        pair_t pa = {a, b};
//        pair_t pb = {b, a};
//        insert_pair(&pset, &pa);
//        insert_pair(&pset, &pb);
//    }
//    printf("pset size=%d\n", pset.n);
//
//    return 0;

    // read inputs
#define GENERATE_TEST_POINTS
#ifndef GENERATE_TEST_POINTS
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
#else
    size_t n_rows = 50000;
    float* xs = (float*) malloc(n_rows * 3 * sizeof(float));
    float* vs = (float*) malloc(n_rows * 3 * sizeof(float));
    float* rs = (float*) malloc(n_rows * sizeof(float));
    for (size_t i = 0; i < n_rows; i++) {
        size_t j = i * 3;
        // positions [-10, 10]
        xs[j + 0] = rand_float(-10, 10);
        xs[j + 1] = rand_float(-10, 10);
        xs[j + 2] = rand_float(-10, 10);
        // velocities [-1, 1] (v_z=0)
        vs[j + 0] = rand_float(-0.1, 0.01);
        vs[j + 1] = rand_float(-0.1, 0.01);
        vs[j + 2] = 0.0f;
        // radii
        rs[i] = 1.0f;
    }
#endif

    // first test
    fast_test(n_rows, xs, vs, rs);

    // second test
    slow_test(n_rows, xs, vs, rs);

    return 0;
//
//    // timings
//    printf("Beginning timing tests\n");
//    float time_acc = 0;
//    const int n_tests = 50;
//    for (size_t i = 0; i < n_tests; i++) {
//        printf("Test %lu\n", i);
//        time_acc += time_test();
//    }
//    time_acc /= ((float) n_tests);
//    printf("\n\nDid %d tests in %.4f ms (avg).\n", n_tests, time_acc * 1000.0);

    free(rs);
    free(vs);
    free(xs);

    return 0;
}
