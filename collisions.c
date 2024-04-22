#include "collisions.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

float signf(float x) {
    return x > 0 ? 1 : -1;
}

void unravel_index(size_t idx, size_t* arr_shape, size_t* idxs) {
    for (int i = 0; i < 3; i++) {
        idxs[3 - i] = idx % arr_shape[3 - i];
        idx = idx / arr_shape[3 - i];
    }
}

typedef struct {
    size_t n_elem;
    size_t n_capacity;
    size_t cell_index;
    size_t n_adjacent;
    size_t* adjacents;
    size_t* indexes;
} g2c_element;

const size_t g2c_size_incr = 1024;

void g2c_element_init(g2c_element* elem, size_t idx, size_t* n_cells) {
    elem->indexes = malloc(g2c_size_incr * sizeof(size_t));
    elem->n_capacity = g2c_size_incr;
    elem->n_elem = 0;
    elem->cell_index = idx;

    // figure out current cell coords
//    size_t m_idxs[3];
//    unravel_index(idx, n_cells, &m_idxs);
}

void g2c_element_emplace(g2c_element* elem, size_t index) {
    // make sure there is size
    if (elem->n_elem == elem->n_capacity) {
        elem->n_capacity += g2c_size_incr;
        size_t* upd_arr = realloc(elem->indexes,
            elem->n_capacity * sizeof(size_t));
        if (upd_arr) {
            elem->indexes = upd_arr;
        } else {
            printf("Realloc failed.\n");
            exit(2);
        }
    }
    // assign stuff
    elem->indexes[elem->n_elem] = index;
    elem->n_elem++;
};

bool pairs_equal(pair_t* A, pair_t* B) {
    return ((A->a == B->a) && (A->b == B->b)) ||
        ((A->b == B->a) && (A->a == B->b));
}

bool is_colliding(size_t A, size_t B, float* xs, float* vs, float* rs,
    float* t) {
    // get r sum
    float r_sum = rs[A] + rs[B];
    // convert from index to offset
    A *= 3; B *= 3;
    // compute \vec{1} = c_{A} - c_{B}
    float one[3] = {
        xs[A + 0] - xs[B + 0],
        xs[A + 1] - xs[B + 1],
        xs[A + 2] - xs[B + 2]
    };
    // compute v_{AB} = v_{A} - v_{B}
    float vAB[3] = {
        vs[A + 0] - vs[B + 0],
        vs[A + 1] - vs[B + 1],
        vs[A + 2] - vs[B + 2]
    };
    // compute a, b, c
    float a = (vAB[0] * vAB[0]) + (vAB[1] * vAB[1]) + (vAB[2] * vAB[2]);
    float b = 2.0f *
        ((one[0] * vAB[0]) + (one[1] * vAB[1]) + (one[2] * vAB[2]));
    float c = ((one[0] * one[0]) + (one[1] * one[1]) + (one[2] * one[2])) -
        (r_sum * r_sum);
    // static check
    if (c < 0) {
        *t = -1.0;
        return true;
    }
    // compute q
    float q = -0.5f * (b + signf(b) * sqrtf((b * b) - (4.0 * a * c)));
    // compute t0, t1
    float t0 = q / a;
    float t1 = c / q;
    // compute t
    *t = fmin(t0, t1);
    *t = ((*t) < 0) ? NAN : *t;
    *t = ((*t) > 1) ? NAN : *t;
    // check if collided
    return isfinite(*t);
}

void make_collision_pairs_naiive(size_t n_pts, float* xs, float* vs, float* rs,
    size_t* n_pairs, pair_t** pairs) {
    (*n_pairs) = ((n_pts * n_pts) / 2);
    (*pairs) = malloc((*n_pairs) * sizeof(pair_t));
    size_t idx = 0;
    for (size_t i = 0; i < n_pts; i++) {
        for (size_t j = i + 1; j < n_pts; j++) {
            (*pairs)[idx].a = i;
            (*pairs)[idx].a = j;
            idx++;
        }
    }
}

void make_collision_pairs(size_t n_pts, float* xs, float* vs, float* rs,
    size_t* n_pairs, pair_t** pairs) {
    // generate offset array
    int* offs = malloc(27 * 3 * sizeof(int));
    for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
            for (int k = -1; k <= 1; k++) {
                int idx = (i + 1) + ((j + 1) * 2) + ((k + 1) * 9);
                offs[(idx * 3) + 0] = i;
                offs[(idx * 3) + 1] = j;
                offs[(idx * 3) + 2] = k;
            }
        }
    }

    // find max velocity and max radius
    float max_vel, max_rad; max_vel = -1e99; max_rad = -1e99;
    float min_xs[3] = {1e99};
    for (size_t i = 0; i < n_pts; i++) {
        // compute offset index
        size_t idx = i * 3;
        // compute velocity
        float vel = sqrtf( (vs[idx] * vs[idx]) + (vs[idx + 1] * vs[idx + 1]) +
            (vs[idx + 2] * vs[idx + 2]) );
        // compute max velocity
        max_vel = vel     > max_vel ? vel     : max_vel;
        // find max radius
        max_rad = rs[i]   > max_rad ? rs[i]   : max_rad;
        // compute min position(s)
        for (int j = 0; j < 3; j++) {
            if (xs[(i * 3) + j] < min_xs[j]) {
                min_xs[j] = xs[(i * 3) + j];
            }
        }
    }
    printf("max_vel=%f max_rad=%f\n", max_vel, max_rad);
    printf("min xs: %f %f %f\n", min_xs[0], min_xs[1], min_xs[2]);

    // compute grid cell size
//    float cube_root_of_2 = powf(2, -3);
//    printf("cube_root_of_2=%f\n", cube_root_of_2);
    float L = fmax(max_vel * 2 , max_rad * 2);
    printf("cell_size=%f\n", L);

    // create cell coordinate array and find max values
    int* Xs = malloc(n_pts * 3 * sizeof(int));
    int n_cells[4] = {0};
    for (size_t i = 0; i < n_pts; i++) {
        for (size_t j = 0; j < 3; j++) {
            // compute integer coordinate
            Xs[(i * 3) + j] = round(xs[(i * 3) + j] - min_xs[j]) / L;

            // find maximums;
            if (Xs[(i * 3) + j] > n_cells[j + 1]) {
                n_cells[j + 1] = Xs[(i * 3) + j];
            }
        }
    }

    // convert n_cells from max_idx to sizes
    size_t n_cells_total = 1;
    for (int i = 0; i < 4; i++) {
        n_cells[i]++;
        n_cells_total *= n_cells[i];
    }
    printf("n_cells: %d %d %d\n", n_cells[1], n_cells[2], n_cells[3]);

    // cell_index --> pointer to array of indexes to points
    // size: (total cells, )
    g2c_element* g2c_table = malloc(n_cells_total * sizeof(g2c_element));
    for (size_t i = 0; i < n_cells_total; i++) {
        g2c_element_init(&g2c_table[i], i, &n_cells[0]);
    }

    // build g2c table and compute point cell indexes
    size_t* cell_indexes = malloc(n_pts * sizeof(size_t));
    for (size_t i = 0; i < n_pts; i++) {
        cell_indexes[i] = 0;
        for (int j = 0; j < 3; j++) {
            cell_indexes[i] += Xs[(i * 3) + j] * n_cells[j];
        }
        g2c_element_emplace(&g2c_table[cell_indexes[i]], i);
    }

    // start generating pairs
    size_t pair_capacity = 1024;
    *n_pairs = 0;
    (*pairs) = malloc(1024 * sizeof(pair_t));
    for (size_t i = 0; i < n_cells_total; i++) {
        if (g2c_table[i].n_elem == 0) continue;

        for (size_t j = 0; j < g2c_table[i].n_elem; j++) {
            for (size_t k = j + 1; k < g2c_table[i].n_elem; k++) {
                // check size
                if ((*n_pairs) == pair_capacity) {
                    pair_capacity += 1024;
                    size_t* upd_arr = realloc(*pairs,
                        pair_capacity * 2 * sizeof(size_t));
                    if (upd_arr) {
                        *pairs = upd_arr;
                    } else {
                        printf("Pair realloc failed.\n");
                        exit(2);
                    }
                }

                // insert the pair
//                (*pairs)[(*n_pairs) * 2 + 0] = g2c_table[i].indexes[j];
//                (*pairs)[(*n_pairs) * 2 + 1] = g2c_table[i].indexes[k];
                (*pairs)[(*n_pairs)].a = g2c_table[i].indexes[j];
                (*pairs)[(*n_pairs)].b = g2c_table[i].indexes[k];

                // increment pair counter
                (*n_pairs)++;
            }
        }

        free(g2c_table[i].indexes);
    }

    // check for duplicates
//    size_t n_dupl = 0;
//    for (size_t i = 0; i < *n_pairs; i++) {
//        // look for occurrences of this element elsewhere
//        for (size_t j = i + 1; j < *n_pairs; j++) {
//            // if pairs[i] == pairs[j], swap the jth with the last and
//            // repeat until they do not match, decrementing n_pairs each time
//            while (pairs_equal(&(*pairs)[i], &(*pairs)[j])) {
////                printf("match\n");
//                // swap first element
//                pair_t tmp = {(*pairs)[i].a, (*pairs)[i].b}; // = pairs[i];
//                (*pairs)[i].a = (*pairs)[j].a;
//                (*pairs)[i].b = (*pairs)[j].b;
//                (*pairs)[j].a = tmp.a;
//                (*pairs)[j].b = tmp.b;
//
//                // decrement size
//                (*n_pairs)--;
//
//                // count dupl
//                n_dupl++;
//            }
//        }
//    }
//    printf("found %d dupl\n", n_dupl);

    free(g2c_table);
    free(Xs);
    free(cell_indexes);
    free(offs);
}

size_t count_collisions(size_t n_pts, float* xs, float* vs, float* rs,
    size_t n_pairs, pair_t* pairs) {
    size_t dy_colls = 0;
    size_t st_colls = 0;

    float t;
    for (size_t i = 0; i < n_pairs; i++) {
        if (is_colliding(pairs[i].a, pairs[i].b, xs, vs, rs, &t)) {
            if (t < 0) {
                st_colls++;
            } else {
                dy_colls++;
            }
        } else {
            continue;
        }
//        printf("%lu %lu %f\n", pairs[i].a, pairs[i].b, t);
    }
    //printf("\n");

    printf("static: %d\ndynamic: %d\ntotal: %d\n", st_colls, dy_colls,
        st_colls + dy_colls);

    return dy_colls + st_colls;
}
