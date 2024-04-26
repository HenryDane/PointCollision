#include "collisions.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pair_set.h"

float signf(float x) {
    return x > 0 ? 1 : -1;
}

void unravel_index(size_t idx, int* arr_shape, int* idxs) {
//    printf("arr_shape=%d %d %d\n", arr_shape[0], arr_shape[1], arr_shape[2]);
    for (int i = 0; i < 3; i++) {
//        idxs[2 - i] = idx % arr_shape[2 - i];
        idxs[i] = idx % arr_shape[i];
//        idx = idx / arr_shape[2 - i];
        idx = idx / arr_shape[i];
    }
}

typedef struct {
    size_t n_elem;
    size_t n_capacity;
    size_t cell_index;
    size_t m_idxs[3];
    size_t* indexes;
} g2c_element;

const size_t g2c_size_incr = 1024;

void g2c_element_init(g2c_element* elem, size_t idx,
    size_t x, size_t y, size_t z) {
    elem->indexes = malloc(g2c_size_incr * sizeof(size_t));
    elem->n_capacity = g2c_size_incr;
    elem->n_elem = 0;
    elem->cell_index = idx;

    // figure out current cell coords
    elem->m_idxs[0] = x;
    elem->m_idxs[1] = y;
    elem->m_idxs[2] = z;
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
}

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

typedef struct {
    size_t a, b, c;
} triple_t ;

void make_collision_pairs(size_t n_pts, float* xs, float* vs, float* rs,
    size_t* n_pairs, pair_t** pairs) {
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
    const float sqrt3 = 1.73205080757;
    float L = fmax(max_vel * 3 * sqrt3, max_rad * 3 * sqrt3);
    printf("cell_size=%f\n", L);

    // find grid size
    size_t nc_x = 0;
    size_t nc_y = 0;
    size_t nc_z = 0;


}

void make_collision_pairs_old(size_t n_pts, float* xs, float* vs, float* rs,
    size_t* n_pairs, pair_t** pairs) {
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
    const float sqrt3 = 1.73205080757;
    float L = fmax(max_vel * 3 * sqrt3, max_rad * 3 * sqrt3);
//    float L = (max_vel + max_rad) * 2 * sqrt3;
//    float L = max_vel * 2;
    printf("cell_size=%f\n", L);

    // create cell coordinate array and find max values
    int* Xs = malloc(n_pts * 3 * sizeof(int));
    int n_cells[4] = {0};
    for (size_t i = 0; i < n_pts; i++) {
        for (size_t j = 0; j < 3; j++) {
            // compute integer coordinate
            Xs[(i * 3) + j] = (xs[(i * 3) + j] - min_xs[j]) / L;

            // find maximums;
            if (Xs[(i * 3) + j] > n_cells[j + 1]) {
                n_cells[j + 1] = Xs[(i * 3) + j];
            }
        }
    }

    // convert n_cells from max_idx to sizes
    size_t n_cells_total = 1;
    for (int i = 0; i < 4; i++) {
//        printf("i=%d n_cells[i]=%d\n", i, n_cells[i]);
        n_cells[i]++;
//        n_cells[i] += 2;
        n_cells_total *= n_cells[i];
    }
    printf("n_cells: %d %d %d\n", n_cells[1], n_cells[2], n_cells[3]);

    // cell_index --> pointer to array of indexes to points
    // size: (total cells, )
    g2c_element* g2c_table = malloc(n_cells_total * sizeof(g2c_element));
    for (size_t x = 0; x < n_cells[1]; x++) {
        for (size_t y = 0; y < n_cells[2]; y++) {
            for (size_t z = 0; z < n_cells[3]; z++) {
                size_t idx = x + (y * n_cells[1]) +
                    (z * n_cells[1] * n_cells[2]);
                g2c_element_init(&g2c_table[idx], idx, x, y, z);
            }
        }
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
    pair_set_t p_set;
    init_set(&p_set);
    for (size_t i = 0; i < n_cells_total; i++) {
//        if (g2c_table[i].n_elem == 0) continue;
        // count how many points we need to process
        size_t n_pts_tbl = 0;

        int X = g2c_table[i].m_idxs[0];
        int Y = g2c_table[i].m_idxs[1];
        int Z = g2c_table[i].m_idxs[2];

        // check neighbors to find count of points
        const int search_size = 1;
        for (int x = -search_size; x <= search_size; x++) {
            for (int y = -search_size; y <= search_size; y++) {
                for (int z = -search_size; z <= search_size; z++) {
                    int gx = X + x;
                    int gy = Y + y;
                    int gz = Z + z;
                    if (gx < 0 || gx >= n_cells[1]) continue;
                    if (gy < 0 || gy >= n_cells[2]) continue;
                    if (gz < 0 || gz >= n_cells[3]) continue;
                    size_t idx = gx + (gy * n_cells[1]) +
                        (gz * n_cells[1] * n_cells[2]);
//                    if (idx >= n_cells_total) continue;
                    n_pts_tbl += g2c_table[idx].n_elem;
                }
            }
        }

        size_t* points = malloc(n_pts_tbl * sizeof(size_t));
        size_t pidx = 0;
        for (int x = -search_size; x <= search_size; x++) {
            for (int y = -search_size; y <= search_size; y++) {
                for (int z = -search_size; z <= search_size; z++) {
                    int gx = X + x;
                    int gy = Y + y;
                    int gz = Z + z;
                    if (gx < 0 || gx >= n_cells[1]) continue;
                    if (gy < 0 || gy >= n_cells[2]) continue;
                    if (gz < 0 || gz >= n_cells[3]) continue;
                    size_t idx = gx + (gy * n_cells[1]) +
                        (gz * n_cells[1] * n_cells[2]);
//                    if (idx >= n_cells_total) continue;
                    // copy
                    for (size_t j = 0; j < g2c_table[idx].n_elem; j++) {
                        if (g2c_table[idx].indexes[j] >= n_pts) {
                            printf("g2c idx[j]=%d n_pts=%d\n",
                                g2c_table[idx].indexes[j], n_pts);
                            break; //exit(3);
                        }
                        points[pidx++] = g2c_table[idx].indexes[j];
                    }
                }
            }
        }

        printf("pidx=%d n_pts_tbl=%d\n", pidx, n_pts_tbl);
        printf("inserting %d pairs for i=%d\n", pidx * pidx / 2, i);
        for (size_t j = 0; j < pidx; j++) {
            for (size_t k = j + 1; k < pidx; k++) {
                if (points[j] >= n_pts) {
                    printf("points[j]=%d, n_pts=%d\n", points[j], n_pts);
                    exit(3);
                }
                if (points[k] >= n_pts) {
                    printf("points[k]=%d, n_pts=%d\n", points[k], n_pts);
                    exit(3);
                }
                pair_t p = {points[j], points[k]};
                insert_pair(&p_set, &p);
            }
        }
        printf("    set size=%d\n", p_set.n);

        free(g2c_table[i].indexes);
    }

    to_flat_array(&p_set, n_pairs, pairs);

    free(g2c_table);
    free(Xs);
    free(cell_indexes);
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
