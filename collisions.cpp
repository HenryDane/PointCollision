#include "collisions.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <set>
#include <vector>
#include <algorithm>

float signf(float x) {
    return x > 0 ? 1 : -1;
}

typedef struct {
    size_t n_elem;
    size_t n_capacity;
    size_t idx_x, idx_y, idx_z;
    size_t* indexes;
} g2c_element;

const size_t g2c_size_incr = 1024;

void g2c_element_init(g2c_element* elem, size_t x, size_t y, size_t z) {
    elem->indexes = (size_t*) malloc(g2c_size_incr * sizeof(size_t));
    elem->n_capacity = g2c_size_incr;
    elem->n_elem = 0;

    // figure out current cell coords
    elem->idx_x = x;
    elem->idx_y = y;
    elem->idx_z = z;
}

void g2c_element_emplace(g2c_element* elem, size_t index) {
    // make sure there is size
    if (elem->n_elem == elem->n_capacity) {
        elem->n_capacity += g2c_size_incr;
        size_t* upd_arr = (size_t*) realloc(elem->indexes,
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
    //(*n_pairs) = ((n_pts * n_pts) / 2);
    (*n_pairs) = (n_pts * (n_pts - 1)) / 2;
    (*pairs) = (pair_t*) malloc((*n_pairs) * sizeof(pair_t));
    size_t idx = 0;
    for (size_t i = 0; i < n_pts; i++) {
        for (size_t j = i + 1; j < n_pts; j++) {
            (*pairs)[idx].a = i;
            (*pairs)[idx].b = j;
            idx++;
        }
    }
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

    printf("static: %lu\ndynamic: %lu\ntotal: %lu\n", st_colls, dy_colls,
        st_colls + dy_colls);

    return dy_colls + st_colls;
}

size_t count_collisions_fast(size_t n_pts, float* xs, float* vs, float* rs) {
        // find max velocity and max radius
    float max_vel, max_rad; max_vel = -1e99; max_rad = -1e99;
    float min_x = 1e99;
    float min_y = 1e99;
    float min_z = 1e99;
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
        min_x = xs[idx + 0] < min_x ? xs[idx + 0] : min_x;
        min_y = xs[idx + 1] < min_y ? xs[idx + 1] : min_y;
        min_z = xs[idx + 2] < min_z ? xs[idx + 2] : min_z;
    }

    // compute grid cell size
    const float sqrt3 = 1.73205080757;
    float L = fmax(max_vel * 1.0 * sqrt3, max_rad * 1.0 * sqrt3);

    // prepare grid size
    size_t nc_x = 0;
    size_t nc_y = 0;
    size_t nc_z = 0;

    // prepare integer coordinate array
    int* Xs = (int*) malloc(n_pts * 3 * sizeof(int));

    // find grid size and generate integer coord array
    for (size_t i = 0; i < n_pts; i++) {
        size_t j = i * 3;

        // compute integer coordinates
        Xs[j + 0] = (xs[j + 0] - min_x) / L;
        Xs[j + 1] = (xs[j + 1] - min_y) / L;
        Xs[j + 2] = (xs[j + 2] - min_z) / L;

        // find grid cells
        nc_x = Xs[j + 0] > nc_x ? Xs[j + 0] : nc_x;
        nc_y = Xs[j + 1] > nc_y ? Xs[j + 1] : nc_y;
        nc_z = Xs[j + 2] > nc_z ? Xs[j + 2] : nc_z;
    }

    // sizes are currently the max valid index along axis, so incr by one to
    // convert them to sizes
    nc_x++; nc_y++; nc_z++;

    // calculate total number of cells
    size_t n_cells_total = nc_x * nc_y * nc_z;

    // create grid table -- this maps grid cell index to a list of points
    // that are in that grid cell
    g2c_element* g2c_table = (g2c_element*) malloc(n_cells_total * sizeof(g2c_element));

    // initialize grid table
    for (size_t x = 0; x < nc_x; x++) {
        for (size_t y = 0; y < nc_y; y++) {
            for (size_t z = 0; z < nc_z; z++) {
                size_t idx = x + (y * nc_x) +
                    (z * nc_x * nc_y);
                g2c_element_init(&g2c_table[idx], x, y, z);
            }
        }
    }

    // build g2c table and compute point cell indexes
    for (size_t i = 0; i < n_pts; i++) {
        size_t idx = Xs[(i * 3) + 0] + (Xs[(i * 3) + 1] * nc_x) +
            (Xs[(i * 3) + 2] * nc_x * nc_y);
        g2c_element_emplace(&g2c_table[idx], i);
    }

    // start generating pairs
    std::set<pair_t> p_set;
    printf("max set size=%lu\n", p_set.max_size());
//    std::vector<pair_t> p_set;
//    p_set.reserve(n_pts / 10); // guess some initial size
    for (size_t i = 0; i < n_cells_total; i++) {
        // dont process empty grid cells
        if (g2c_table[i].n_elem == 0) continue;

        // count how many points we need to process
        size_t n_pts_tbl = 0;

        int X = g2c_table[i].idx_x;
        int Y = g2c_table[i].idx_y;
        int Z = g2c_table[i].idx_z;

        // check neighbors to find count of points
        const int search_size = 1;
        for (int x = -search_size; x <= search_size; x++) {
            for (int y = -search_size; y <= search_size; y++) {
                for (int z = -search_size; z <= search_size; z++) {
                    int gx = X + x;
                    int gy = Y + y;
                    int gz = Z + z;
                    if (gx < 0 || gx >= nc_x) continue;
                    if (gy < 0 || gy >= nc_y) continue;
                    if (gz < 0 || gz >= nc_z) continue;
                    size_t idx = gx + (gy * nc_x) + (gz * nc_x * nc_y);
                    n_pts_tbl += g2c_table[idx].n_elem;
                }
            }
        }

        size_t* points = (size_t*) malloc(n_pts_tbl * sizeof(size_t));
        size_t pidx = 0;
        for (int x = -search_size; x <= search_size; x++) {
            for (int y = -search_size; y <= search_size; y++) {
                for (int z = -search_size; z <= search_size; z++) {
                    int gx = X + x;
                    int gy = Y + y;
                    int gz = Z + z;
                    if (gx < 0 || gx >= nc_x) continue;
                    if (gy < 0 || gy >= nc_y) continue;
                    if (gz < 0 || gz >= nc_z) continue;
                    size_t idx = gx + (gy * nc_x) + (gz * nc_x * nc_y);

                    // copy
                    for (size_t j = 0; j < g2c_table[idx].n_elem; j++) {
                        points[pidx++] = g2c_table[idx].indexes[j];
                    }
                }
            }
        }

        for (size_t j = 0; j < pidx; j++) {
            for (size_t k = j + 1; k < pidx; k++) {
                pair_t p = {points[j], points[k]};
                float t;
                if (!is_colliding(points[j], points[k], xs, vs, rs, &t)) {
                    continue;
                }
//                insert_pair(&p_set, &p);
                p_set.insert(p);
//                p_set.emplace_back(std::move(p));
            }
        }

//        p_set.erase(std::unique(p_set.begin(), p_set.end()), p_set.end());
//        std::sort(p_set.begin(), p_set.end());
//        p_set.erase(std::unique(p_set.begin(), p_set.end()), p_set.end());

        free(points);
    }

    return p_set.size();
}
