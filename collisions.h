#ifndef _COLLISIONS_H
#define _COLLISIONS_H

#include <stddef.h>
#include <stdbool.h>

bool is_colliding(size_t A, size_t B, float* xs, float* vs, float* rs,
    float* t);

void make_collision_pairs(size_t n_pts, float* xs, float* vs, float* rs,
    size_t* n_pairs, size_t** pairs);

size_t count_collisions(size_t n_pts, float* xs, float* vs, float* rs,
    size_t n_pairs, size_t* pairs);

#endif //_COLLISIONS_H
