#include "pair_set.h"
#include <stdlib.h>
#include <stdio.h>

size_t hash_pair(pair_t* p) {
    size_t a, b;
    if (p->a < p->b) {
        a = p->a;
        b = p->b;
    } else {
        a = p->b;
        b = p->a;
    }
    return (a << 32) | b;
}

bool hash_compare(pair_t a, pair_t b) {
    if (a.a < b.a) {
        return true;
    } else if (a.a > b.a) {
        return false;
    } else {
        return (a.b < b.b);
    }
}

bool pair_equal(pair_t a, pair_t b) {
    return ((a.a == b.a) && (a.b == b.b)) ||
        ((a.a == b.b) && (a.b == b.a));
}

//bool pair_equal(pair_t* a, pair_t* b) {
//    return ((a->a == b->a) && (a->b == b->b)) ||
//        ((a->a == b->b) && (a->b == b->a));
//}

void init_set(pair_set_t* set) {
    set->nodes = malloc(NUM_HASH_BUCKETS * sizeof(pair_set_node_t*));
    for (size_t i = 0; i < NUM_HASH_BUCKETS; i++) {
        set->nodes[i] = (pair_set_node_t*) malloc(sizeof(pair_set_node_t));
        set->nodes[i]->capacity = INIT_HASH_ARR_SIZE;
        set->nodes[i]->index = 0;
        set->nodes[i]->pairs =
            (pair_t*) malloc(set->nodes[i]->capacity * sizeof(pair_t));
    }
    set->n = 0;
}

void clear_set(pair_set_t* set) {
    // TODO
}

void free_set(pair_set_t* set) {
    // TODO
}

void insert_pair(pair_set_t* set, pair_t* pair) {
    // get current bucket
    pair_set_node_t* b = set->nodes[hash_pair(pair) % NUM_HASH_BUCKETS];

    // check if we need to realloc
    if (b->index >= b->capacity) {
        // grow exponentially
        b->capacity *= 2;
        pair_t* newpairs = realloc(b->pairs, b->capacity * sizeof(pair_t));
        // make sure new pointer is useable
        if (newpairs == NULL) {
            printf("set realloc failed.\n");
            exit(2);
        }
        // update this bucket
        b->pairs = newpairs;
    }

    // check if pair exists
    for (size_t i = 0; i < b->index; i++) {
        if (pair_equal(b->pairs[i], *pair)) {
            return;
        }
    }

    // insert the pair
    b->pairs[b->index].a = (*pair).a;
    b->pairs[b->index].b = (*pair).b;

    // increment size
    b->index++;
    set->n++;
}

void to_flat_array(pair_set_t* set, size_t* n_pairs, pair_t** pairs) {
    (*n_pairs) = set->n;
    (*pairs) = (pair_t*) malloc(set->n * sizeof(pair_t));
    size_t i = 0;
    for (size_t b = 0; b < NUM_HASH_BUCKETS; b++) {
        pair_set_node_t* node = set->nodes[b];
        for (size_t j = 0; j < node->index; j++) {
            (*pairs)[i] = node->pairs[j];
            i++;
        }
    }
//    printf("i=%lu n=%lu\n", i, set->n);
}
