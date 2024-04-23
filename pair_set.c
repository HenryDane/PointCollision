#include "pair_set.h"
#include <stdlib.h>

size_t hash_pair(pair_t* p) {
//    return (p->a) ^ (p->b);
    return (p->a << 32) | p->b;
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

#ifdef USE_ARRAY_HASH_SET
void init_set(pair_set_t* set) {
    for (size_t i = 0; i < NUM_HASH_BUCKETS; i++) {
        set->nodes[i] = (pair_set_node_t*) malloc(sizeof(pair_set_node_t));
        set->nodes[i]->capacity = INIT_HASH_ARR_SIZE;
        set->nodes[i]->index = 0;
        set->nodes[i]->pairs =
            (pair_t*) malloc(set->nodes[i]->capacity * sizeof(pair_t));
    }
    set->n = 0;
}
#else
void init_set(pair_set_t* set) {
    for (size_t i = 0; i < NUM_HASH_BUCKETS; i++) {
        set->nodes[i] = NULL;
    }
    set->n = 0;
}
#endif
void clear_set(pair_set_t* set) {
    // TODO
}

void free_set(pair_set_t* set) {
    // TODO
}

#ifdef USE_ARRAY_HASH_SET
void insert_pair(pair_set_t* set, pair_t* pair) {
    // get current bucket
    pair_set_node_t* b = set->nodes[hash_pair(pair) % NUM_HASH_BUCKETS];

    // check if we need to realloc
    if (b->index >= b->capacity) {
//        printf("reallocating\n");
        // grow exponentially
        b->capacity += INIT_HASH_ARR_SIZE; // 274 ms avg
//        b->capacity *= 2; // 283ms avg
        // perform reallocation
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
//    printf("checking duplicates\n");
    for (size_t i = 0; i < b->index; i++) {
        if (pair_equal(b->pairs[i], *pair)) {
//            printf("found duplicate\n");
            return;
        }
    }

    // insert the pair
//    printf("inserting pair\n");
    b->pairs[b->index].a = (*pair).a;
    b->pairs[b->index].b = (*pair).b;

    // increment size
    b->index++;
    set->n++;
}

void to_flat_array(pair_set_t* set, size_t* n_pairs, pair_t** pairs) {
    *n_pairs = set->n;
    (*pairs) = (pair_t*) malloc((*n_pairs) * sizeof(pair_t));
    size_t i = 0;
    pair_set_node_t* node = NULL;
    for (size_t b = 0; b < NUM_HASH_BUCKETS; b++) {
        node = set->nodes[b];
        for (size_t j = 0; j < node->index; j++) {
            (*pairs)[i] = node->pairs[j];
            i++;
        }
    }
}
#else
void insert_pair(pair_set_t* set, pair_t* pair) {
    set->n++;
    size_t i = hash_pair(pair) % NUM_HASH_BUCKETS;
    if (set->nodes[i] == NULL) {
//        printf("first insert\n");
        // insert first into this bucket
        set->nodes[i] = (pair_set_node_t*) malloc(sizeof(pair_set_node_t));
        set->nodes[i]->data.a = pair->a;
        set->nodes[i]->data.b = pair->b;
        set->nodes[i]->next = NULL;
        return;
    }

    // find and insert
    pair_set_node_t* current = set->nodes[i];
    while (current->next != NULL) {
        if (pair_equal(current->data, *pair)) {
//            printf("skipping duplicate\n");
            return;
        }
//        printf("looking...\n");
        current = current->next;
    }
//    printf("done...\n");

    current->next = (pair_set_node_t*) malloc(sizeof(pair_set_node_t));
    current->next->data.a = pair->a;
    current->next->data.b = pair->b;
    current->next->next = NULL;
}

void to_flat_array(pair_set_t* set, size_t* n_pairs, pair_t** pairs) {
    *n_pairs = set->n;
    (*pairs) = (pair_t*) malloc((*n_pairs) * sizeof(pair_t));
    size_t i = 0;
    pair_set_node_t* node = NULL;
    for (size_t b = 0; b < NUM_HASH_BUCKETS; b++) {
        node = set->nodes[b];
        while (node != NULL) {
            (*pairs)[i].a = node->data.a;
            (*pairs)[i].b = node->data.b;
            i++;
            node = node->next;
        }
    }
}
#endif
