#include "pair_set.h"
#include <stdlib.h>

size_t hash_pair(pair_t* p) {
    return (p->a) ^ (p->b);
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

void init_set(pair_set_t* set) {
    for (size_t i = 0; i < NUM_HASH_BUCKETS; i++) {
        set->nodes[i] = NULL;
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
    set->n++;
    size_t i = hash_pair(pair) % NUM_HASH_BUCKETS;
    if (set->nodes[i] == NULL) {
        printf("first insert\n");
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
        if (hash_compare(current->data, *pair)) {
            printf("skipping duplicate\n");
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
