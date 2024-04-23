#include "pair_set.h"

size_t hash_pair(pair_t* p) {
    return (p->a << 8) ^ (p->b);
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
}

void clear_set(pair_set_t* set) {
    // TODO
}

void free_set(pair_set_t* set) {
    // TODO
}

void insert_pair(pair_set_t* set, pair_t* pair) {
    size_t i = hash_pair(pair) % NUM_HASH_BUCKETS;
    if (set->nodes[i] == NULL) {
        // insert first into this bucket
        set->nodes[i] = malloc(sizeof(pair_set_node_t));
        set->nodes[i]->data.a = pair->a;
        set->nodes[i]->data.b = pair->b;
        set->nodes[i]->next = NULL;
        return;
    }

    // find and insert
    pair_set_node_t* current = set->nodes[i];
    while (current->next != NULL) {
        current = current->next;
    }

    current->next = malloc(sizeof(pair_set_node_t));
    current->next->data.a = pair->a;
    current->next->data.b = pair->b;
    current->next->next = NULL;
}
