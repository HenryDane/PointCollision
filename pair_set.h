#ifndef _PAIR_SET_H
#define _PAIR_SET_H

typedef struct {
    size_t a, b;
} pair_t;

size_t hash_pair(pair_t* p);
bool hash_compare(pair_t a, pair_t b);

#define NUM_HASH_BUCKETS 1024

typedef struct {
    pair_set_node_t* next;
    pair_t data;
} pair_set_node_t;

typedef struct {
    pair_set_node_t* nodes[NUM_HASH_BUCKETS];
} pair_set_t;

void init_set(pair_set_t* set);
void clear_set(pair_set_t* set);
void free_set(pair_set_t* set);

void insert_pair(pair_set_t* set, pair_t* pair);

#endif
