#ifndef _FILES_H
#define _FILES_H

#include <stdbool.h>
#include <stddef.h>

bool read_file(const char* path, float** data, size_t* n_rows, size_t* n_cols);

#endif // _FILES_H
