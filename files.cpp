#include "files.hpp"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

bool read_file(const char* path, float** data,
    size_t* n_rows, size_t* n_cols) {
    const int max_len = 256;

    // open file
    FILE* fp = fopen(path, "r");
    if (!fp) {
        printf("Error: Could not open file: %s\n", path);
        return false;
    }

    // read header
    if (fscanf(fp, "# %lu %lu\n", n_rows, n_cols) != 2) {
        printf("Failed to parse header for file: %s", path);
        return false;
    }

    // allocate memory
    printf("n_rows=%lu n_cols=%lu\n", *n_rows, *n_cols);
    if (*n_cols <= 1) *n_cols = 1;
    (*data) = (float*) malloc( (*n_rows) * (*n_cols) * sizeof(float) );
    size_t idx = 0;

    char buffer[max_len];
    while (fgets(buffer, max_len, fp)) {
        //printf("new row\n");
        buffer[strcspn(buffer, "\n")] = 0;

        // try to find n_cols of floats
        char* token = strtok(buffer, " ");
        float f;
        if (sscanf(token, "%f", &f) != 1) {
            printf("Failed parsing float\n");
            return false;
        }
        (*data)[idx] = f;
        idx++;
        for (int i = 1; i < *n_cols; i++) {
            token = strtok(NULL, " ");
            if (!token) {
                printf("Token is null\n");
                return false;
            }
            if (sscanf(token, "%f", &f) != 1) {
                printf("Failed parsing float\n");
                return false;
            }
            (*data)[idx] = f;
            idx++;
        }
    }

    printf("idx=%lu\n", idx);

    fclose(fp);

    return true;
}
