#ifndef __TEST_FILE_H
#define __TEST_FILE_H

int file_exists(const char *fname);

const char* file_ext(const char *fname);

void file_read_vector(const char* fname, int* len, void** data, int eleSize);

void file_write_vector(const char* fname, int len, void* data, int eleSize);

#endif
