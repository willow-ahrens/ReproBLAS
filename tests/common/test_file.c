#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "test_file.h"
int file_exists(const char *fname)
{
  FILE *file = fopen(fname, "r");
  if (file)
  {
    fclose(file);
    return 1;
  }
  return 0;
}

const char* file_ext(const char* fname) {
  int i;
  for(i = strlen(fname) - 1; i >= 0; i--){
    if(fname[i] == '.'){
      return fname + i;
    }
  }
  return fname + strlen(fname) - 1;
}

void file_read_vector(const char* fname, int* len, void** data, int eleSize) {
  int i;
  FILE* file;

  //check extension
  if(strcmp(file_ext(fname), ".dat") != 0){
    fprintf(stderr, "error: invalid file type %s\n", file_ext(fname));
    exit(125);
  }

  file = fopen(fname, "rb");
  if(file == NULL) {
    fprintf(stderr,"error: failed to open file %s to read.\n", fname);
    exit(125);
  }
  if(fread(len, sizeof(int), 1, file) != 1){
    fprintf(stderr,"error: failed to read file %s header\n", fname);
    exit(125);
  }
  *data = malloc(*len * eleSize);
  char* ptr = *data;
  for(i = 0; i < *len; i++, ptr += eleSize) {
    if(fread(ptr, eleSize, 1, file) != 1){
      fprintf(stderr,"error: failed to read file %s (incorrect element size %d?)\n", fname, eleSize);
      exit(125);
    }
  }
  if(fread(ptr, eleSize, 1, file) != 0){
    fprintf(stderr,"error: failed to read file %s (incorrect element size %d?)\n", fname, eleSize);
    exit(125);
  }
  fclose(file);
}

void file_write_vector(const char* fname, int len, void* data, int eleSize) {
  int i;
  FILE* file;

  //check extension
  if(strcmp(file_ext(fname), ".dat") != 0){
    fprintf(stderr, "error: invalid file type %s\n", file_ext(fname));
    exit(125);
  }

  file = fopen(fname, "wb");
  if(file == NULL) {
    fprintf(stderr,"error: failed to open file %s to write.\n", fname);
    exit(125);
  }
  fwrite(&len, sizeof(int), 1, file);
  char* ptr = data;
  for(i = 0; i < len; i++, ptr += eleSize) {
    fwrite(ptr, eleSize, 1, file);
  }
  fclose(file);

  //check written data
  void *retrieved_data;
  int retrieved_len;
  file_read_vector(fname, &retrieved_len, &retrieved_data, eleSize);
  if (retrieved_len != len) {
    fprintf(stderr, "error: failed to write data correctly\n");
    exit(125);
  }
  if (memcmp(data, retrieved_data, len * eleSize) != 0){
    fprintf(stderr, "error: failed to write data correctly\n");
    exit(125);
  }
  free(retrieved_data);
}
