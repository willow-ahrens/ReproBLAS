#include <stdio.h>
#include "test_opt.h"

#include "test_header.h"

const char *file_name(int argc, char** argv);
int file_test(int argc, char** argv, char *fname);

const char* name(int argc, char** argv){
  opt_option input_file;
  static char name_buffer[MAX_LINE];

  input_file.header.type       = opt_string;
  input_file.header.short_name = 'i';
  input_file.header.long_name  = "input";
  input_file.header.help       = "input file name";
  input_file._string.required  = 1;
  input_file._string.value     = "";

  if(help._flag.exists){
    opt_show_option(input_file);
    file_name(argc, argv);
    return "";
  }

  opt_eval_option(argc, argv, &input_file);
  snprintf(name_buffer, MAX_LINE, "%s (%s)", file_name(argc, argv), input_file._string.value);
  return name_buffer;
}

int test(int argc, char** argv){
  opt_option input_file;
  int rc;

  input_file.header.type       = opt_string;
  input_file.header.short_name = 'i';
  input_file.header.long_name  = "input";
  input_file.header.help       = "input file name";
  input_file._string.required  = 1;
  input_file._string.value     = "";

  opt_eval_option(argc, argv, &input_file);
  rc = file_test(argc, argv, input_file._string.value);
  return rc;
}
