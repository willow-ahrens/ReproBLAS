#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "test_opt.h"

static void opt_fprintf_option(FILE *f, opt_option option) {
  char* token;
  token = strtok(option.header.help, "\n");
  if(token){
    fprintf(f, "-%c, --%8s: %64s\n", option.header.short_name, option.header.long_name, token);
    while(token){
      fprintf(f, "                %64s\n", token);
      token = strtok(NULL, "\n");
    }
  }else{
    fprintf(f, "-%c, --%8s:", option.header.short_name, option.header.long_name);
  }
  token = strtok(NULL, "\n");
  fprintf(f, "                (");
  switch(option.type){
    case opt_flag:
      fprintf(f, "flag)\n");
      break;
    case opt_int:
      if(option._int.required){
        fprintf(f, "required ");
      }
      fprintf(f, "int) values=[%d, %d] default=%d\n", option._int.min, option._int.max, option._int.value);
      break;
    case opt_string:
      if(option._string.required){
        fprintf(f, "required ");
      }
      fprintf(f, "string) default=%s\n", option._string.value);
      break;
    case opt_double:
      if(option._double.required){
        fprintf(f, "required ");
      }
      fprintf(f, "double) values=[%f, %f] default=%f\n", option._double.min, option._double.max, option._double.value);
      break;
    case opt_float:
      if(option._float.required){
        fprintf(f, "required ");
      }
      fprintf(f, "float) values=[%f, %f] default=%f\n", option._float.min, option._float.max, option._float.value);
      break;
    case opt_named:
      if(option._named.required){
        fprintf(f, "required ");
      }
      fprintf(f, "enum) values=[\n");
      {
        int i;
        char namedesc[63];
        for(i = 0; i < option._named.n_names - 1; i++){
          snprintf(namedesc, 63 * sizeof(char), "%s: %s", option._named.names[i], option._named.descs[i]);
          fprintf(f, "                  %s\n", namedesc);
        }
        snprintf(namedesc, 62 * sizeof(char), "%s: %s", option._named.names[option._named.n_names - 1], option._named.descs[option._named.n_names - 1]);
        fprintf(f, "                  %s]\n", namedesc);
        fprintf(f, "                  default=%54s\n", option._named.names[option._named.value]);
      }
      break;
  }
}

static int opt_find_flag(int argc, char **argv, char *flag)
{
  int i;
  for( i = 1; i < argc; i++ )
    if( strcmp( argv[i], flag ) == 0 )
      return i;
  return -1;
}

void opt_show_option(opt_option option) {
  opt_fprintf_option(stdout, option);
}

void opt_eval_option(int argc, char **argv, opt_option *option){
  char flag[11];
  int i_flag;
  snprintf(flag, 11 * sizeof(char), "-%c", option->header.short_name);
  i_flag = opt_find_flag(argc, argv, flag);
  if(i_flag < 0){
    snprintf(flag, 11 * sizeof(char), "--%s", option->header.long_name);
    i_flag = opt_find_flag(argc, argv, flag);
  }
  switch(option->type){
    case opt_flag:
      option->_flag.exists = i_flag >= 0;
      break;
    case opt_int:
      if(i_flag >= 0){
        if(i_flag + 1 >= argc){
          fprintf(stderr, "error: option takes an argument\n", argv[i_flag + 1]);
          opt_fprintf_option(stderr, *option);
        }
        char *end;
        option->_int.value = (int)strtol(argv[i_flag + 1], &end, 0);
        if(end != argv[i_flag + 1] + strlen(argv[i_flag + 1])){
          fprintf(stderr, "error: invalid integer format %s\n", argv[i_flag + 1]);
          opt_fprintf_option(stderr, *option);
        }
        if(option->_int.value > option->_int.max || option->_int.value < option->_int.min){
          fprintf(stderr, "error: integer out of bounds %d\n", option->_int.value);
        }
      }else{
        if(i_flag < 0 && option->_int.required){
          fprintf(stderr, "error: required argument not found\n");
          opt_fprintf_option(stderr, *option);
          exit(125);
        }
      }
      break;
    case opt_string:
      if(i_flag >= 0){
        if(i_flag + 1 >= argc){
          fprintf(stderr, "error: option takes an argument\n", argv[i_flag + 1]);
          opt_fprintf_option(stderr, *option);
        }
        option->_string.value = argv[i_flag + 1];
      }else{
        if(i_flag < 0 && option->_string.required){
          fprintf(stderr, "error: required argument not found\n");
          opt_fprintf_option(stderr, *option);
          exit(125);
        }
      }
      break;
    case opt_double:
      if(i_flag >= 0){
        if(i_flag + 1 >= argc){
          fprintf(stderr, "error: option takes an argument\n", argv[i_flag + 1]);
          opt_fprintf_option(stderr, *option);
        }
        char *end;
        option->_double.value = strtod(argv[i_flag + 1], &end);
        if(end != argv[i_flag + 1] + strlen(argv[i_flag + 1])){
          fprintf(stderr, "error: invalid double format %s\n", argv[i_flag + 1]);
          opt_fprintf_option(stderr, *option);
        }
        if(option->_double.value > option->_double.max || option->_double.value < option->_double.min){
          fprintf(stderr, "error: integer out of bounds %f\n", option->_double.value);
        }
      }else{
        if(i_flag < 0 && option->_double.required){
          fprintf(stderr, "error: required argument not found\n");
          opt_fprintf_option(stderr, *option);
          exit(125);
        }
      }
      break;
    case opt_float:
      if(i_flag >= 0){
        if(i_flag + 1 >= argc){
          fprintf(stderr, "error: option takes an argument\n", argv[i_flag + 1]);
          opt_fprintf_option(stderr, *option);
        }
        char *end;
        option->_float.value = (float)strtod(argv[i_flag + 1], &end);
        if(end != argv[i_flag + 1] + strlen(argv[i_flag + 1])){
          fprintf(stderr, "error: invalid float format %s\n", argv[i_flag + 1]);
          opt_fprintf_option(stderr, *option);
        }
        if(option->_float.value > option->_float.max || option->_float.value < option->_float.min){
          fprintf(stderr, "error: integer out of bounds %f\n", option->_float.value);
        }
      }else{
        if(i_flag < 0 && option->_float.required){
          fprintf(stderr, "error: required argument not found\n");
          opt_fprintf_option(stderr, *option);
          exit(125);
        }
      }
      break;
    case opt_named:
      if(i_flag >= 0){
        if(i_flag + 1 >= argc){
          fprintf(stderr, "error: option takes an argument\n", argv[i_flag + 1]);
          opt_fprintf_option(stderr, *option);
        }
        int i;
        for(i = 0; i < option->_named.n_names; i++){
          if(strcmp(option->_named.names[i], argv[i_flag + 1]) == 0){
            option->_named.value = i;
            return;
          }
        }
        fprintf(stderr, "error: invalid argument %s\n", argv[i_flag + 1]);
        opt_fprintf_option(stderr, *option);
      }else{
        if(i_flag < 0 && option->_named.required){
          fprintf(stderr, "error: required argument not found\n");
          opt_fprintf_option(stderr, *option);
          exit(125);
        }
      }
      break;
  }
}
