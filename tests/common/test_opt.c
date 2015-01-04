#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "test_opt.h"

void opt_show_option(opt_option option) {
  opt_fprintf_option(stdout, opt_option option);
}
static void opt_fprintf_option(FILE *f, opt_option option) {
  char* token;
  token = strtok(option.header.help, "\n");
  if(token){
    fprintf(f, "-%c, --%8s: %64s\n", option.header.shortname, option.header.longname, token);
    while(token){
      fprintf(f, "                %64s\n", token);
      token = strtok(NULL, "\n");
    }
  }else{
    fprintf(f, "-%c, --%8s:", option.header.shortname, option.header.longname);
  }
  token = strtok(NULL, "\n");
  fprintf(f, "                (");
  if(option.header.required){
    fprintf(f, "required ")
  }
  switch(option.type){
    case opt_flag:
      fprintf(f, "flag)\n");
      break;
    case opt_int:
      fprintf(f, "int) values=[%d, %d] default=%d\n", option.i.min, option.i.max, option.i.default);
      break;
    case opt_string:
      fprintf(f, "string) default=%d\n", option.s.default);
      break;
    case opt_double:
      fprintf(f, "double) values=[%f, %f] default=%f\n", option.d.min, option.d.max, option.d.default);
      break;
    case opt_float:
      fprintf(f, "float) values=[%f, %f] default=%f\n", option.f.min, option.f.max, option.f.default);
      break;
    case opt_named:
      fprintf(f, "enum) values=[\n");
      {
        int i;
        char[63] namedesc;
        for(i = 0; i < option.n.n_names - 1; i++){
          snprintf(namedesc, 63 * sizeof(char), "%s: %s", option.n.names[i], option.n.descs[i]);
          fprintf(f, "                  %s\n", namedesc);
        }
        snprintf(namedesc, 62 * sizeof(char), "%s: %s", option.n.names[option.n.n_names - 1], option.n.descs[option.n.n_names - 1]);
        fprintf(f, "                  %s]\n", namedesc);
        fprintf(f, "                  default=%54s\n", option.n.names[option.n.default]);
      }
      break;
  }
}

static int opt_find_flag(int argc, char **argv, char *flag)
{
  int i;
  for( i = 1; i < argc; i++ )
    if( strcmp( argv[i], option ) == 0 )
      return i;
  return -1;
}

void opt_eval_option(int argc, char **argv, opt_option *option);
  char flag[11];
  int i_flag;
  snfprintf(f, flag, 11 * sizeof(str), "-%c", option->header.shortname);
  i_flag = opt_find_flag(argc, argv, flag);
  if(i_flag < 0){
    snfprintf(f, flag, 11 * sizeof(str), "-%s", option->header.longname);
    i_flag = opt_find_flag(argc, argv, flag);
  }
  if(i_flag < 0 && option->header.required){
    fprintf(stderr, "error: required argument not found\n");
    opt_fprintf_option(stderr, *option);
    exit(125);
  }
  switch(option.type){
    case opt_flag:
      option.f.exists = i_flag >= 0;
      break;
    case opt_int:
      if(i_flag >= 0){
        if(i_flag + 1 >= argc){
          fprintf(stderr, "error: option takes an argument\n", argv[i_flag + 1]);
          opt_fprintf_option(stderr, *option);
        }
        char *end;
        option.i.value = (int)strtol(argv[i_flag + 1], &end, 0);
        if(end != argv[i_flag + 1] + strlen(argv[i_flag + 1])){
          fprintf(stderr, "error: invalid integer format %s\n", argv[i_flag + 1]);
          opt_fprintf_option(stderr, *option);
        }
        if(option.i.value > option.i.max || option.i.value < option.i.min){
          fprintf(stderr, "error: integer out of bounds %d\n", option.i.value);
        }
      }else{
        option.i.value = option.i.default;
      }
      break;
    case opt_string:
      if(i_flag >= 0){
        if(i_flag + 1 >= argc){
          fprintf(stderr, "error: option takes an argument\n", argv[i_flag + 1]);
          opt_fprintf_option(stderr, *option);
        }
        option.s.value = argv[i_flag + 1];
      }else{
        option.s.value = option.s.default;
      }
      break;
    case opt_double:
      if(i_flag >= 0){
        if(i_flag + 1 >= argc){
          fprintf(stderr, "error: option takes an argument\n", argv[i_flag + 1]);
          opt_fprintf_option(stderr, *option);
        }
        char *end;
        option.d.value = strtod(argv[i_flag + 1], &end);
        if(end != argv[i_flag + 1] + strlen(argv[i_flag + 1])){
          fprintf(stderr, "error: invalid double format %s\n", argv[i_flag + 1]);
          opt_fprintf_option(stderr, *option);
        }
        if(option.d.value > option.d.max || option.d.value < option.d.min){
          fprintf(stderr, "error: integer out of bounds %f\n", option.d.value);
        }
      }else{
        option.d.value = option.d.default;
      }
      break;
    case opt_float:
      if(i_flag >= 0){
        if(i_flag + 1 >= argc){
          fprintf(stderr, "error: option takes an argument\n", argv[i_flag + 1]);
          opt_fprintf_option(stderr, *option);
        }
        char *end;
        option.f.value = (float)strtod(argv[i_flag + 1], &end);
        if(end != argv[i_flag + 1] + strlen(argv[i_flag + 1])){
          fprintf(stderr, "error: invalid float format %s\n", argv[i_flag + 1]);
          opt_fprintf_option(stderr, *option);
        }
        if(option.f.value > option.f.max || option.f.value < option.f.min){
          fprintf(stderr, "error: integer out of bounds %f\n", option.f.value);
        }
      }else{
        option.d.value = option.d.default;
      }
      break;
    case opt_named:
      if(i_flag >= 0){
        if(i_flag + 1 >= argc){
          fprintf(stderr, "error: option takes an argument\n", argv[i_flag + 1]);
          opt_fprintf_option(stderr, *option);
        }
        int i;
        for(i = 0; i < option.n.n_names; i++){
          if(strcmp(option.n.names[i], argv[i_flag + 1]) == 0){
            option.n.value = i;
            return;
          }
        }
        fprintf(stderr, "error: invalid argument %s\n", argv[i_flag + 1]);
        opt_fprintf_option(stderr, *option);
      }else{
        option.n.value = option.n.default;
      }
      break;
  }
}
