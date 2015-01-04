#include <stdio.h>
#include <limits.h>
#include "test_limits.h"
#include "test_opt.h"

opt_option help;

const char *name(int argc, char** argv);
int test(int argc, char** argv);

int main(int argc, char** argv){
  opt_option print;

  help.header.type       = opt_flag;
  help.header.short_name = 'h';
  help.header.long_name  = "help";
  help.header.help       = "show help";

  print.header.type       = opt_flag;
  print.header.short_name = 'p';
  print.header.long_name  = "print";
  print.header.help       = "print parameters but do not execute";

  opt_eval_option(argc, argv, &help);
  if(help._flag.exists){
    opt_show_option(help);
    opt_show_option(print);
    name(argc, argv);
    return 0;
  }

  opt_eval_option(argc, argv, &print);
  if(print._flag.exists){
    printf("%s\n", name(argc, argv));
    return 0;
  } else {
    return test(argc, argv);
  }
}
