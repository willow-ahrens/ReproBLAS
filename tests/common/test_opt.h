#ifndef __TEST_OPT_H
#define __TEST_OPT_H

typedef enum opt_type{
  opt_flag,
  opt_int,
  opt_string,
  opt_double,
  opt_float,
  opt_named
} opt_type;

typedef struct opt_option_header{
  opt_type    type;
  char        short_name;
  const char  *long_name;
  const char  *help;
  int         required;
} opt_option_header;

typedef struct opt_option_flag{
  opt_option_header header;
  int               exists;
} opt_option_flag;

typedef struct opt_option_int{
  opt_option_header header;
  int               min;
  int               max;
  int               default;
  int               value;
} opt_option_int;

typedef struct opt_option_string{
  opt_option_header header;
  const char        *default;
  char              *value;
} opt_option_string;

typedef struct opt_option_double{
  opt_option_header header;
  double            min;
  double            max;
  double            default;
  double            value;
} opt_option_double;

typedef struct opt_option_float{
  opt_option_header header;
  float             min;
  float             max;
  float             default;
  float             value;
} opt_option_float;

typedef struct opt_option_named{
  opt_option_header header;
  char              **names;
  char              **descs;
  int               n_names;
  int               default;
  int               value;
} opt_option_named;

typedef union opt_option{
  opt_type          type;
  opt_option_header header;
  opt_option_flag   f;
  opt_option_int    i;
  opt_option_string s;
  opt_option_double d;
  opt_option_float  f;
  opt_option_named  n;
}

void opt_show_option(opt_option option);

void opt_eval_option(int argc, char **argv, opt_option *option);

#endif
