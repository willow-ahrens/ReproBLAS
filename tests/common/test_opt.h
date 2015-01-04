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
  opt_type type;
  char     short_name;
  char     *long_name;
  char     *help;
} opt_option_header;

typedef struct opt_option_flag{
  opt_option_header header;
  int               exists;
} opt_option_flag;

typedef struct opt_option_int{
  opt_option_header header;
  int               required;
  int               min;
  int               max;
  int               value;
} opt_option_int;

typedef struct opt_option_string{
  opt_option_header header;
  int               required;
  char              *value;
} opt_option_string;

typedef struct opt_option_double{
  opt_option_header header;
  int               required;
  double            min;
  double            max;
  double            value;
} opt_option_double;

typedef struct opt_option_float{
  opt_option_header header;
  int               required;
  float             min;
  float             max;
  float             value;
} opt_option_float;

typedef struct opt_option_named{
  opt_option_header header;
  int               required;
  char              **names;
  char              **descs;
  int               n_names;
  int               value;
} opt_option_named;

typedef union opt_option{
  opt_type          type;
  opt_option_header header;
  opt_option_flag   _flag;
  opt_option_int    _int;
  opt_option_string _string;
  opt_option_double _double;
  opt_option_float  _float;
  opt_option_named  _named;
} opt_option;

void opt_show_option(opt_option option);

void opt_eval_option(int argc, char **argv, opt_option *option);

#endif
