void opt_show_option( char* name, char* desc );

void opt_show_option_extended( char* desc );

int opt_find_option( int argc, char **argv, const char *option );

int opt_read_int( int argc, char **argv, const char *option, int default_value );

double opt_read_double( int argc, char **argv, const char *option, double default_value );

float opt_read_float( int argc, char **argv, const char *option, float default_value );

char *opt_read_string( int argc, char **argv, const char *option, char *default_value );

int opt_read_irange( int argc, char **argv, const char *option,
    int* pstart, int* pstop, int* pstep );
