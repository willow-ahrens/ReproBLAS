#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "test_opt.h"

void opt_show_option(char* name, char* desc) {
  printf("  %8s : %s\n", name, desc);
}

void opt_show_option_extended(char* desc) {
  printf("           : %s\n", desc);
}

int opt_find_option( int argc, char **argv, const char *option )
{
  int i;
  for( i = 1; i < argc; i++ )
    if( strcmp( argv[i], option ) == 0 )
      return i;
  return -1;
}

int opt_read_int( int argc, char **argv, const char *option, int default_value )
{
  int iplace = opt_find_option( argc, argv, option );
  if( iplace >= 0 && iplace < argc-1 )
    return atoi( argv[iplace+1] );
  return default_value;
}

double opt_read_double( int argc, char **argv, const char *option, double default_value )
{
  int iplace = opt_find_option( argc, argv, option );
  if( iplace >= 0 && iplace < argc-1 )
    return atof( argv[iplace+1] );
  return default_value;
}

float opt_read_float( int argc, char **argv, const char *option, float default_value )
{
  int iplace = opt_find_option( argc, argv, option );
  if( iplace >= 0 && iplace < argc-1 )
    return atof( argv[iplace+1] );
  return default_value;
}

char *opt_read_string( int argc, char **argv, const char *option, char *default_value )
{
  int iplace = opt_find_option( argc, argv, option );
  if( iplace >= 0 && iplace < argc-1 )
    return argv[iplace+1];
  return default_value;
}

int opt_read_irange(int argc, char **argv, const char *option,
    int* pstart, int* pstop, int* pstep) {

  int nbp = opt_find_option(argc, argv, option);
  if (nbp < 0 || nbp + 1 >= argc) {
    return  -1;
  }
  int start, stop, step;
  char* nbs = argv[nbp + 1];
  char* sep = NULL;
  sep = strchr(nbs, ':');
  if (sep == NULL) {
    start = stop = atoi(nbs);
    step = 1;
  }
  else {
    sep[0] = 0;
    start = atoi(nbs);
    nbs = sep + 1;
    sep = strchr(nbs, ':');
    if (sep == NULL) {
      stop = atoi(nbs);
      step = 1;
    }
    else {
      sep[0] = 0;
      step = atoi(nbs);
      nbs = sep + 1;
      stop = atoi(nbs);
    }
  }
  if (stop < start) stop = start;
  if (step < 1) step = 1;
  *pstart = start;
  *pstop = stop;
  *pstep = step;
  return 0;
}
