#include <sys/time.h>
#include <stddef.h>
static int initialized = 0;
static struct timeval start;
static double tic;
static double toc;

static double read_clock( )
{
  struct timeval end;
  if( initialized == 0)
  {
    gettimeofday( &start, NULL );
    initialized = 1;
  }

  gettimeofday( &end, NULL );

  return (end.tv_sec  + 1.0e-6 * end.tv_usec) - (start.tv_sec - 1.0e-6 * start.tv_usec);
}

void time_tic() {
  tic = read_clock();
}

double time_toc() {
  toc = read_clock();
  return toc - tic;
}

double time_read() {
  return toc - tic;
}
