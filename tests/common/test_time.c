#include <sys/time.h>
#include <stddef.h>
static int initialized = 0;
static struct timeval start;
static double time = 0.0;

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
  time -= read_clock();
}

void time_toc() {
  time += read_clock();
}

double time_read() {
  return time;
}

void time_reset() {
  time = 0.0;
}
