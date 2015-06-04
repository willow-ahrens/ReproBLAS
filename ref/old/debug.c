#include "debug.h"
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>

static int initialized = 0;
static struct timeval start;
static double tic_;
static double timing_breakdown_[BREAKDOWN_N];

double read_timer( )
{
  struct timeval end;
  if( initialized == 0)
    {
      gettimeofday( &start, NULL );
      initialized = 1;
    }

  gettimeofday( &end, NULL );

  return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

void show_opt(char* name, char* desc) {
	printf("  %8s : %s\n", name, desc);
}

void timing_breakdown(double* timing) {
	memcpy(timing, timing_breakdown_, BREAKDOWN_N * sizeof(double));
}

void clock_start() {
	int i;
	tic_ = read_timer();
	for (i = 0; i < BREAKDOWN_N; i++)
		timing_breakdown_[i] = 0.0;
}
void tictoc(int ind) {
	double toc_ = read_timer();
	timing_breakdown_[ind] += toc_ - tic_;
	tic_ = toc_;
}

int find_option( int argc, char **argv, const char *option )
{
  int i;
  for( i = 1; i < argc; i++ )
    if( strcmp( argv[i], option ) == 0 )
      return i;
  return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
  int iplace = find_option( argc, argv, option );
  if( iplace >= 0 && iplace < argc-1 )
    return atoi( argv[iplace+1] );
  return default_value;
}
double read_double( int argc, char **argv, const char *option, double default_value )
{
  int iplace = find_option( argc, argv, option );
  if( iplace >= 0 && iplace < argc-1 )
    return atof( argv[iplace+1] );
  return default_value;
}

int rsign() {
	double t1 = drand48();
	if (t1 > 0.5) return 1;
	else return -1;
}

float read_float( int argc, char **argv, const char *option, float default_value )
{
  int iplace = find_option( argc, argv, option );
  if( iplace >= 0 && iplace < argc-1 )
    return atof( argv[iplace+1] );
  return default_value;
}

int read_irange(int argc, char **argv, const char *option,
	int* pstart, int* pstop, int* pstep) {

	int nbp = find_option(argc, argv, option);
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

void drandm(double* A, int M, int N, int lda, double from, double to) {
	int i, j;
	for (i = 0; i < N; i++) {
		for (j = 0; j < M; j++)
			A[j] = drand48() * (to - from) + from;
		A += lda;
	}
}

double ddiff(int N, double* v1, double* v2, int* ind) {
	int i;
	double diff = 0.0;
	int mi = -1;
	double tmp;
	for (i = 0; i < N; i++) {
		tmp = (v1[i] - v2[i]);
		if (tmp < 0) tmp = -tmp;
		if (diff < tmp) {
			diff  = tmp;
			mi = i;
		}
	}

	if (ind != NULL)
		*ind = mi;

	return diff;
}

double derrs(int N, double* v1, double* v2, int* ind, double* abs) {
	int i;
	double diff = 0.0;
	double rel  = 0.0;
	int mi = -1;
	double tmp, tmpm;
	for (i = 0; i < N; i++) {
		tmp = (v1[i] - v2[i]);
		if (tmp < 0) tmp = -tmp;
		tmpm = v1[i];
		if (tmpm < 0.0) tmpm = -tmpm;
		if (tmpm < v2[i]) tmpm = v2[i];
		else if (tmpm < -v2[i]) tmpm = -v2[i];

		if (tmpm > 0.0) {
			tmpm = tmp / tmpm;
			if (rel < tmpm) {
				rel = tmpm;
				mi  = i;
			}
		}

		if (diff < tmp) {
			diff  = tmp;
		}
	}

	if (ind != NULL)
		*ind = mi;

	if (abs != NULL)
		*abs = diff;

	return rel;
}

int fileExists(const char *fname)
{
    FILE *file = fopen(fname, "r");
    if (file)
    {
        fclose(file);
        return 1;
    }
    return 0;
}

void readFile(const char* fname, int* len, void** data, int eleSize) {
	int i;
	FILE* file = fopen(fname, "rb");
	fread(len, sizeof(int), 1, file);
	*data = malloc(*len * eleSize);
	char* ptr = *data;
	for(i = 0; i < *len; i++, ptr += eleSize) {
		fread(ptr, eleSize, 1, file);
	}
	fclose(file);
}

void writeFile(const char* fname, int len, void* data, int eleSize) {
	int i;
	FILE* file = fopen(fname, "wb");
	if(file == NULL) {
		fprintf(stderr,"FAILED to open file %s to write.\n", fname);
		return;
	}
	fwrite(&len, sizeof(int), 1, file);
	char* ptr = data;
	for(i = 0; i < len; i++, ptr += eleSize) {
		fwrite(ptr, eleSize, 1, file);
	}
	fclose(file);
}

