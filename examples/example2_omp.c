#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>

#include <IndexedFP.h>
#include <Iblas.h>

#define element(i,n) sin(2 * M_PI * (i / (double)n - 0.5))

static int initialized = 0;
static struct timeval start;

inline void start_timer( )
{
    gettimeofday( &start, NULL );
}

inline double read_timer( )
{
    struct timeval end;
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

int main (int argc, char** args) {
    int n = 100000;
    int i;
    I_double I_s;
    double s;
    int chunk = 128;
    double* v;
    int nbthreads;

    if (argc > 1) {
        nbthreads = atoi(args[1]);
        omp_set_num_threads(nbthreads);
    }
    if (argc > 2) {
        n = atoi(args[2]);
    }

    v = malloc(n * sizeof(double));
    for (i = 0; i < n; i++) {
        v[i] = element(i,n);
    }

    // multi-threaded
    start_timer();
    double t = 0;
    int niters = 0;
    while (t < 1.0) {
    niters++;
    dISetZero(I_s);
    I_double buffer[16];
    #pragma omp parallel
    {
        int me;
        int NB, n_start, n_end;
        nbthreads = omp_get_num_threads();
        me = omp_get_thread_num();
        NB = (n / nbthreads);
        NB = NB - NB % 64;
        if (NB * nbthreads < n) {
            NB += 64;
        }
        n_start = me * NB;
        n_end = me * NB + NB;
        if (n_start > n) n_start = n;
        if (n_end > n) n_end = n;
        dISetZero(buffer[me]);
	    dsumI1((n_end - n_start), (v + n_start), 1, DEFAULT_FOLD, 0, buffer[me].m, buffer[me].c);
    }
    
    for (i = 0; i < nbthreads; i++) {
        dIAdd(&I_s, buffer[i]);
    }
    
    s = Iconv2d(I_s);
    t = read_timer();
    }
	printf("%2d theads.  ", nbthreads);
    printf("Result: %g, Elapased: %g, effective performance: %.3g GFLOPs \n", s, t, (n / (1e9 * t)) * niters);
    free(v);
}
