#include <complex.h>
#include <stdio.h>


int main(int argc, char** argv){
  double complex foo;
  double complex bar;
  ((double*)&foo)[0] = 1.0/0.0;
  ((double*)&foo)[1] = 0.0;
  ((double*)&bar)[0] = 1.0/0.0;
  ((double*)&bar)[1] = 1.0/0.0;
  double complex qux = foo * bar;
  printf("(%g + %gi)(%g + %gi) = %g + %gi\n", creal(foo), cimag(foo), creal(bar), cimag(bar), creal(qux), cimag(qux));
}
