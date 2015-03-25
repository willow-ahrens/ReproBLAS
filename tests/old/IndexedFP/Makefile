include ../../make.inc

TESTS = tst1

all : $(TESTS)

RBLAS   = ../../build/librblas.a
COMMON = ../debug.o ../dgenvec.o

TESTOBJ := $(patsubst %,%.o,$(TESTS))

${TESTOBJ}: %.o: %.c ../../config.h
	$(MPCC) $(CFLAGS) -c $< -o $@

$(TESTS): $(COMMON) $(TESTOBJ) $(RBLAS)
	$(MPCC) -o $(patsubst %,run-%,$@) $(COMMON) $@.o $(BLAS) $(RBLAS) $(MPI_LIB) $(LDFLAGS)

clean:
	rm -rf *.o run-test* run-rep*

