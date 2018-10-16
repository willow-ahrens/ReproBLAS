CC = gcc
LDC = gcc

#MPICC = mpicc
MPILDC = $(MPICC)

#OMP_FLAG = -fopenmp

examples = sum_sine
omp_examples = 
mpi_examples = 

CFLAGS += -O3
LDFLAGS += -lm

ifneq ($(OMP_FLAG),)
	omp_examples += 
endif

ifneq ($(MPICC),)
	mpi_examples += MPI_sum_sine
endif

.PHONY: all clean

all: $(examples) $(omp_examples) $(mpi_examples)

OBJS   = $(patsubst %,%.o,$(examples))
TARGET = $(patsubst %,%.o,$(examples))

%.o :%.c
	$(CC) $(CFLAGS) -c -o $@ $<

$(examples): % : %.o
	$(LDC) -o $@ $(patsubst %,%.o,$@) $(LDFLAGS) -lbinned -lbinnedblas -lreproblas

$(omp_examples): % : %.o
	$(CC) $(CFLAGS) $(OMP_FLAG) -c $(patsubst %,%.c,$@)
	$(LDC) $(OMP_FLAG) -o run-$@ $(patsubst %,%.o,$@) $(LDFLAGS) -lbinned -lbinnedblas -lreproblas

$(mpi_examples): % : %.o
	$(MPICC) $(CFLAGS) $(INCR) -c $(patsubst %,%.c,$@)
	$(MPILDC) -o $@ $(patsubst %,%.o,$@) $(LDFLAGS) -lbinned -lbinnedblas -lreproblas -lbinnedMPI

clean:
	rm -rf *.o
	rm -rf $(examples) $(omp_examples) $(mpi_examples)
