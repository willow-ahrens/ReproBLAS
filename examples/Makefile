ROOT_DIR=..
include $(ROOT_DIR)/Make.local

examples = example1 example1s example1z example1c
omp_examples = 
mpi_examples = 

#OMP_FLAG = -fopenmp

ifneq (,$(OMP_FLAG))
omp_examples += example1_omp
endif

ifneq (,$(MPIC_COMPILER))
mpi_examples += example1_mpi 
endif

all: $(examples) $(omp_examples) $(mpi_examples)

OBJS   = $(patsubst %,%.o,$(examples))
TARGET = $(patsubst %,%.o,$(examples))

RBLAS   = ../$(BUILD)/libs/libreproblas_seq.a ../$(BUILD)/libs/libiblas.a ../$(BUILD)/libs/libindexedfp.a
PRBLAS   = ../$(BUILD)/libs/libreproblas_mpi.a ../$(BUILD)/libs/libmpindexedfp.a

INCR = -I../$(BUILD)/include/
LDFLAGS += -L$(realpath ../build/libs) -lm

%.o :%.c
	$(C_COMPILER) $(INCR) $(CFLAGS) -c -o $@ $<

$(examples): % :%.o
	$(LDC) -o run-$@ $(patsubst %,%.o,$@) $(LDFLAGS) -lreproblas_seq -liblas -lindexedfp

$(omp_examples): 
	$(C_COMPILER) $(CFLAGS) $(INCR) $(OMP_FLAG) -c $(patsubst %,%.c,$@)
	$(LDC) $(OMP_FLAG) -o run-$@ $(patsubst %,%.o,$@) $(LDFLAGS) -lindexedfp

$(mpi_examples): 
	$(MPIC_COMPILER) $(CFLAGS) $(INCR) -c $(patsubst %,%.c,$@)
	$(MPILDC)  -o run-$@ $(patsubst %,%.o,$@) $(LDFLAGS) -lreproblas_mpi -lreproblas_seq -liblas -lmpindexedfp -lindexedfp 

clean:
	rm -rf *.o run-*
	rm -rf *.o
	rm -rf run-*

