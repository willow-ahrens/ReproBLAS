ROOT_DIR=../..
include $(ROOT_DIR)/Make.local

examples = funarc funarc_omp funarc_I funarc_I2 funarc_I_vec
omp_examples = 
mpi_examples = 

#OMP_FLAG = -fopenmp

ifneq (,$(OMP_FLAG))
omp_examples += funarc_omp
endif

ifneq (,$(MPIC_COMPILER))
mpi_examples += funarc_mpi funarc_I_mpi funarc_I2_mpi
endif

all: $(examples) $(omp_examples) $(mpi_examples)

OBJS   = $(patsubst %,%.o,$(examples))
TARGET = $(patsubst %,%.o,$(examples))

INCR = -I$(ROOT_DIR)/$(BUILD)/include/
LDFLAGS += -L$(realpath $(ROOT_DIR)/$(BUILD)/libs) -lm

%.o :%.c
	$(C_COMPILER) $(INCR) $(CFLAGS) -c -o $@ $<

$(examples): % :%.o
	$(LDC) -o run-$@ $(patsubst %,%.o,$@) $(LDFLAGS) -lreproblas_mpi -liblas -lindexedfp

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

