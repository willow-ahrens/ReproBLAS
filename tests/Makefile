ROOT_DIR=..
include $(ROOT_DIR)/Make.local

COMMON = dgenvec.o zgenvec.o sgenvec.o cgenvec.o debug.o

MPITESTS =
TESTS = indexed

ifneq ($(strip $(MPIC_COMPILER)),)
MPTESTS += prdblas1 prdblas1_break prsblas1 prsblas1_break przblas1 prcblas1
endif

TESTTARGET := $(patsubst %,run-%,$(TESTS))
TESTOBJ := $(patsubst %,%.o,$(TESTS))

MPTESTTARGET := $(patsubst %,run-%,$(MPTESTS))
MPTESTOBJ := $(patsubst %,%.o,$(MPTESTS))

all : $(MPTESTTARGET) $(TESTTARGET)

RBLAS   = ../$(BUILD)/libs/libreproblas_seq.a ../$(BUILD)/libs/libiblas.a ../$(BUILD)/libs/libindexedfp.a
PRBLAS   = ../$(BUILD)/libs/libreproblas_mpi.a ../$(BUILD)/libs/libmpindexedfp.a

$(COMMON) : %.o :%.c
	@$(C_COMPILER) -c -o $@ $<

CFLAGS += -I$(ROOT_DIR)/$(BUILD)/include/

${MPTESTOBJ}: %.o: %.c debug.h
	$(MPIC_COMPILER) $(CFLAGS) -c $< -o $@

$(MPTESTTARGET): run-%: %.o $(COMMON) $(PRBLAS) $(RBLAS)
	$(MPILDC) -o $@ $(COMMON) $< $(PRBLAS) $(RBLAS) $(LDFLAGS)

$(TESTTARGET): run-%: %.o $(COMMON) $(RBLAS)
	$(LDC) -o $@ $(COMMON) $< $(RBLAS) $(LDFLAGS)

clean:
	rm -rf *.o run-*
	rm -rf *.o
	rm -rf run-*

include check.make

include accuracy.make

include benchmark.make
