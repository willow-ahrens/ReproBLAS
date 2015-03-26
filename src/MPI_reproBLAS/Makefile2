ROOT_DIR=../..
include $(ROOT_DIR)/Make.local

all: LIB

ifneq ($(MPIC_COMPILER),)
libs += MPI_reproBLAS
endif

LIB : $(libs)
	@cd ../indexed && $(MAKE) --no-print-directory

MPOBJ = prdasum.o  prdsum.o prddot.o prdnrm2.o \
	prsasum.o  prssum.o prsdot.o prsnrm2.o \
	prdzasum.o przsum.o przdot.o prdznrm2.o \
	prscasum.o prcsum.o prcdot.o prscnrm2.o

$(MPOBJ): %.o : %.c
	$(MPIC_COMPILER) $(CFLAGS) -c $< -o $@

.PHONY: indexed MPI_indexed

indexed:
	cd ../indexed && $(MAKE) indexed

MPI_indexed: indexed
	cd ../MPI_indexed && $(MAKE) MPI_indexed

indexedBLAS: indexed
	cd ../indexedBLAS && $(MAKE) indexedBLAS

MPI_reproBLAS: MPI_indexed indexed indexedBLAS $(ROOT_DIR)/$(BUILD)/libs/libmpireproblas.a
$(ROOT_DIR)/$(BUILD)/libs/libmpireproblas.a: $(MPOBJ)
	-ar -crs $(ROOT_DIR)/$(BUILD)/libs/libmpireproblas.a $(MPOBJ)

clean:
	-rm -rf *.o

