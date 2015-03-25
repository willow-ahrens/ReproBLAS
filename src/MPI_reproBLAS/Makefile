ROOT_DIR=../..
include $(ROOT_DIR)/Make.local

all: LIB

ifneq ($(MPIC_COMPILER),)
libs += MPI_reproBLAS
endif

LIB : $(libs)
	@cd ../Indexed && $(MAKE) --no-print-directory

MPOBJ = prdasum.o  prdsum.o prddot.o prdnrm2.o \
	prsasum.o  prssum.o prsdot.o prsnrm2.o \
	prdzasum.o przsum.o przdot.o prdznrm2.o \
	prscasum.o prcsum.o prcdot.o prscnrm2.o

$(MPOBJ): %.o : %.c
	$(MPIC_COMPILER) $(CFLAGS) -c $< -o $@

.PHONY: Indexed MPI_Indexed

Indexed:
	cd ../Indexed && $(MAKE) Indexed

MPI_Indexed: Indexed
	cd ../MPI_Indexed && $(MAKE) MPI_Indexed

IndexedBLAS: Indexed
	cd ../IndexedBLAS && $(MAKE) IndexedBLAS

MPI_reproBLAS: MPI_Indexed Indexed IndexedBLAS $(ROOT_DIR)/$(BUILD)/libs/libmpireproblas.a
$(ROOT_DIR)/$(BUILD)/libs/libmpireproblas.a: $(MPOBJ)
	-ar -crs $(ROOT_DIR)/$(BUILD)/libs/libmpireproblas.a $(MPOBJ)

clean:
	-rm -rf *.o

