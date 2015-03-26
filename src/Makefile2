ROOT_DIR=..
include $(ROOT_DIR)/Make.local

all: LIB

libs = indexed reproBLAS indexedBLAS
ifneq ($(MPIC_COMPILER),)
libs += MPI_indexed MPI_reproBLAS
endif

LIB : $(libs)
	@cd indexed && $(MAKE) --no-print-directory

indexed:
	cd indexed && $(MAKE) indexed

MPI_indexed: indexed
	cd MPI_indexed && $(MAKE) MPI_indexed

indexedBLAS: indexed
	cd indexedBLAS && $(MAKE) indexedBLAS

reproBLAS: indexedBLAS indexed
	cd reproBLAS && $(MAKE) reproBLAS

MPI_reproBLAS: indexedBLAS MPI_indexed indexed
	cd MPI_reproBLAS && $(MAKE) MPI_reproBLAS

clean:
	-rm -rf *.o
	-rm -f $(GENERATED)
	cd indexed && make clean
	cd MPI_indexed && make clean
	cd indexedBLAS && make clean
	cd MPI_reproBLAS && make clean
	cd reproBLAS && make clean
	cd $(ROOT_DIR)/$(BUILD)/libs && rm -f *.a

