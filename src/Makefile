ROOT_DIR=..
include $(ROOT_DIR)/Make.local

all: LIB

libs = Indexed reproBLAS IndexedBLAS
ifneq ($(MPIC_COMPILER),)
libs += MPI_Indexed MPI_reproBLAS
endif

LIB : $(libs)
	@cd Indexed && $(MAKE) --no-print-directory

Indexed:
	cd Indexed && $(MAKE) Indexed

MPI_Indexed: Indexed
	cd MPI_Indexed && $(MAKE) MPI_Indexed

IndexedBLAS: Indexed
	cd IndexedBLAS && $(MAKE) IndexedBLAS

reproBLAS: IndexedBLAS Indexed
	cd reproBLAS && $(MAKE) reproBLAS

MPI_reproBLAS: IndexedBLAS MPI_Indexed Indexed
	cd MPI_reproBLAS && $(MAKE) MPI_reproBLAS

clean:
	-rm -rf *.o
	-rm -f $(GENERATED)
	cd Indexed && make clean
	cd MPI_Indexed && make clean
	cd IndexedBLAS && make clean
	cd MPI_reproBLAS && make clean
	cd reproBLAS && make clean
	cd $(ROOT_DIR)/$(BUILD)/libs && rm -f *.a

