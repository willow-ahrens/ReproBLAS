ROOT_DIR=../..
include $(ROOT_DIR)/Make.local

all: reproBLAS

indexed:
	@cd ../indexed && $(MAKE) --no-print-directory

RBLAS1_OBJ = rdzasum.o rzsum.o rdznrm2.o rzdot.o \
	rcdot.o rcsum.o rscasum.o rscnrm2.o \
	rdsum.o rdasum.o rdnrm2.o rddot.o \
	rssum.o rsasum.o rsnrm2.o rsdot.o

$(MPOBJ) : prblas.h rblas1.h 
.PHONY: indexed

reproBLAS: indexed $(ROOT_DIR)/$(BUILD)/libs/libreproblas.a
$(ROOT_DIR)/$(BUILD)/libs/libreproblas.a: $(RBLAS1_OBJ)
	ar -crs $(ROOT_DIR)/$(BUILD)/libs/libreproblas.a $(RBLAS1_OBJ)

clean:
	-rm -rf *.o
	cd $(ROOT_DIR)/$(BUILD)/libs && rm -f *.a

