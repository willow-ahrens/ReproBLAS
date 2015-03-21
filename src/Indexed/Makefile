ROOT_DIR=../..
include $(ROOT_DIR)/Make.local

OBJ = dIConv.o dIndexed.o dIAdd.o dIRenorm.o dIUpdate.o \
	sIndexed.o sIAdd.o sIConv.o sIRenorm.o sIUpdate.o \

MPOBJ = MPIndexedFP.o dIndexedMPI.o sIndexedMPI.o

libs := Indexed
ifneq ($(MPIC_COMPILER),)
libs += MPIndexedFP
else
endif

all : $(libs)
	@:


.PHONY: Indexed MPIndexedFP

%o:%c
	$(C_COMPILER) $(FPFLAGS) $(CFLAGS) -c $< -o $@

$(MPOBJ): %.o : %.c
	$(MPIC_COMPILER) $(CFLAGS) -c $< -o $@

IndexedFP_mod = types.mod dIndexed.mod sIndexed.mod \
	dIWrapper.mod sIWrapper.mod

$(Indexed_mod) :
	cat $< >> IndexedFP_tmp.h

Indexed_header:
	@cp Indexed.h $(ROOT_DIR)/$(BUILD)/include/IndexedFP.h

Indexed: $(ROOT_DIR)/$(BUILD)/libs/libindexedfp.a Indexed_header

$(ROOT_DIR)/$(BUILD)/libs/libindexedfp.a: $(OBJ)
	ar -crs $(ROOT_DIR)/$(BUILD)/libs/libindexedfp.a $(OBJ)

MPIndexedFP: Indexed $(ROOT_DIR)/$(BUILD)/libs/libmpindexedfp.a
	@cp MPIndexedFP.h $(ROOT_DIR)/$(BUILD)/include/

$(ROOT_DIR)/$(BUILD)/libs/libmpindexedfp.a: $(MPOBJ)
	-ar -crs $(ROOT_DIR)/$(BUILD)/libs/libmpindexedfp.a $(MPOBJ)

clean:
	rm -rf *.o *.s

