ROOT_DIR=../..
include $(ROOT_DIR)/Make.local

OBJ = dIConv.o dIndexed.o dIAdd.o dIRenorm.o dIUpdate.o \
	sIndexed.o sIAdd.o sIConv.o sIRenorm.o sIUpdate.o \

MPOBJ = MPIndexedFP.o dIndexedMPI.o sIndexedMPI.o

libs := IndexedFP
ifneq ($(MPIC_COMPILER),)
libs += MPIndexedFP
else
endif

all : $(libs)
	@:


.PHONY: IndexedFP MPIndexedFP

%o:%c
	$(C_COMPILER) $(FPFLAGS) $(CFLAGS) -c $< -o $@

$(MPOBJ): %.o : %.c
	$(MPIC_COMPILER) $(CFLAGS) -c $< -o $@

IndexedFP_mod = config.inc types.mod dIndexed.mod sIndexed.mod \
	dIWrapper.mod sIWrapper.mod

$(Indexed_mod) :
	cat $< >> IndexedFP_tmp.h

IndexedFP_header:
	@echo '#ifndef _INDEXED_FP__H' >  IndexedFP_tmp.h
	@echo '#define _INDEXED_FP__H' >> IndexedFP_tmp.h
	@$(foreach mod, $(IndexedFP_mod),cat $(mod) >> IndexedFP_tmp.h;)
	@echo '#endif' >> IndexedFP_tmp.h
	@mv IndexedFP_tmp.h $(ROOT_DIR)/$(BUILD)/include/IndexedFP.h

IndexedFP: $(ROOT_DIR)/$(BUILD)/libs/libindexedfp.a IndexedFP_header

$(ROOT_DIR)/$(BUILD)/libs/libindexedfp.a: $(OBJ)
	ar -crs $(ROOT_DIR)/$(BUILD)/libs/libindexedfp.a $(OBJ)

MPIndexedFP: IndexedFP $(ROOT_DIR)/$(BUILD)/libs/libmpindexedfp.a
	@cp MPIndexedFP.h $(ROOT_DIR)/$(BUILD)/include/

$(ROOT_DIR)/$(BUILD)/libs/libmpindexedfp.a: $(MPOBJ)
	-ar -crs $(ROOT_DIR)/$(BUILD)/libs/libmpindexedfp.a $(MPOBJ)

clean:
	rm -rf *.o *.s

