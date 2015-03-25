ROOT_DIR=../..
include $(ROOT_DIR)/Make.local

MPOBJ = MPI_Indexed.o MPI_dIndexed.o MPI_sIndexed.o

libs := 
ifneq ($(MPIC_COMPILER),)
libs += MPI_Indexed
else
endif

all : $(libs)
	@:


.PHONY: MPI_Indexed

%o:%c
	$(C_COMPILER) $(FPFLAGS) $(CFLAGS) -c $< -o $@

$(MPOBJ): %.o : %.c
	$(MPIC_COMPILER) $(CFLAGS) -c $< -o $@

MPI_Indexed: $(ROOT_DIR)/$(BUILD)/libs/libmpiindexed.a
	@cp MPI_Indexed.h $(ROOT_DIR)/$(BUILD)/include/

$(ROOT_DIR)/$(BUILD)/libs/libmpiindexed.a: $(MPOBJ)
	-ar -crs $(ROOT_DIR)/$(BUILD)/libs/libmpiindexed.a $(MPOBJ)

clean:
	rm -rf *.o *.s
