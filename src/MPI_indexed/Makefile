ROOT_DIR=../..
include $(ROOT_DIR)/Make.local

MPOBJ = MPI_indexed.o MPI_dIndexed.o MPI_sIndexed.o

libs := 
ifneq ($(MPIC_COMPILER),)
libs += MPI_indexed
else
endif

all : $(libs)
	@:


.PHONY: MPI_indexed

%o:%c
	$(C_COMPILER) $(FPFLAGS) $(CFLAGS) -c $< -o $@

$(MPOBJ): %.o : %.c
	$(MPIC_COMPILER) $(CFLAGS) -c $< -o $@

MPI_indexed: $(ROOT_DIR)/$(BUILD)/libs/libmpiindexed.a

$(ROOT_DIR)/$(BUILD)/libs/libmpiindexed.a: $(MPOBJ)
	-ar -crs $(ROOT_DIR)/$(BUILD)/libs/libmpiindexed.a $(MPOBJ)

clean:
	rm -rf *.o *.s
