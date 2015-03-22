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

$(ROOT_DIR)/$(BUILD)/libs/libindexedfp.a:
	cd ../Indexed
	$(Make) Indexed

MPI_Indexed: $(ROOT_DIR)/$(BUILD)/libs/libmpindexedfp.a
	@cp MPI_Indexed.h $(ROOT_DIR)/$(BUILD)/include/

$(ROOT_DIR)/$(BUILD)/libs/libmpindexedfp.a: $(MPOBJ)
	-ar -crs $(ROOT_DIR)/$(BUILD)/libs/libmpindexedfp.a $(MPOBJ)

clean:
	rm -rf *.o *.s
