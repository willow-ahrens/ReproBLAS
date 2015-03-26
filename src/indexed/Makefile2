ROOT_DIR=../..
include $(ROOT_DIR)/Make.local

OBJ = dIConv.o dIndexed.o dIAdd.o dIRenorm.o dIUpdate.o \
	sIndexed.o sIAdd.o sIConv.o sIRenorm.o sIUpdate.o \

libs := indexed

all : $(libs)
	@:


.PHONY: indexed

%o:%c
	$(C_COMPILER) $(FPFLAGS) $(CFLAGS) -c $< -o $@

indexed: $(ROOT_DIR)/$(BUILD)/libs/libindexed.a

$(ROOT_DIR)/$(BUILD)/libs/libindexed.a: $(OBJ)
	ar -crs $(ROOT_DIR)/$(BUILD)/libs/libindexed.a $(OBJ)

clean:
	rm -rf *.o *.s

