ROOT_DIR=../..
include $(ROOT_DIR)/Make.local

all: LIB

libs = IndexedBLAS 

LIB : $(libs)

DIBLAS1 = dasumI2.o dsumI2.o dnrm2I2.o ddotI2.o
SIBLAS1 = sasumI2.o ssumI2.o snrm2I2.o sdotI2.o
ZIBLAS1 = zsumI2.o dzasumI2.o dznrm2I2.o zdotcI2.o zdotuI2.o
CIBLAS1 = csumI2.o scasumI2.o scnrm2I2.o cdotcI2.o cdotuI2.o

IBLAS1  = ssumI.o sasumI.o sdotI.o snrm2I.o \
	dsumI.o dasumI.o ddotI.o dnrm2I.o	\
	zdotI.o zsumI.o dzasumI.o dznrm2I.o \
	csumI.o scasumI.o cdotI.o scnrm2I.o

SSEO = samax.o samaxm.o damax.o damaxm.o zamax.o zamaxm.o camax.o camaxm.o	\
	dIAccum.o sIAccum.o zIAccum.o cIAccum.o

CORE_OBJ = $(SSEO) $(DIBLAS1) $(SIBLAS1) $(ZIBLAS1) $(CIBLAS1) $(IBLAS1)

DIBLAS1_REF = ref/dasumI2.o ref/dsumI2.o ref/dnrm2I2.o ref/ddotI2.o
SIBLAS1_REF = ref/sasumI2.o ref/ssumI2.o ref/snrm2I2.o ref/sdotI2.o
ZIBLAS1_REF = ref/zsumI2.o ref/dzasumI2.o ref/dznrm2I2.o ref/zdotcI2.o ref/zdotuI2.o
CIBLAS1_REF = ref/csumI2.o ref/scasumI2.o ref/scnrm2I2.o ref/cdotcI2.o ref/cdotuI2.o

CORE_OBJ_REF = $(SSEO) $(DIBLAS1_REF) $(SIBLAS1_REF) $(ZIBLAS1_REF) $(CIBLAS1_REF) $(IBLAS1)

%o:%c
	PYTHONPATH=../gen python -m cogapp -r -D args=../default_args.json -D params=../params.json $<
	$(C_COMPILER) $(CFLAGS) -c $< -o $@


$(SSEO): %o:%c ../Env/env.h
	PYTHONPATH=../gen python -m cogapp -r -D args=../default_args.json -D params=../params.json $<
	$(C_COMPILER) $(CFLAGS) -c $< -o $@

$(OBJ): config.h
$(MPOBJ): config.h

.PHONY: Indexed MPI_Indexed IndexedBLAS

Indexed:
	cd ../Indexed && $(MAKE) Indexed

MPI_Indexed: Indexed
	cd ../MPI_Indexed && $(MAKE) MPI_Indexed

IndexedBLAS: Indexed $(ROOT_DIR)/$(BUILD)/libs/libiblas.a

$(ROOT_DIR)/$(BUILD)/libs/libiblas.a: $(CORE_OBJ)
	ar -crs $(ROOT_DIR)/$(BUILD)/libs/libiblas.a $(CORE_OBJ)

iblas_ref: Indexed $(ROOT_DIR)/$(BUILD)/libs/libiblas_ref.a Iblas_header rblas_header

$(ROOT_DIR)/$(BUILD)/libs/libiblas_ref.a: $(CORE_OBJ)
	cd ref && $(MAKE) --no-print-directory
	ar -crs $(ROOT_DIR)/$(BUILD)/libs/libiblas_ref.a $(CORE_OBJ_REF)

clean:
	-rm -rf *.o
	cd ref && make clean
	cd $(ROOT_DIR)/$(BUILD)/libs && rm *.a

