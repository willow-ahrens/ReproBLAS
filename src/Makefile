ROOT_DIR=..
include $(ROOT_DIR)/Make.local

all: LIB

libs = iblas reproblas_seq
ifneq ($(MPIC_COMPILER),)
libs += reproblas_mpi
endif

LIB : $(libs)
	@cd IndexedFP && $(MAKE) --no-print-directory

GENERATED = dasumI2.c dnrm2I2.c ddotI2.c\
            sasumI2.c ssumI2.c snrm2I2.c sdotI2.c\
            zsumI2.c dzasumI2.c dznrm2I2.c zdotcI2.c zdotuI2.c\
            csumI2.c scasumI2.c scnrm2I2.c cdotcI2.c cdotuI2.c\


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

RBLAS1_OBJ = rdzasum.o rzsum.o rdznrm2.o rzdot.o \
	rcdot.o rcsum.o rscasum.o rscnrm2.o \
	rdsum.o rdasum.o rdnrm2.o rddot.o \
	rssum.o rsasum.o rsnrm2.o rsdot.o

MPOBJ = prdasum.o  prdsum.o prddot.o prdnrm2.o \
	prsasum.o  prssum.o prsdot.o prsnrm2.o \
	prdzasum.o przsum.o przdot.o prdznrm2.o \
	prscasum.o prcsum.o prcdot.o prscnrm2.o

$(GENERATED): %.c: %.py %.tune
	python $<

%o:%c
	python -m gen.cogapp -r -D args=default_args.json -D params=params.json $<
	$(C_COMPILER) $(CFLAGS) -c $< -o $@

$(MPOBJ) : prblas.h rblas1.h 
$(RBLAS1_OBJ) : rblas1.h

$(SSEO): %o:%c Env/env.h
	python -m gen.cogapp -r -D args=default_args.json -D params=params.json $<
	$(C_COMPILER) $(CFLAGS) -c $< -o $@

$(OBJ): config.h
$(MPOBJ): config.h

Compile : $(OBJ) $(MPOBJ)

$(MPOBJ): %.o : %.c
	$(MPIC_COMPILER) $(CFLAGS) -c $< -o $@

.PHONY: IndexedFP MPIndexedFP

IndexedFP:
	cd IndexedFP && $(MAKE) IndexedFP

MPIndexedFP: IndexedFP
	cd IndexedFP && $(MAKE) MPIndexedFP

iblas: IndexedFP $(ROOT_DIR)/$(BUILD)/libs/libiblas.a Iblas_header rblas_header
$(ROOT_DIR)/$(BUILD)/libs/libiblas.a: $(CORE_OBJ)
	ar -crs $(ROOT_DIR)/$(BUILD)/libs/libiblas.a $(CORE_OBJ)

iblas_ref: IndexedFP $(ROOT_DIR)/$(BUILD)/libs/libiblas_ref.a Iblas_header rblas_header
$(ROOT_DIR)/$(BUILD)/libs/libiblas_ref.a: $(CORE_OBJ)
	cd ref && $(MAKE) --no-print-directory
	ar -crs $(ROOT_DIR)/$(BUILD)/libs/libiblas_ref.a $(CORE_OBJ_REF)

rblas_header : rblas1.mod
	@echo '#ifndef _REPRODUCIBLE_BLAS__H_' > rblas_tmp.h
	@echo '#define _REPRODUCIBLE_BLAS__H_' >> rblas_tmp.h
	@echo '#include "IndexedFP.h"' >> rblas_tmp.h
	@echo '#include "Iblas.h"' >> rblas_tmp.h
	@cat rblas1.mod >> rblas_tmp.h
	@echo '#endif' >> rblas_tmp.h
	@mv rblas_tmp.h $(ROOT_DIR)/$(BUILD)/include/rblas.h

rblas_mpi.header : prblas1.h
	@echo '#ifndef _REPRODUCIBLE_PAR_BLAS__H_' > reproblas_mpi_tmp.h
	@echo '#define _REPRODUCIBLE_PAR_BLAS__H_' >> reproblas_mpi_tmp.h
	@echo '#include "rblas.h"' >> reproblas_mpi_tmp.h
	@echo '#include "MPIndexedFP.h"' >> reproblas_mpi_tmp.h
	@cat prblas1.h >> reproblas_mpi_tmp.h
	@echo '#endif' >> reproblas_mpi_tmp.h
	@mv reproblas_mpi_tmp.h $(ROOT_DIR)/$(BUILD)/include/rblas_mpi.h

Iblas_header : Iblas1.mod
	@echo '#ifndef _INDEXED_BLAS__H_' > Iblas_tmp.h
	@echo '#define _INDEXED_BLAS__H_' >> Iblas_tmp.h
	@echo '#include "IndexedFP.h"' >> Iblas_tmp.h
	@cat Iblas1.mod >> Iblas_tmp.h
	@cat accumulator.mod >> Iblas_tmp.h
	@echo '#endif' >> Iblas_tmp.h
	@mv Iblas_tmp.h $(ROOT_DIR)/$(BUILD)/include/Iblas.h

reproblas_seq: iblas $(ROOT_DIR)/$(BUILD)/libs/libreproblas_seq.a rblas_header
$(ROOT_DIR)/$(BUILD)/libs/libreproblas_seq.a: $(RBLAS1_OBJ)
	ar -crs $(ROOT_DIR)/$(BUILD)/libs/libreproblas_seq.a $(RBLAS1_OBJ)

reproblas_mpi: MPIndexedFP reproblas_seq $(ROOT_DIR)/$(BUILD)/libs/libreproblas_mpi.a rblas_mpi.header
$(ROOT_DIR)/$(BUILD)/libs/libreproblas_mpi.a: $(MPOBJ)
	-ar -crs $(ROOT_DIR)/$(BUILD)/libs/libreproblas_mpi.a $(MPOBJ)

clean:
	-rm -rf *.o
	-rm -f $(GENERATED)
	cd IndexedFP && make clean
	cd ref && make clean
	cd $(ROOT_DIR)/$(BUILD)/include && rm *.h
	cd $(ROOT_DIR)/$(BUILD)/libs && rm *.a

