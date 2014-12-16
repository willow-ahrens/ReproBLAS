ROOT_DIR=..
include $(ROOT_DIR)/Make.local

all: LIB

libs = iblas reproblas_seq
ifneq ($(MPIC_COMPILER),)
libs += reproblas_mpi
endif

LIB : $(libs)
	@cd IndexedFP && $(MAKE) --no-print-directory

DIBLAS1 = DASUMI2.o DSUMI2.o DNRM2I2.o DDOTI2.o
SIBLAS1 = SASUMI2.o SSUMI2.o SNRM2I2.o SDOTI2.o
ZIBLAS1 = zsumI2.o dzasumI2.o dznrm2I2.o zdotcI2.o zdotuI2.o
CIBLAS1 = csumI2.o scasumI2.o scnrm2I2.o cdotcI2.o cdotuI2.o

IBLAS1  = ssumI.o sasumI.o sdotI.o snrm2I.o \
	dsumI.o dasumI.o ddotI.o dnrm2I.o	\
	zdotI.o zsumI.o dzasumI.o dznrm2I.o \
	csumI.o scasumI.o cdotI.o scnrm2I.o

SSEO = samax.o samaxm.o damax.o damaxm.o zamax.o zamaxm.o camax.o camaxm.o	\
	dIAccum.o sIAccum.o zIAccum.o cIAccum.o

CORE_OBJ = $(SSEO) $(DIBLAS1) $(SIBLAS1) $(ZIBLAS1) $(CIBLAS1) $(IBLAS1)

RBLAS1_OBJ = rdzasum.o rzsum.o rdznrm2.o rzdot.o \
	rcdot.o rcsum.o rscasum.o rscnrm2.o \
	rdsum.o rdasum.o rdnrm2.o rddot.o \
	rssum.o rsasum.o rsnrm2.o rsdot.o

MPOBJ = prdasum.o  prdsum.o prddot.o prdnrm2.o \
	prsasum.o  prssum.o prsdot.o prsnrm2.o \
	prdzasum.o przsum.o przdot.o prdznrm2.o \
	prscasum.o prcsum.o prcdot.o prscnrm2.o

$(MPOBJ) : prblas.h rblas1.h 
$(RBLAS1_OBJ) : rblas1.h

$(SSEO): %o:%c Env/env.h
	$(C_COMPILER) $(CFLAGS) -c $< -o $@

$(OBJ): config.h
$(MPOBJ): config.h

Compile : $(OBJ) $(MPOBJ)

%o:%c
	$(C_COMPILER) $(CFLAGS) -c $< -o $@

$(MPOBJ): %.o : %.c
	$(MPIC_COMPILER) $(CFLAGS) -c $< -o $@

TMPL_RDBLAS1 = dIBlas1.c
TMPL_RSBLAS1 = sIBlas1.c
TMPL_ZSUMI2 = zsumI2.c
TMPL_ZDOTI2 = zdotI2.c
TMPL_CSUMI2 = csumI2.c
TMPL_CDOTI2 = cdotI2.c

ifeq ($(vect),SISD)
TMPL_RDBLAS1 = dIBlas1_SISD.c
TMPL_RSBLAS1 = sIBlas1_SISD.c
TMPL_ZSUMI2  = zsumI2_SISD.c
TMPL_ZDOTI2  = zdotI2_SISD.c
TMPL_CSUMI2  = csumI2_SISD.c
TMPL_CDOTI2  = cdotI2_SISD.c
else
ifeq ($(C_COMPILERNAME),pgcc)
TMPL_RDBLAS1 = dIBlas1_SISD.c
TMPL_RSBLAS1 = sIBlas1_SISD.c
TMPL_ZSUMI2  = zsumI2_SISD.c
TMPL_ZDOTI2  = zdotI2_SISD.c
TMPL_CSUMI2 = csumI2_SISD.c
TMPL_CDOTI2  = cdotI2_SISD.c
endif
endif

${DIBLAS1} : $(TMPL_RDBLAS1)
	@echo $(patsubst %.o,'#define %',$@) > dIBlas1.h
	$(C_COMPILER) $(CFLAGS) $(FPFLAGS) -c $(TMPL_RDBLAS1) -o $@
	@rm dIBlas1.h

${SIBLAS1} : $(TMPL_RSBLAS1)
	@echo $(patsubst %.o,'#define %',$@) > sIBlas1.h
	$(C_COMPILER) $(CFLAGS) $(FPFLAGS) -c $(TMPL_RSBLAS1) -o $@
	@rm sIBlas1.h

zsumI2.o : $(TMPL_ZSUMI2)
	@echo '#define ZSUMI2' > zIBlas1.h
	$(C_COMPILER) $(CFLAGS) $(FPFLAGS) -c $< -o $@
	@rm zIBlas1.h
dzasumI2.o : $(TMPL_ZSUMI2)
	@echo '#define DZASUMI2' > zIBlas1.h
	$(C_COMPILER) $(CFLAGS) $(FPFLAGS) -c $< -o $@
	@rm zIBlas1.h
dznrm2I2.o : $(TMPL_ZSUMI2)
	@echo '#define DZNRM2I2' > zIBlas1.h
	$(C_COMPILER) $(CFLAGS) $(FPFLAGS) -c $< -o $@
	@rm zIBlas1.h
zdotcI2.o : $(TMPL_ZDOTI2)
	@echo '#define ZDOTCI2' > zIBlas1.h
	$(C_COMPILER) $(CFLAGS) $(FPFLAGS) -c $< -o $@
	@rm zIBlas1.h
zdotuI2.o : $(TMPL_ZDOTI2)
	@echo '#define ZDOTUI2' > zIBlas1.h
	$(C_COMPILER) $(CFLAGS) $(FPFLAGS) -c $< -o $@
	@rm zIBlas1.h


csumI2.o : $(TMPL_CSUMI2)
	@echo '#define CSUMI2' > cIBlas1.h
	$(C_COMPILER) $(CFLAGS) $(FPFLAGS) -c $< -o $@
	@rm cIBlas1.h
scasumI2.o : $(TMPL_CSUMI2)
	@echo '#define SCASUMI2' > cIBlas1.h
	$(C_COMPILER) $(CFLAGS) $(FPFLAGS) -c $< -o $@
	@rm cIBlas1.h
scnrm2I2.o : $(TMPL_CSUMI2)
	@echo '#define SCNRM2I2' > cIBlas1.h
	$(C_COMPILER) $(CFLAGS) $(FPFLAGS) -c $< -o $@
	@rm cIBlas1.h
cdotcI2.o : $(TMPL_CDOTI2)
	@echo '#define CDOTCI2' > cIBlas1.h
	$(C_COMPILER) $(CFLAGS) $(FPFLAGS) -c $< -o $@
	@rm cIBlas1.h
cdotuI2.o : $(TMPL_CDOTI2)
	@echo '#define CDOTUI2' > cIBlas1.h
	$(C_COMPILER) $(CFLAGS) $(FPFLAGS) -c $< -o $@
	@rm cIBlas1.h

.PHONY: IndexedFP MPIndexedFP

IndexedFP:
	cd IndexedFP && $(MAKE) IndexedFP

MPIndexedFP: IndexedFP
	cd IndexedFP && $(MAKE) MPIndexedFP

iblas: IndexedFP $(ROOT_DIR)/$(BUILD)/libs/libiblas.a Iblas_header
$(ROOT_DIR)/$(BUILD)/libs/libiblas.a: $(CORE_OBJ)
	ar -crs $(ROOT_DIR)/$(BUILD)/libs/libiblas.a $(CORE_OBJ)

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
	@cp $(ROOT_DIR)/$(BUILD)/include/rblas_mpi.h $(ROOT_DIR)/$(BUILD)/include/prblas.h

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
	cd IndexedFP && make clean

