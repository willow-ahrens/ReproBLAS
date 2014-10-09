ROOT_DIR=../..
include $(ROOT_DIR)/Make.local

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

all: $(SIBLAS1) $(DIBLAS1) $(CIBLAS1) $(ZIBLAS1) 

TMPL_RDBLAS1 = dIBlas1_AVX.c
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

clean:
	-rm -rf *.o
