TARGETS := libreproblas.a
SUBDIRS := 

INSTALL_LIB := $(TARGETS)

LIBREPROBLAS := $(OBJPATH)/libreproblas.a

libreproblas.a_DEPS = $$(LIBINDEXED) $$(LIBINDEXEDBLAS)   \
                      rdsum.o rdasum.o rdnrm2.o rddot.o   \
                      rdzasum.o rzsum.o rdznrm2.o rzdot.o \
                      rssum.o rsasum.o rsnrm2.o rsdot.o   \
                      rcdot.o rcsum.o rscasum.o rscnrm2.o \
                      rdgemv.o
