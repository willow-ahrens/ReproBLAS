TARGETS := libreproblas.a
SUBDIRS := 

INSTALL_LIB := $(TARGETS)

RBLAS1 = rdsum.o rdasum.o rdnrm2.o rddot.o   \
         rdzasum.o rzsum.o rdznrm2.o rzdot.o \
         rssum.o rsasum.o rsnrm2.o rsdot.o   \
         rcdot.o rcsum.o rscasum.o rscnrm2.o

libreproblas.a_DEPS = $(RBLAS1) $$(LIBINDEXED) $$(LIBINDEXEDBLAS)
