TARGETS := libreproblas.a
SUBDIRS := 

INSTALL_LIB := $(TARGETS)

LIBREPROBLAS := $(OBJPATH)/libreproblas.a

libreproblas.a_DEPS = $$(LIBINDEXED) $$(LIBINDEXEDBLAS)                        \
                      rdsum.o rdasum.o rdnrm2.o rddot.o rdgemv.o rdgemm.o              \
                      rzsum_sub.o rdzasum.o rdznrm2.o rzdotc_sub.o rzdotu_sub.o rzgemv.o \
                      rssum.o rsasum.o rsnrm2.o rsdot.o                        \
                      rcsum_sub.o rscasum.o rscnrm2.o rcdotc_sub.o rcdotu_sub.o\
