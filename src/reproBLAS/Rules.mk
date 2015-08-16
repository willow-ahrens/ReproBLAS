TARGETS := libreproblas.a
SUBDIRS :=

INSTALL_LIB := $(TARGETS)

LIBREPROBLAS := $(OBJPATH)/libreproblas.a

libreproblas.a_DEPS = $$(LIBIDXD) $$(LIBIDXDBLAS)                        \
                      rdsum.o rdasum.o rdnrm2.o rddot.o rdgemv.o rdgemm.o              \
                      rzsum_sub.o rdzasum.o rdznrm2.o rzdotc_sub.o rzdotu_sub.o rzgemv.o rzgemm.o \
                      rssum.o rsasum.o rsnrm2.o rsdot.o                        \
                      rcsum_sub.o rscasum.o rscnrm2.o rcdotc_sub.o rcdotu_sub.o\
                      dsum.o dasum.o dnrm2.o ddot.o dgemv.o dgemm.o              \
                      zsum_sub.o dzasum.o dznrm2.o zdotc_sub.o zdotu_sub.o zgemv.o zgemm.o \
                      ssum.o sasum.o snrm2.o sdot.o                        \
                      csum_sub.o scasum.o scnrm2.o cdotc_sub.o cdotu_sub.o\
