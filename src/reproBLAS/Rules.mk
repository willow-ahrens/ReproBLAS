TARGETS := libreproblas.a
SUBDIRS :=

INSTALL_LIB := $(TARGETS)

LIBREPROBLAS := $(OBJPATH)/libreproblas.a

libreproblas.a_DEPS = $$(LIBIDXD) $$(LIBIDXDBLAS)                          \
                      rdsum.o rdasum.o rdnrm2.o rddot.o                    \
                      rzsum_sub.o rdzasum.o rdznrm2.o rzdotc_sub.o         \
                        rzdotu_sub.o                                       \
                      rssum.o rsasum.o rsnrm2.o rsdot.o                    \
                      rcsum_sub.o rscasum.o rscnrm2.o rcdotc_sub.o         \
                        rcdotu_sub.o                                       \
                      rdgemv.o rdgemm.o                                    \
                      rzgemv.o rzgemm.o                                    \
                      rsgemv.o rsgemm.o                                    \
                      rcgemv.o rcgemm.o                                    \
                      dsum.o dasum.o dnrm2.o ddot.o                        \
                      zsum_sub.o dzasum.o dznrm2.o zdotc_sub.o zdotu_sub.o \
                      ssum.o sasum.o snrm2.o sdot.o                        \
                      csum_sub.o scasum.o scnrm2.o cdotc_sub.o cdotu_sub.o \
                      dgemv.o dgemm.o                                      \
                      zgemv.o zgemm.o                                      \
                      sgemv.o sgemm.o                                      \
                      cgemv.o cgemm.o                                      \
