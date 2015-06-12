#TARGETS := acc_rcdotc$(EXE) acc_rcdotu$(EXE) acc_rcsum$(EXE) acc_rscasum$(EXE) acc_rscnrm2$(EXE) \
#           acc_rdasum$(EXE) acc_rddot$(EXE) acc_rdnrm2$(EXE) acc_rdsum$(EXE) \
#           acc_rsasum$(EXE) acc_rsdot$(EXE) acc_rsnrm2$(EXE) acc_rssum$(EXE) \
#           acc_rzdotc$(EXE) acc_rzdotu$(EXE) acc_rzsum$(EXE) acc_rdzasum$(EXE) acc_rdznrm2$(EXE) \

TARGETS := acc_rdsum$(EXE) acc_rzsum$(EXE) acc_rssum$(EXE) acc_rcsum$(EXE)

ifneq ($(BLAS),)
TARGETS += #acc_dasum$(EXE) acc_ddot$(EXE) acc_dnrm2$(EXE) \
           #acc_zdotc$(EXE) acc_zdotu$(EXE) acc_dzasum$(EXE) acc_dznrm2$(EXE) \
           #acc_sasum$(EXE) acc_sdot$(EXE) acc_snrm2$(EXE) \
           #acc_cdotc$(EXE) acc_cdotu$(EXE) acc_scasum$(EXE) acc_scnrm2$(EXE) \
           #acc_idamax$(EXE) acc_izamax$(EXE) acc_isamax$(EXE) acc_icamax$(EXE)
endif

SUBDIRS :=

#acc_cdotc$(EXE)_DEPS = $$(LIBTEST) acc_cdotc.o
#acc_cdotu$(EXE)_DEPS = $$(LIBTEST) acc_cdotu.o
#acc_dasum$(EXE)_DEPS = $$(LIBTEST) acc_dasum.o
#acc_ddot$(EXE)_DEPS = $$(LIBTEST) acc_ddot.o
#acc_dnrm2$(EXE)_DEPS = $$(LIBTEST) acc_dnrm2.o
#acc_dzasum$(EXE)_DEPS = $$(LIBTEST) acc_dzasum.o
#acc_dznrm2$(EXE)_DEPS = $$(LIBTEST) acc_dznrm2.o
#acc_rcdotc$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_rcdotc.o
#acc_rcdotu$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_rcdotu.o
acc_rcsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_rcsum.o
#acc_rdasum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_rdasum.o
#acc_rddot$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_rddot.o
#acc_rdnrm2$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_rdnrm2.o
acc_rdsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_rdsum.o
#acc_rdzasum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_rdzasum.o
#acc_rdznrm2$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_rdznrm2.o
#acc_rsasum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_rsasum.o
#acc_rscasum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_rscasum.o
#acc_rscnrm2$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_rscnrm2.o
#acc_rsdot$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_rsdot.o
#acc_rsnrm2$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_rsnrm2.o
acc_rssum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_rssum.o
#acc_rzdotc$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_rzdotc.o
#acc_rzdotu$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_rzdotu.o
acc_rzsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_rzsum.o
#acc_sasum$(EXE)_DEPS = $$(LIBTEST) acc_sasum.o
#acc_scasum$(EXE)_DEPS = $$(LIBTEST) acc_scasum.o
#acc_scnrm2$(EXE)_DEPS = $$(LIBTEST) acc_scnrm2.o
#acc_sdot$(EXE)_DEPS = $$(LIBTEST) acc_sdot.o
#acc_snrm2$(EXE)_DEPS = $$(LIBTEST) acc_snrm2.o
#acc_zdotc$(EXE)_DEPS = $$(LIBTEST) acc_zdotc.o
#acc_zdotu$(EXE)_DEPS = $$(LIBTEST) acc_zdotu.o

#acc_rdgemv$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_rdgemv.o
#acc_prdgemv$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) $$(LIBMPIREPROBLAS) acc_prdgemv.o
#acc_prbdgemv$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) $$(LIBMPIREPROBLAS) acc_prbdgemv.o
#acc_pdgemv$(EXE)_DEPS = $$(LIBTEST) acc_pdgemv.o
#acc_dgemv$(EXE)_DEPS = $$(LIBTEST) acc_dgemv.o
