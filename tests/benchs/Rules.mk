TARGETS := bench_damax$(EXE) bench_damaxm$(EXE) \
           bench_zamax$(EXE) bench_zamaxm$(EXE) \
           bench_samax$(EXE) bench_samaxm$(EXE) \
           bench_camax$(EXE) bench_camaxm$(EXE) \
           bench_rcdotc$(EXE) bench_rcdotu$(EXE) bench_rcsum$(EXE) bench_rscasum$(EXE) bench_rscnrm2$(EXE) \
           bench_rdasum$(EXE) bench_rddot$(EXE) bench_rdnrm2$(EXE) bench_rdsum$(EXE) \
           bench_rsasum$(EXE) bench_rsdot$(EXE) bench_rsnrm2$(EXE) bench_rssum$(EXE) \
           bench_rzdotc$(EXE) bench_rzdotu$(EXE) bench_rzsum$(EXE) bench_rdzasum$(EXE) bench_rdznrm2$(EXE) \
           bench_rdgemv$(EXE) bench_prdgemv$(EXE)

ifneq ($(BLAS),)
TARGETS += bench_dasum$(EXE) bench_ddot$(EXE) bench_dnrm2$(EXE) \
           bench_zdotc$(EXE) bench_zdotu$(EXE) bench_dzasum$(EXE) bench_dznrm2$(EXE) \
           bench_sasum$(EXE) bench_sdot$(EXE) bench_snrm2$(EXE) \
           bench_cdotc$(EXE) bench_cdotu$(EXE) bench_scasum$(EXE) bench_scnrm2$(EXE) \
           bench_idamax$(EXE) bench_izamax$(EXE) bench_isamax$(EXE) bench_icamax$(EXE)\
	   bench_dgemv$(EXE) bench_pdgemv$(EXE)
endif

SUBDIRS :=

bench_camax$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXEDBLAS) bench_camax.o
bench_camaxm$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXEDBLAS) bench_camaxm.o
bench_cdotc$(EXE)_DEPS = $$(LIBTEST) bench_cdotc.o
bench_cdotu$(EXE)_DEPS = $$(LIBTEST) bench_cdotu.o
bench_damax$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXEDBLAS) bench_damax.o
bench_damaxm$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXEDBLAS) bench_damaxm.o
bench_dasum$(EXE)_DEPS = $$(LIBTEST) bench_dasum.o
bench_ddot$(EXE)_DEPS = $$(LIBTEST) bench_ddot.o
bench_dnrm2$(EXE)_DEPS = $$(LIBTEST) bench_dnrm2.o
bench_dzasum$(EXE)_DEPS = $$(LIBTEST) bench_dzasum.o
bench_dznrm2$(EXE)_DEPS = $$(LIBTEST) bench_dznrm2.o
bench_icamax$(EXE)_DEPS = $$(LIBTEST) bench_icamax.o
bench_idamax$(EXE)_DEPS = $$(LIBTEST) bench_idamax.o
bench_isamax$(EXE)_DEPS = $$(LIBTEST) bench_isamax.o
bench_izamax$(EXE)_DEPS = $$(LIBTEST) bench_izamax.o
bench_rcdotc$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rcdotc.o
bench_rcdotu$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rcdotu.o
bench_rcsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rcsum.o
bench_rdasum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rdasum.o
bench_rddot$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rddot.o
bench_rdnrm2$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rdnrm2.o
bench_rdsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rdsum.o
bench_rdzasum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rdzasum.o
bench_rdznrm2$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rdznrm2.o
bench_rsasum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rsasum.o
bench_rscasum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rscasum.o
bench_rscnrm2$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rscnrm2.o
bench_rsdot$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rsdot.o
bench_rsnrm2$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rsnrm2.o
bench_rssum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rssum.o
bench_rzdotc$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rzdotc.o
bench_rzdotu$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rzdotu.o
bench_rzsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rzsum.o
bench_samax$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXEDBLAS) bench_samax.o
bench_samaxm$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXEDBLAS) bench_samaxm.o
bench_sasum$(EXE)_DEPS = $$(LIBTEST) bench_sasum.o
bench_scasum$(EXE)_DEPS = $$(LIBTEST) bench_scasum.o
bench_scnrm2$(EXE)_DEPS = $$(LIBTEST) bench_scnrm2.o
bench_sdot$(EXE)_DEPS = $$(LIBTEST) bench_sdot.o
bench_snrm2$(EXE)_DEPS = $$(LIBTEST) bench_snrm2.o
bench_zamax$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXEDBLAS) bench_zamax.o
bench_zamaxm$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXEDBLAS) bench_zamaxm.o
bench_zdotc$(EXE)_DEPS = $$(LIBTEST) bench_zdotc.o
bench_zdotu$(EXE)_DEPS = $$(LIBTEST) bench_zdotu.o

bench_rdgemv$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rdgemv.o
bench_prdgemv$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) $$(LIBMPIREPROBLAS) bench_prdgemv.o
bench_pdgemv$(EXE)_DEPS = $$(LIBTEST) bench_pdgemv.o
bench_dgemv$(EXE)_DEPS = $$(LIBTEST) bench_dgemv.o
