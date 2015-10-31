TARGETS := bench_damax$(EXE) bench_damaxm$(EXE)                     \
           bench_zamax$(EXE) bench_zamaxm$(EXE)                     \
           bench_samax$(EXE) bench_samaxm$(EXE)                     \
           bench_camax$(EXE) bench_camaxm$(EXE)                     \
           bench_dsum$(EXE)                                         \
           bench_zsum$(EXE)                                         \
           bench_ssum$(EXE)                                         \
           bench_csum$(EXE)                                         \
           bench_rdasum$(EXE)  bench_rdnrm2$(EXE) bench_rdsum$(EXE) \
             bench_rddot$(EXE)                                      \
           bench_rzsum$(EXE) bench_rdzasum$(EXE) bench_rdznrm2$(EXE)\
             bench_rzdotc$(EXE) bench_rzdotu$(EXE)                  \
           bench_rsasum$(EXE)  bench_rsnrm2$(EXE) bench_rssum$(EXE) \
             bench_rsdot$(EXE)                                      \
           bench_rcsum$(EXE) bench_rscasum$(EXE) bench_rscnrm2$(EXE)\
             bench_rcdotc$(EXE) bench_rcdotu$(EXE)                  \
           bench_rdgemv$(EXE) bench_rdgemm$(EXE)                    \
           bench_rzgemv$(EXE) bench_rzgemm$(EXE)                    \
           bench_rsgemv$(EXE) bench_rsgemm$(EXE)                    \
           bench_rcgemv$(EXE) bench_rcgemm$(EXE)                    \
           bench_ddiconv$(EXE)                                      \
           bench_zziconv$(EXE)                                      \
           bench_ssiconv$(EXE)                                      \
           bench_cciconv$(EXE)                                      \
           bench_didiadd$(EXE)                                      \
           bench_ziziadd$(EXE)                                      \
           bench_sisiadd$(EXE)                                      \
           bench_ciciadd$(EXE)

ifneq ($(BLAS),)

TARGETS += bench_idamax$(EXE)                                  \
           bench_izamax$(EXE)                                  \
           bench_isamax$(EXE)                                  \
           bench_icamax$(EXE)                                  \
           bench_dasum$(EXE) bench_dnrm2$(EXE) bench_ddot$(EXE)\
           bench_dzasum$(EXE) bench_dznrm2$(EXE)               \
             bench_zdotc$(EXE) bench_zdotu$(EXE)               \
           bench_sasum$(EXE) bench_snrm2$(EXE) bench_sdot$(EXE)\
           bench_scasum$(EXE) bench_scnrm2$(EXE)               \
             bench_cdotc$(EXE) bench_cdotu$(EXE)               \
           bench_dgemv$(EXE) bench_dgemm$(EXE)                 \
           bench_zgemv$(EXE) bench_zgemm$(EXE)                 \
           bench_sgemv$(EXE) bench_sgemm$(EXE)                 \
           bench_cgemv$(EXE) bench_cgemm$(EXE)
endif

SUBDIRS :=

bench_camax$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXDBLAS) bench_camax.o
bench_camaxm$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXDBLAS) bench_camaxm.o
bench_cciconv$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXDBLAS) bench_cciconv.o
bench_cdotc$(EXE)_DEPS = $$(LIBTEST) bench_cdotc.o
bench_cdotu$(EXE)_DEPS = $$(LIBTEST) bench_cdotu.o
bench_cgemm$(EXE)_DEPS = $$(LIBTEST) bench_cgemm.o
bench_cgemv$(EXE)_DEPS = $$(LIBTEST) bench_cgemv.o
bench_ciciadd$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXDBLAS) bench_ciciadd.o
bench_csum$(EXE)_DEPS = $$(LIBTEST) bench_csum.o
bench_damax$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXDBLAS) bench_damax.o
bench_damaxm$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXDBLAS) bench_damaxm.o
bench_dasum$(EXE)_DEPS = $$(LIBTEST) bench_dasum.o
bench_ddiconv$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXDBLAS) bench_ddiconv.o
bench_ddot$(EXE)_DEPS = $$(LIBTEST) bench_ddot.o
bench_dgemm$(EXE)_DEPS = $$(LIBTEST) bench_dgemm.o
bench_dgemv$(EXE)_DEPS = $$(LIBTEST) bench_dgemv.o
bench_didiadd$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXDBLAS) bench_didiadd.o
bench_dnrm2$(EXE)_DEPS = $$(LIBTEST) bench_dnrm2.o
bench_dsum$(EXE)_DEPS = $$(LIBTEST) bench_dsum.o
bench_dzasum$(EXE)_DEPS = $$(LIBTEST) bench_dzasum.o
bench_dznrm2$(EXE)_DEPS = $$(LIBTEST) bench_dznrm2.o
bench_icamax$(EXE)_DEPS = $$(LIBTEST) bench_icamax.o
bench_idamax$(EXE)_DEPS = $$(LIBTEST) bench_idamax.o
bench_isamax$(EXE)_DEPS = $$(LIBTEST) bench_isamax.o
bench_izamax$(EXE)_DEPS = $$(LIBTEST) bench_izamax.o
bench_rcdotc$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rcdotc.o
bench_rcdotu$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rcdotu.o
bench_rcgemm$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rcgemm.o
bench_rcgemv$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rcgemv.o
bench_rcsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rcsum.o
bench_rdasum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rdasum.o
bench_rddot$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rddot.o
bench_rdgemm$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rdgemm.o
bench_rdgemv$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rdgemv.o
bench_rdnrm2$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rdnrm2.o
bench_rdsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rdsum.o
bench_rdzasum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rdzasum.o
bench_rdznrm2$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rdznrm2.o
bench_rsasum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rsasum.o
bench_rscasum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rscasum.o
bench_rscnrm2$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rscnrm2.o
bench_rsdot$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rsdot.o
bench_rsgemm$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rsgemm.o
bench_rsgemv$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rsgemv.o
bench_rsnrm2$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rsnrm2.o
bench_rssum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rssum.o
bench_rzdotc$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rzdotc.o
bench_rzdotu$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rzdotu.o
bench_rzgemm$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rzgemm.o
bench_rzgemv$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rzgemv.o
bench_rzsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) bench_rzsum.o
bench_samax$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXDBLAS) bench_samax.o
bench_samaxm$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXDBLAS) bench_samaxm.o
bench_sasum$(EXE)_DEPS = $$(LIBTEST) bench_sasum.o
bench_scasum$(EXE)_DEPS = $$(LIBTEST) bench_scasum.o
bench_scnrm2$(EXE)_DEPS = $$(LIBTEST) bench_scnrm2.o
bench_sdot$(EXE)_DEPS = $$(LIBTEST) bench_sdot.o
bench_sgemm$(EXE)_DEPS = $$(LIBTEST) bench_sgemm.o
bench_sgemv$(EXE)_DEPS = $$(LIBTEST) bench_sgemv.o
bench_sisiadd$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXDBLAS) bench_sisiadd.o
bench_snrm2$(EXE)_DEPS = $$(LIBTEST) bench_snrm2.o
bench_ssiconv$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXDBLAS) bench_ssiconv.o
bench_ssum$(EXE)_DEPS = $$(LIBTEST) bench_ssum.o
bench_zamax$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXDBLAS) bench_zamax.o
bench_zamaxm$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXDBLAS) bench_zamaxm.o
bench_zdotc$(EXE)_DEPS = $$(LIBTEST) bench_zdotc.o
bench_zdotu$(EXE)_DEPS = $$(LIBTEST) bench_zdotu.o
bench_zgemm$(EXE)_DEPS = $$(LIBTEST) bench_zgemm.o
bench_zgemv$(EXE)_DEPS = $$(LIBTEST) bench_zgemv.o
bench_ziziadd$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXDBLAS) bench_ziziadd.o
bench_zsum$(EXE)_DEPS = $$(LIBTEST) bench_zsum.o
bench_zziconv$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXDBLAS) bench_zziconv.o

bench_camax$(EXE)_LIBS = -lm
bench_camaxm$(EXE)_LIBS = -lm
bench_cciconv$(EXE)_LIBS = -lm
bench_cdotc$(EXE)_LIBS = -lm
bench_cdotu$(EXE)_LIBS = -lm
bench_cgemm$(EXE)_LIBS = -lm
bench_cgemv$(EXE)_LIBS = -lm
bench_ciciadd$(EXE)_LIBS = -lm
bench_csum$(EXE)_LIBS = -lm
bench_damax$(EXE)_LIBS = -lm
bench_damaxm$(EXE)_LIBS = -lm
bench_dasum$(EXE)_LIBS = -lm
bench_ddiconv$(EXE)_LIBS = -lm
bench_ddot$(EXE)_LIBS = -lm
bench_dgemm$(EXE)_LIBS = -lm
bench_dgemv$(EXE)_LIBS = -lm
bench_didiadd$(EXE)_LIBS = -lm
bench_dnrm2$(EXE)_LIBS = -lm
bench_dsum$(EXE)_LIBS = -lm
bench_dzasum$(EXE)_LIBS = -lm
bench_dznrm2$(EXE)_LIBS = -lm
bench_icamax$(EXE)_LIBS = -lm
bench_idamax$(EXE)_LIBS = -lm
bench_isamax$(EXE)_LIBS = -lm
bench_izamax$(EXE)_LIBS = -lm
bench_rcdotc$(EXE)_LIBS = -lm
bench_rcdotu$(EXE)_LIBS = -lm
bench_rcgemm$(EXE)_LIBS = -lm
bench_rcgemv$(EXE)_LIBS = -lm
bench_rcsum$(EXE)_LIBS = -lm
bench_rdasum$(EXE)_LIBS = -lm
bench_rddot$(EXE)_LIBS = -lm
bench_rdgemm$(EXE)_LIBS = -lm
bench_rdgemv$(EXE)_LIBS = -lm
bench_rdnrm2$(EXE)_LIBS = -lm
bench_rdsum$(EXE)_LIBS = -lm
bench_rdzasum$(EXE)_LIBS = -lm
bench_rdznrm2$(EXE)_LIBS = -lm
bench_rsasum$(EXE)_LIBS = -lm
bench_rscasum$(EXE)_LIBS = -lm
bench_rscnrm2$(EXE)_LIBS = -lm
bench_rsdot$(EXE)_LIBS = -lm
bench_rsgemm$(EXE)_LIBS = -lm
bench_rsgemv$(EXE)_LIBS = -lm
bench_rsnrm2$(EXE)_LIBS = -lm
bench_rssum$(EXE)_LIBS = -lm
bench_rzdotc$(EXE)_LIBS = -lm
bench_rzdotu$(EXE)_LIBS = -lm
bench_rzgemm$(EXE)_LIBS = -lm
bench_rzgemv$(EXE)_LIBS = -lm
bench_rzsum$(EXE)_LIBS = -lm
bench_samax$(EXE)_LIBS = -lm
bench_samaxm$(EXE)_LIBS = -lm
bench_sasum$(EXE)_LIBS = -lm
bench_scasum$(EXE)_LIBS = -lm
bench_scnrm2$(EXE)_LIBS = -lm
bench_sdot$(EXE)_LIBS = -lm
bench_sgemm$(EXE)_LIBS = -lm
bench_sgemv$(EXE)_LIBS = -lm
bench_sisiadd$(EXE)_LIBS = -lm
bench_snrm2$(EXE)_LIBS = -lm
bench_ssiconv$(EXE)_LIBS = -lm
bench_ssum$(EXE)_LIBS = -lm
bench_zamax$(EXE)_LIBS = -lm
bench_zamaxm$(EXE)_LIBS = -lm
bench_zdotc$(EXE)_LIBS = -lm
bench_zdotu$(EXE)_LIBS = -lm
bench_zgemm$(EXE)_LIBS = -lm
bench_zgemv$(EXE)_LIBS = -lm
bench_ziziadd$(EXE)_LIBS = -lm
bench_zsum$(EXE)_LIBS = -lm
bench_zziconv$(EXE)_LIBS = -lm
