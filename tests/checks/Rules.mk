TARGETS := validate_internal_damax$(EXE) validate_internal_zamax$(EXE) validate_internal_samax$(EXE) validate_internal_camax$(EXE) \
           validate_internal_ufp$(EXE) validate_internal_ufpf$(EXE) \
           validate_internal_dscale$(EXE) validate_internal_sscale$(EXE) \
           validate_internal_dindex$(EXE) validate_internal_sindex$(EXE) \
           validate_internal_dmindex$(EXE) validate_internal_smindex$(EXE) \
           verify_daugsum$(EXE) verify_zaugsum$(EXE) verify_saugsum$(EXE) verify_caugsum$(EXE) \
           validate_internal_daugsum$(EXE) validate_internal_zaugsum$(EXE) validate_internal_saugsum$(EXE) validate_internal_caugsum$(EXE)\
           validate_xblas_ddot$(EXE) validate_xblas_zdot$(EXE)\
           verify_didssq$(EXE) verify_dizssq$(EXE) verify_sisssq$(EXE) verify_sicssq$(EXE) \
           corroborate_rdgemv$(EXE) \
           corroborate_rdgemm$(EXE) \
           corroborate_rzgemv$(EXE) \
           corroborate_rzgemm$(EXE) \
           corroborate_rsgemv$(EXE) \
           corroborate_rsgemm$(EXE) \
           corroborate_rcgemv$(EXE) \
           corroborate_rcgemm$(EXE) \

SUBDIRS :=

LDFLAGS += $(MPILDFLAGS)
CFLAGS += $(MPICFLAGS)

validate_internal_camax$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXDBLAS) validate_internal_camax.o
validate_internal_damax$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXDBLAS) validate_internal_damax.o
validate_internal_samax$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXDBLAS) validate_internal_samax.o
validate_internal_zamax$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXDBLAS) validate_internal_zamax.o
validate_internal_ufp$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXDBLAS) validate_internal_ufp.o
validate_internal_ufpf$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXDBLAS) validate_internal_ufpf.o
validate_internal_dscale$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXD) validate_internal_dscale.o
validate_internal_sscale$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXD) validate_internal_sscale.o
validate_internal_dindex$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXD) validate_internal_dindex.o
validate_internal_sindex$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXD) validate_internal_sindex.o
validate_internal_dmindex$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXD) validate_internal_dmindex.o
validate_internal_smindex$(EXE)_DEPS = $$(LIBTEST) $$(LIBIDXD) validate_internal_smindex.o
verify_daugsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_daugsum.o
verify_zaugsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_zaugsum.o
verify_saugsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_saugsum.o
verify_caugsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_caugsum.o
validate_internal_daugsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_internal_daugsum.o
validate_internal_zaugsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_internal_zaugsum.o
validate_internal_saugsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_internal_saugsum.o
validate_internal_caugsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_internal_caugsum.o
validate_xblas_ddot$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_xblas_ddot.o
validate_xblas_zdot$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_xblas_zdot.o
validate_xblas_sdot$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_xblas_sdot.o
validate_xblas_cdot$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_xblas_cdot.o
verify_didssq$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_didssq.o
verify_dizssq$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_dizssq.o
verify_sisssq$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_sisssq.o
verify_sicssq$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_sicssq.o
corroborate_rdgemv$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) corroborate_rdgemv.o
corroborate_rdgemm$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) corroborate_rdgemm.o
corroborate_rzgemv$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) corroborate_rzgemv.o
corroborate_rzgemm$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) corroborate_rzgemm.o
corroborate_rsgemv$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) corroborate_rsgemv.o
corroborate_rsgemm$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) corroborate_rsgemm.o
corroborate_rcgemv$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) corroborate_rcgemv.o
corroborate_rcgemm$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) corroborate_rcgemm.o

validate_internal_damax$(EXE)_LIBS = -lm
validate_internal_zamax$(EXE)_LIBS = -lm
validate_internal_samax$(EXE)_LIBS = -lm
validate_internal_camax$(EXE)_LIBS = -lm
validate_internal_ufp$(EXE)_LIBS = -lm
validate_internal_ufpf$(EXE)_LIBS = -lm
validate_internal_dscale$(EXE)_LIBS = -lm
validate_internal_sscale$(EXE)_LIBS = -lm
validate_internal_dindex$(EXE)_LIBS = -lm
validate_internal_sindex$(EXE)_LIBS = -lm
validate_internal_dmindex$(EXE)_LIBS = -lm
validate_internal_smindex$(EXE)_LIBS = -lm
verify_daugsum$(EXE)_LIBS = -lm
verify_zaugsum$(EXE)_LIBS = -lm
verify_saugsum$(EXE)_LIBS = -lm
verify_caugsum$(EXE)_LIBS = -lm
validate_internal_daugsum$(EXE)_LIBS = -lm
validate_internal_zaugsum$(EXE)_LIBS = -lm
validate_internal_saugsum$(EXE)_LIBS = -lm
validate_internal_caugsum$(EXE)_LIBS = -lm
verify_didssq$(EXE)_LIBS = -lm
verify_dizssq$(EXE)_LIBS = -lm
verify_sisssq$(EXE)_LIBS = -lm
verify_sicssq$(EXE)_LIBS = -lm
corroborate_rdgemv$(EXE)_LIBS = -lm
corroborate_rdgemm$(EXE)_LIBS = -lm
corroborate_rzgemv$(EXE)_LIBS = -lm
corroborate_rzgemm$(EXE)_LIBS = -lm
corroborate_rsgemv$(EXE)_LIBS = -lm
corroborate_rsgemm$(EXE)_LIBS = -lm
corroborate_rcgemv$(EXE)_LIBS = -lm
corroborate_rcgemm$(EXE)_LIBS = -lm
