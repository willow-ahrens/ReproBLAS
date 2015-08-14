TARGETS := validate_internal_damax$(EXE) validate_internal_zamax$(EXE) validate_internal_samax$(EXE) validate_internal_camax$(EXE) \
           validate_internal_ufp$(EXE) validate_internal_ufpf$(EXE) \
           validate_internal_dscale$(EXE) validate_internal_sscale$(EXE) \
           validate_internal_dindex$(EXE) validate_internal_sindex$(EXE) \
           validate_internal_dmindex$(EXE) validate_internal_smindex$(EXE) \
           verify_daugsum$(EXE) verify_zaugsum$(EXE) verify_saugsum$(EXE) verify_caugsum$(EXE) \
           validate_internal_daugsum$(EXE) validate_internal_zaugsum$(EXE) validate_internal_saugsum$(EXE) validate_internal_caugsum$(EXE)\
           verify_didssq$(EXE) verify_dizssq$(EXE) verify_sisssq$(EXE) verify_sicssq$(EXE) \
           corroborate_rdgemv$(EXE) \
           corroborate_rdgemm$(EXE) \
           corroborate_rzgemv$(EXE) \

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
verify_didssq$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_didssq.o
verify_dizssq$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_dizssq.o
verify_sisssq$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_sisssq.o
verify_sicssq$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_sicssq.o
corroborate_rdgemv$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) corroborate_rdgemv.o
corroborate_rdgemm$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) corroborate_rdgemm.o
corroborate_rzgemv$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) corroborate_rzgemv.o
