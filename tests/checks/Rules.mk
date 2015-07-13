TARGETS := validate_external_rdblas1$(EXE) validate_external_rzblas1$(EXE) validate_external_rsblas1$(EXE) validate_external_rcblas1$(EXE) \
           validate_internal_damax$(EXE) validate_internal_zamax$(EXE) validate_internal_samax$(EXE) validate_internal_camax$(EXE) \
           validate_internal_ufp$(EXE) validate_internal_ufpf$(EXE) \
           verify_rdblas1$(EXE) verify_rzblas1$(EXE) verify_rsblas1$(EXE) verify_rcblas1$(EXE) \
           validate_internal_dscale$(EXE) validate_internal_sscale$(EXE) \
           validate_internal_dindex$(EXE) validate_internal_sindex$(EXE) \
           validate_internal_dmindex$(EXE) validate_internal_smindex$(EXE) \
           verify_daugsum$(EXE) verify_zaugsum$(EXE) verify_saugsum$(EXE) verify_caugsum$(EXE) \
           verify_didssq$(EXE) verify_dizssq$(EXE) verify_sisssq$(EXE) verify_sicssq$(EXE) \
           validate_internal_daugsum$(EXE) validate_internal_zaugsum$(EXE) validate_internal_saugsum$(EXE) validate_internal_caugsum$(EXE) \
           verify_rdgemv$(EXE)

SUBDIRS :=

LDFLAGS += $(MPILDFLAGS)
CFLAGS += $(MPICFLAGS)

validate_external_rcblas1$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_external_rcblas1.o
validate_external_rdblas1$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_external_rdblas1.o
validate_external_rsblas1$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_external_rsblas1.o
validate_external_rzblas1$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_external_rzblas1.o
validate_internal_camax$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXEDBLAS) validate_internal_camax.o
validate_internal_damax$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXEDBLAS) validate_internal_damax.o
validate_internal_samax$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXEDBLAS) validate_internal_samax.o
validate_internal_ufp$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXEDBLAS) validate_internal_ufp.o
validate_internal_ufpf$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXEDBLAS) validate_internal_ufpf.o
validate_internal_zamax$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXEDBLAS) validate_internal_zamax.o
verify_rcblas1$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_rcblas1.o
verify_rdblas1$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_rdblas1.o
verify_rsblas1$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_rsblas1.o
verify_rzblas1$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_rzblas1.o
validate_internal_dscale$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXED) validate_internal_dscale.o
validate_internal_sscale$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXED) validate_internal_sscale.o
validate_internal_dindex$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXED) validate_internal_dindex.o
validate_internal_sindex$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXED) validate_internal_sindex.o
validate_internal_dmindex$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXED) validate_internal_dmindex.o
validate_internal_smindex$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXED) validate_internal_smindex.o
verify_daugsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_daugsum.o
verify_zaugsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_zaugsum.o
verify_saugsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_saugsum.o
verify_caugsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_caugsum.o
verify_didssq$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_didssq.o
verify_dizssq$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_dizssq.o
verify_sisssq$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_sisssq.o
verify_sicssq$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_sicssq.o
validate_internal_daugsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_internal_daugsum.o
validate_internal_zaugsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_internal_zaugsum.o
validate_internal_saugsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_internal_saugsum.o
validate_internal_caugsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_internal_caugsum.o
verify_rdgemv$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_rdgemv.o

#verify_prdgemv$(EXE)_DEPS = $$(LIBTEST) $$(LIBMPIREPROBLAS) verify_prdgemv.o
