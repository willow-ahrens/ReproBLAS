TARGETS := validate_external_rdblas1$(EXE) validate_external_rzblas1$(EXE) validate_external_rsblas1$(EXE) validate_external_rcblas1$(EXE) \
           validate_internal_damax$(EXE) validate_internal_zamax$(EXE) validate_internal_samax$(EXE) validate_internal_camax$(EXE) \
           validate_internal_rdblas1$(EXE) validate_internal_rzblas1$(EXE) validate_internal_rsblas1$(EXE) validate_internal_rcblas1$(EXE) \
           validate_internal_ufp$(EXE) validate_internal_ufpf$(EXE) \
           verify_rdblas1$(EXE) verify_rzblas1$(EXE) verify_rsblas1$(EXE) verify_rcblas1$(EXE) \
           verify_dindex$(EXE) verify_sindex$(EXE) \
           verify_dmindex$(EXE) verify_smindex$(EXE) \
           verify_daugsum$(EXE) verify_zaugsum$(EXE) verify_saugsum$(EXE) verify_caugsum$(EXE) \
           validate_internal_daugsum$(EXE)

SUBDIRS :=

LDFLAGS += $(MPILDFLAGS)
CFLAGS += $(MPICFLAGS)

validate_external_rcblas1$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_external_rcblas1.o
validate_external_rdblas1$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_external_rdblas1.o
validate_external_rsblas1$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_external_rsblas1.o
validate_external_rzblas1$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_external_rzblas1.o
validate_internal_camax$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXEDBLAS) validate_internal_camax.o
validate_internal_damax$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXEDBLAS) validate_internal_damax.o
validate_internal_rcblas1$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_internal_rcblas1.o
validate_internal_rdblas1$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_internal_rdblas1.o
validate_internal_rsblas1$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_internal_rsblas1.o
validate_internal_rzblas1$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_internal_rzblas1.o
validate_internal_samax$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXEDBLAS) validate_internal_samax.o
validate_internal_ufp$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXEDBLAS) validate_internal_ufp.o
validate_internal_ufpf$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXEDBLAS) validate_internal_ufpf.o
validate_internal_zamax$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXEDBLAS) validate_internal_zamax.o
verify_rcblas1$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_rcblas1.o
verify_rdblas1$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_rdblas1.o
verify_rsblas1$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_rsblas1.o
verify_rzblas1$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_rzblas1.o
verify_dindex$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXED) verify_dindex.o
verify_sindex$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXED) verify_sindex.o
verify_dmindex$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXED) verify_dmindex.o
verify_smindex$(EXE)_DEPS = $$(LIBTEST) $$(LIBINDEXED) verify_smindex.o
verify_daugsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_daugsum.o
verify_zaugsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_zaugsum.o
verify_saugsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_saugsum.o
verify_caugsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_caugsum.o
validate_internal_daugsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) validate_internal_daugsum.o

#verify_rdgemv$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_rdgemv.o
#verify_prdgemv$(EXE)_DEPS = $$(LIBTEST) $$(LIBMPIREPROBLAS) verify_prdgemv.o
