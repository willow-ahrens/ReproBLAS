TARGETS := validate_external_rdblas1$(EXE) validate_external_rzblas1$(EXE) validate_external_rsblas1$(EXE) validate_external_rcblas1$(EXE) \
           validate_internal_damax$(EXE) validate_internal_zamax$(EXE) validate_internal_samax$(EXE) validate_internal_camax$(EXE) \
           validate_internal_rdblas1$(EXE) validate_internal_rzblas1$(EXE) validate_internal_rsblas1$(EXE) validate_internal_rcblas1$(EXE) \
           validate_internal_ufp$(EXE) validate_internal_ufpf$(EXE) \
           verify_rdblas1$(EXE) verify_rzblas1$(EXE) verify_rsblas1$(EXE) verify_rcblas1$(EXE)\
           verify_rdgemv$(EXE)

SUBDIRS :=

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
verify_rdgemv$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) verify_rdgemv.o
