TARGETS := acc_rdsum$(EXE) acc_rzsum$(EXE) acc_rssum$(EXE) acc_rcsum$(EXE)

ifneq ($(BLAS),)
TARGETS += acc_dsum$(EXE) acc_zsum$(EXE) acc_ssum$(EXE) acc_csum$(EXE)
endif

SUBDIRS :=

acc_rcsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_rcsum.o
acc_rdsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_rdsum.o
acc_rssum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_rssum.o
acc_rzsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_rzsum.o
acc_csum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_csum.o
acc_dsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_dsum.o
acc_ssum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_ssum.o
acc_zsum$(EXE)_DEPS = $$(LIBTEST) $$(LIBREPROBLAS) acc_zsum.o

acc_rcsum$(EXE)_LIBS = -lm
acc_rdsum$(EXE)_LIBS = -lm
acc_rssum$(EXE)_LIBS = -lm
acc_rzsum$(EXE)_LIBS = -lm
acc_csum$(EXE)_LIBS = -lm
acc_dsum$(EXE)_LIBS = -lm
acc_ssum$(EXE)_LIBS = -lm
acc_zsum$(EXE)_LIBS = -lm
