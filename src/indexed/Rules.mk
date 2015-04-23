TARGETS := libindexed.a addtest$(EXE)
SUBDIRS := 

INSTALL_LIB := $(TARGETS)

LIBINDEXED := $(OBJPATH)/libindexed.a

libindexed.a_DEPS = dconv.o ufp.o ufpf.o dbound.o dIAdd.o dIRenorm.o dupdate.o sbound.o sIAdd.o sconv.o sIRenorm.o supdate.o diprint.o ziprint.o siprint.o ciprint.o

addtest$(EXE)_DEPS = addtest.o $$(LIBINDEXED)
