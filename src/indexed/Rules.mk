TARGETS := libindexed.a
SUBDIRS := 

INSTALL_LIB := $(TARGETS)

LIBINDEXED := $(OBJPATH)/libindexed.a

libindexed.a_DEPS = dconv.o ufp.o ufpf.o dbound.o dIAdd.o dIRenorm.o dIUpdate.o sbound.o sIAdd.o sconv.o sIRenorm.o sIUpdate.o diprint.o ziprint.o siprint.o ciprint.o
