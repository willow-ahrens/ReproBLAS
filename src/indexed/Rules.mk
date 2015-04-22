TARGETS := libindexed.a
SUBDIRS := 

INSTALL_LIB := $(TARGETS)

LIBINDEXED := $(OBJPATH)/libindexed.a

libindexed.a_DEPS = dconv.o ufp.o ufpf.o dIndexed.o dIAdd.o dIRenorm.o dIUpdate.o sIndexed.o sIAdd.o sconv.o sIRenorm.o sIUpdate.o
