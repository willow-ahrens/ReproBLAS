TARGETS := libindexed.a
SUBDIRS := 

INSTALL_LIB := $(TARGETS)

LIBINDEXED := $(OBJPATH)/libindexed.a

libindexed.a_DEPS = dIConv.o dIndexed.o dIAdd.o dIRenorm.o dIUpdate.o sIndexed.o sIAdd.o sIConv.o sIRenorm.o sIUpdate.o
