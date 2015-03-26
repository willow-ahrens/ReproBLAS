TARGETS := libmpiindexed.a
SUBDIRS := 

INSTALL_LIB := $(TARGETS)

LIBMPIINDEXED := $(OBJPATH)/libmpiindexed.a

libmpiindexed.a_DEPS = $$(LIBINDEXED) MPI_indexed.o MPI_dIndexed.o MPI_sIndexed.o
