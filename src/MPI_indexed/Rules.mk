TARGETS := libmpiindexed.a
SUBDIRS := 

INSTALL_LIB := $(TARGETS)

LIBMPIINDEXED := $(OBJPATH)/libmpiindexed.a

LDFLAGS += $(MPILDFLAGS)
CFLAGS += $(MPICFLAGS)

libmpiindexed.a_DEPS = $$(LIBINDEXED) #MPI_indexed.o MPI_dIndexed.o MPI_sIndexed.o
