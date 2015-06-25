TARGETS := libmpiindexed.a
SUBDIRS :=

INSTALL_LIB := $(TARGETS)

LIBMPIINDEXED := $(OBJPATH)/libmpiindexed.a

LDFLAGS += $(MPILDFLAGS)
CFLAGS += $(MPICFLAGS)

COGGED = DIDIADD.ccog \
         ZIZIADD.ccog \
         DIDIADDSQ.ccog \
         SISIADD.ccog \
         CICIADD.ccog \
         SISIADDSQ.ccog

libmpiindexed.a_DEPS = $$(LIBINDEXED) DOUBLE_INDEXED.o \
                                      DOUBLE_COMPLEX_INDEXED.o \
                                      DOUBLE_INDEXED_SCALED.o \
                                      FLOAT_INDEXED.o \
                                      FLOAT_COMPLEX_INDEXED.o \
                                      FLOAT_INDEXED_SCALED.o \
                                      DIDIADD.o \
                                      ZIZIADD.o \
                                      DIDIADDSQ.o \
                                      SISIADD.o \
                                      CICIADD.o \
                                      SISIADDSQ.o
