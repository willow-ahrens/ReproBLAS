TARGETS := libmpireproblas.a
SUBDIRS := 

INSTALL_LIB := $(TARGETS)

LDFLAGS += $(MPILDFLAGS)
CFLAGS += $(MPICFLAGS)

LIBMPIREPROBLAS := $(OBJPATH)/libmpireproblas.a

libmpireproblas.a_DEPS = $$(LIBINDEXED) $$(LIBINDEXEDBLAS) $$(LIBMPIINDEXED) \
                         prdasum.o  prdsum.o prddot.o prdnrm2.o              \
                         prdzasum.o przsum.o przdot.o prdznrm2.o             \
                         prsasum.o  prssum.o prsdot.o prsnrm2.o              \
                         prscasum.o prcsum.o prcdot.o prscnrm2.o             \
                         prdgemv.o
