TARGETS := libindexed.a foo$(EXE)
SUBDIRS :=

INSTALL_LIB := $(TARGETS)

LIBINDEXED := $(OBJPATH)/libindexed.a

libindexed.a_DEPS = diadd.o diconv.o diindex.o dimem.o dinegate.o diprint.o direnorm.o diupdate.o siadd.o siconv.o siindex.o simem.o sinegate.o siprint.o sirenorm.o siupdate.o discale.o dmdrescale.o ufp.o ufpf.o

foo$(EXE)_DEPS = foo.o $(LIBINDEXED)
