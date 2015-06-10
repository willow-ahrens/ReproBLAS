TARGETS := getter$(EXE)
SUBDIRS :=

COGGED += getter.ccog

GETTER := $(OBJPATH)/getter$(EXE)

getter$(EXE)_DEPS = getter.o
