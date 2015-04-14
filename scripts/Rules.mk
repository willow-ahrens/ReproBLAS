TARGETS := get_vectorization$(EXE)
SUBDIRS :=

COGGED += get_vectorizations.c

get_vectorization$(EXE)_DEPS = get_vectorization.o
