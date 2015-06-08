TARGETS := get_vectorization$(EXE) get_max_fold$(EXE) get_default_fold$(EXE)
SUBDIRS :=

COGGED += get_vectorization.ccog

get_vectorization$(EXE)_DEPS = get_vectorization.o
get_max_fold$(EXE)_DEPS = get_max_fold.o
get_default_fold$(EXE)_DEPS = get_default_fold.o
