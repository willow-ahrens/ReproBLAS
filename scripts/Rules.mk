TARGETS := get_vectorization$(EXE) \
           get_max_fold$(EXE) \
           get_default_fold$(EXE) \
           get_diendurance$(EXE) \
           get_siendurance$(EXE)
SUBDIRS :=

COGGED += get_vectorization.ccog

get_vectorization$(EXE)_DEPS = get_vectorization.o
get_max_fold$(EXE)_DEPS = get_max_fold.o
get_default_fold$(EXE)_DEPS = get_default_fold.o
get_diendurance$(EXE)_DEPS = get_siendurance.o
get_siendurance$(EXE)_DEPS = get_siendurance.o
