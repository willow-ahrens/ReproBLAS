TARGETS := libtest.a test_data$(EXE)
SUBDIRS := xblas

test_data$(EXE)_DEPS = libtest.a test_data.o
test_data$(EXE)_LIBS = -lm

LIBTEST := $(OBJPATH)/libtest.a

libtest.a_DEPS = test_file.o test_opt.o test_time.o test_util.o test_metric.o $$(LIBTESTXBLAS)
