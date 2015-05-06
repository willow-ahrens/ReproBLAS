TARGETS := libtest.a test_data$(EXE)
SUBDIRS := 

test_data$(EXE)_DEPS = libtest.a test_data.o

LIBTEST := $(OBJPATH)/libtest.a

libtest.a_DEPS = test_file.o test_opt.o test_time.o test_util.o test_metric.o
