ROOT_DIR=../..
include $(ROOT_DIR)/Make.local

COMMON = test_vec.o test_opt.o test_file.o test_time.o

all : libtest.a test_data blas_inc.h

test_data: test_data.c libtest.a
	$(C_COMPILER) $(CFLAGS) test_data.c libtest.a -o test_data

%.o: %.c
	$(C_COMPILER) $(CFLAGS) -c -o $@ $<

CFLAGS += -I$(ROOT_DIR)/$(BUILD)/include/ #POSSIBLY UNNECCESARY

libtest.a: $(COMMON)
	ar -crs libtest.a $(COMMON)

blas_inc.h:
	echo "#define $(CMP_BLAS)\n" > tmp_config.h
	cat tmp_config.h _blas_inc.h > blas_inc.h
	rm -f tmp_config.h

clean:
	rm -f *.o
	rm -f *.a
	rm -f test_data
	rm -f tmp_config.h blas_inc.h
