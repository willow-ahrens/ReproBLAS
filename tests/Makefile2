ROOT_DIR=..
include $(ROOT_DIR)/Make.local

CFLAGS += -I$(ROOT_DIR)/$(BUILD)/include/

common/libtest.a:
	@cd common; make

check: common/libtest.a
	cd checks; make check

bench: common/libtest.a
	cd benchs; make bench

clean:
	rm -rf *.o run-*
	rm -rf *.o
	rm -rf run-*
	cd checks; make clean
	cd benchs; make clean
	cd common; make clean
