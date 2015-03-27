TARGETS :=
SUBDIRS := src tests #examples

INCLUDES += $(TOP)/include
CFLAGS += -march=core-avx-i

INSTALL_BIN := $(TARGETS)
INSTALL_DOC := Readme.txt
