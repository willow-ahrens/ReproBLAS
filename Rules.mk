TARGETS :=
SUBDIRS := src tests #examples

COGFLAGS += -D args=$(ARGS) -D params=$(TOP)/src/params.json -D mode=generate

INCLUDES += $(TOP)/include

INSTALL_BIN := $(TARGETS)
INSTALL_DOC := Readme.txt
