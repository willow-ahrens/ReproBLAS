TARGETS :=
SUBDIRS := src tests scripts include

COGFLAGS += -D args=$(ARGS) -D params=$(TOP)/src/params.json -D mode=generate

INCLUDES += $(TOP)/include $(TOP)/scripts

INSTALL_BIN := $(TARGETS)
INSTALL_DOC := README
