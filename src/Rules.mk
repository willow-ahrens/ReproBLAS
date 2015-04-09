TARGETS := params.json
SUBDIRS := indexed indexedBLAS reproBLAS MPI_indexed MPI_reproBLAS

COGFLAGS += -D args=$(ARGS) -D params=$(OBJPATH)/params.json -D mode=generate

params.json_DEPS = $$(call get_subtree,COGGED,$(TOP))
params.json_CMD = $(if $(PYTHON), rm -f $(OBJPATH)/params.json; $(foreach SOURCE, $(call get_subtree,COGGED,$(TOP)), $(COG) -D mode=params $(SOURCE);))

INCLUDES += $(d)/gen
