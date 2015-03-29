TARGETS :=
SUBDIRS := indexed indexedBLAS reproBLAS MPI_indexed MPI_reproBLAS

INCLUDES += $(d)/gen
ReproBLAS_COGFLAGS := -D args=$(d)/default_args.json -D params=$(d)/params.json
COGFLAGS += $(ReproBLAS_COGFLAGS)
