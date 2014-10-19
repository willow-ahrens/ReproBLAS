ROOT_DIR=../..
include $(ROOT_DIR)/Make.local

CFLAGS += -I$(ROOT_DIR)/$(BUILD)/include/
RBLAS     = $(ROOT_DIR)/$(BUILD)/libs/libreproblas_seq.a $(ROOT_DIR)/$(BUILD)/libs/libiblas.a $(ROOT_DIR)/$(BUILD)/libs/libindexedfp.a
RBLAS_REF = $(ROOT_DIR)/$(BUILD)/libs/libreproblas_seq.a $(ROOT_DIR)/$(BUILD)/libs/libiblas_ref.a $(ROOT_DIR)/$(BUILD)/libs/libindexedfp.a
PRBLAS    = $(ROOT_DIR)/$(BUILD)/libs/libreproblas_mpi.a $(ROOT_DIR)/$(BUILD)/libs/libmpindexedfp.a

FILE_CHECKS = validate_external_rdblas1 validate_external_rsblas1 validate_external_rzblas1 validate_external_rcblas1

VECVEC_CHECKS = validate_internal_ufp validate_internal_ufpf validate_internal_rdblas1 validate_internal_rsblas1 validate_internal_rcblas1 validate_internal_rzblas1

VECVEC_FILL_CHECKS = verify_rdblas1 verify_rsblas1 verify_rcblas1 verify_rzblas1

CHECKS = $(FILE_CHECKS) $(VECVEC_CHECKS) $(VECVEC_FILL_CHECKS)

../common/libtest.a:
	@cd ../common; make

$(CHECKS): $(patsubst %,%.c,$@) ../common/libtest.a
	@$(C_COMPILER) $(CFLAGS) -c $(patsubst %,%.c,$@) -o $(patsubst %,%.o,$@)
	@$(LDC) $(patsubst %,%.o,$@) -o $@ $(RBLAS)  ../common/libtest.a $(LDFLAGS)

check:
	python ../harness.py -s check.suite

reference:
	python ../harness.py -s data/create.suite
	python ../harness.py -s reference.suite

clean:
	-rm -f *.o
	-rm -f $(CHECKS)
