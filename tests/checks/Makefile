ROOT_DIR=../..
include $(ROOT_DIR)/Make.local

CFLAGS += -I$(ROOT_DIR)/$(BUILD)/include/
RBLAS     = $(ROOT_DIR)/$(BUILD)/libs/libreproblas.a $(ROOT_DIR)/$(BUILD)/libs/libindexedblas.a $(ROOT_DIR)/$(BUILD)/libs/libindexed.a
PRBLAS    = $(ROOT_DIR)/$(BUILD)/libs/libmpireproblas.a $(ROOT_DIR)/$(BUILD)/libs/libmpiindexed.a

FILE_CHECKS = validate_external_rdblas1 validate_external_rsblas1 validate_external_rzblas1 validate_external_rcblas1

VECVEC_CHECKS = validate_internal_ufp validate_internal_ufpf validate_internal_rdblas1 validate_internal_rsblas1 validate_internal_rcblas1 validate_internal_rzblas1 validate_internal_damax validate_internal_samax validate_internal_zamax validate_internal_camax

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
