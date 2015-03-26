include Make.local

all: LIB

CHECK_DIR:
	@test -d $(BUILD) || mkdir $(BUILD)
	@test -d $(BUILD)/include || mkdir $(BUILD)/include
	@test -d $(BUILD)/libs || mkdir $(BUILD)/libs

LIB : CHECK_DIR
	@cd src && ${MAKE}
	
TEST: LIB
	@cd tests && ${MAKE} all

check:
	@cd tests && ${MAKE} check --no-print-director

accuracy:
	@cd tests && ${MAKE} accuracy --no-print-director

bench:
	@cd tests && ${MAKE} bench

Indexed:
	cd src && $(MAKE) Indexed

MPIndexedFP:
	cd src && $(MAKE) MPIndexedFP

iblas:
	cd src && $(MAKE) iblas

iblas_ref:
	cd src && $(MAKE) iblas_ref

reproblas_seq:
	cd src && $(MAKE) reproblas_seq

reproblas_mpi:
	cd src && $(MAKE) reproblas_mpi

clean:
	cd src && make clean
	cd tests && make clean
