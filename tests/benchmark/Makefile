ROOT_DIR=../..
include $(ROOT_DIR)/Make.local

CFLAGS += -I$(ROOT_DIR)/$(BUILD)/include/
RBLAS     = $(ROOT_DIR)/$(BUILD)/libs/libreproblas_seq.a $(ROOT_DIR)/$(BUILD)/libs/libiblas.a $(ROOT_DIR)/$(BUILD)/libs/libindexedfp.a
RBLAS_REF = $(ROOT_DIR)/$(BUILD)/libs/libreproblas_seq.a $(ROOT_DIR)/$(BUILD)/libs/libiblas_ref.a $(ROOT_DIR)/$(BUILD)/libs/libindexedfp.a
PRBLAS    = $(ROOT_DIR)/$(BUILD)/libs/libreproblas_mpi.a $(ROOT_DIR)/$(BUILD)/libs/libmpindexedfp.a

BENCHMARKS = bench_rdsum bench_rdasum bench_rdnrm2 bench_rddot                 \
             bench_rssum bench_rsasum bench_rsnrm2 bench_rsdot                 \
             bench_rzsum bench_rdzasum bench_rdznrm2 bench_rzdotu bench_rzdotc \
             bench_rcsum bench_rscasum bench_rscnrm2 bench_rcdotu bench_rcdotc

../common/libtest.a:
	@cd ../common; make

$(BENCHMARKS): $(patsubst %,%.c,$@) ../common/libtest.a
	@$(C_COMPILER) $(CFLAGS) -c $(patsubst %,%.c,$@) -o $(patsubst %,%.o,$@)
	@$(LDC) $(patsubst %,%.o,$@) -o $@ $(RBLAS)  ../common/libtest.a $(LDFLAGS)

bench:
	python ../harness.py -s bench.suite

clean:
	-rm -f *.o
	-rm -f $(BENCHMARKS)
