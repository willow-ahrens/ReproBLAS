ROOT_DIR=../..
include $(ROOT_DIR)/Make.local

CFLAGS += -I$(ROOT_DIR)/$(BUILD)/include/
RBLAS     = $(ROOT_DIR)/$(BUILD)/libs/libreproblas.a $(ROOT_DIR)/$(BUILD)/libs/libindexedblas.a $(ROOT_DIR)/$(BUILD)/libs/libindexed.a
PRBLAS    = $(ROOT_DIR)/$(BUILD)/libs/libmpireproblas.a $(ROOT_DIR)/$(BUILD)/libs/libmpiindexed.a

BENCHMARKS = bench_rdsum bench_rdasum bench_rdnrm2 bench_rddot                 \
             bench_rssum bench_rsasum bench_rsnrm2 bench_rsdot                 \
             bench_rzsum bench_rdzasum bench_rdznrm2 bench_rzdotu bench_rzdotc \
             bench_rcsum bench_rscasum bench_rscnrm2 bench_rcdotu bench_rcdotc \
             bench_damax bench_damaxm                                          \
             bench_samax bench_samaxm                                          \
             bench_zamax bench_zamaxm                                          \
             bench_camax bench_camaxm                                          \
             bench_dasum bench_dnrm2 bench_ddot                                \
             bench_sasum bench_snrm2 bench_sdot                                \
             bench_dzasum bench_dznrm2 bench_zdotu bench_zdotc                 \
             bench_scasum bench_scnrm2 bench_cdotu bench_cdotc                 \
             bench_idamax                                                      \
             bench_isamax                                                      \
             bench_izamax                                                      \
             bench_icamax                                                      \

../common/libtest.a:
	@cd ../common; make

$(BENCHMARKS): $(patsubst %,%.c,$@) ../common/libtest.a
	@$(C_COMPILER) -O0 $(CFLAGS) -c $(patsubst %,%.c,$@) -o $(patsubst %,%.o,$@)
	@$(LDC) $(patsubst %,%.o,$@) -o $@ $(RBLAS)  ../common/libtest.a $(LDFLAGS)

bench:
	python ../harness.py -s bench.suite

bench-compare:
	python ../harness.py -s compare.suite

bench-flops:
	python ../harness.py -s flops.suite

bench-peak:
	python ../harness.py -s peak.suite

bench-record:
	python ../harness.py -s bench.suite -o csv > out.csv

bench-compare-record:
	python ../harness.py -s compare.suite -o csv > out.csv

bench-flops-record:
	python ../harness.py -s flops.suite -o csv > out.csv

bench-peak-record:
	python ../harness.py -s peak.suite -o csv > out.csv

clean:
	-rm -f *.o
	-rm -f $(BENCHMARKS)
