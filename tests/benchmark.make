ifeq ($(shell uname), Darwin)
system_info:
	@echo '  CPU FREQUENCY: ' ${patsubst hw.cpufrequency:,,$(shell sysctl hw.cpufrequency)}
	@echo '  L1 CACHE     : ' ${patsubst hw.l1dcachesize:,,$(shell sysctl hw.l1dcachesize)}
	@echo '  L2 CACHE     : ' ${patsubst hw.l2cachesize:,,$(shell sysctl hw.l2cachesize)}
	@echo '  L3 CACHE     : ' ${patsubst hw.l3cachesize:,,$(shell sysctl hw.l3cachesize)}
else
ifeq ($(shell uname), Linux)
system_info:
	@echo ' ' $(shell grep -m 1 ^'cpu MHz' /proc/cpuinfo)
	@echo ' L1 CACHE: ' $(shell cat /sys/devices/system/cpu/cpu0/cache/index1/size)
	@echo ' L2 CACHE: ' $(shell cat /sys/devices/system/cpu/cpu0/cache/index2/size)
	@echo ' L3 CACHE: ' $(shell cat /sys/devices/system/cpu/cpu0/cache/index3/size)
else
system_info:
endif
endif

CMP_BLAS=NO_BLAS

ifneq (,$(BLAS))
	CMP_BLAS = $(BLAS)
endif

benchmark_header:
	@echo '+----------------------------------------------------+'
	@echo '| BENCHMARKING REPRODUCIBLE BLAS (MHz)               |'
	@echo '+----------------------------------------------------+'


BENCHMARK_RBLAS = rdblas1 rsblas1 rzblas1 rcblas1

$(BENCHMARK_RBLAS): debug.o dgenvec.o sgenvec.o
	@echo '#define ' $(CMP_BLAS) > tmp_config.h
	@$(C_COMPILER) $(CFLAGS) -c $(patsubst %,%.c,$@) -o $@.o
	@$(LDC) $@.o debug.o dgenvec.o sgenvec.o -o benchmark-rblas $(RBLAS) $(LDFLAGS)
	@./benchmark-rblas --header -n 128:128:256 --cmp
	@./benchmark-rblas -n 512 --cmp
	@./benchmark-rblas -n 1024:1024:4096 --cmp
	@-rm $@.o benchmark-rblas

benchmark: benchmark_header system_info $(BENCHMARK_RBLAS)
	@-rm debug.o dgenvec.o sgenvec.o
