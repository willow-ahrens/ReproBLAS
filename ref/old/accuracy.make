RBLAS_ACC = rdblas1_acc rsblas1_acc rzblas1_acc rcblas1_acc

ifeq (,$(N))
	N = 1024
endif

CMP_BLAS=NO_BLAS

ifneq (,$(BLAS))
	CMP_BLAS = $(BLAS)
endif

$(RBLAS_ACC): debug.o dgenvec.o sgenvec.o
	@echo '#define ' $(CMP_BLAS) > tmp_config.h
	@echo ''
	@$(C_COMPILER) $(CFLAGS) -c $(patsubst %,%.c,$@) -o acc_check.o
	@$(LDC) acc_check.o dgenvec.o sgenvec.o debug.o -o acc_check.out $(RBLAS) $(LDFLAGS)
	@./acc_check.out -n $(N) -d 0 --print-header
	@./acc_check.out -n $(N) -d 1
	@./acc_check.out -n $(N) -d 2 
	@./acc_check.out -n $(N) -d 3 
	@./acc_check.out -n $(N) -d 4 
	@./acc_check.out -n $(N) -K 1e2
	@./acc_check.out -n $(N) -K 1e4 
	@./acc_check.out -n $(N) -K 1e6 
	@./acc_check.out -n $(N) -K 1e8 
	@./acc_check.out -n $(N) -K 1e12 
	@./acc_check.out -n $(N) -K 1e15 
	@-rm acc_check.out acc_check.o

acc_header :
	@echo '+--------------------------------------------------------------------+'
	@echo '| CHECK ACCURACY OF RBLAS                                            |'
	@echo '+--------------------------------------------------------------------+'

accuracy: acc_header $(RBLAS_ACC)
	@-rm debug.o dgenvec.o sgenvec.o

