DSANITY = sanity_check0 sanity_check1 sanity_check2 sanity_check3 \
	sanity_check4	\
	sanity_check5 	\
	sanity_check6	\
	sanity_check7

ZSANITY = sanity_check_z4 sanity_check_z5 sanity_check_z6 sanity_check_z7

CSANITY = sanity_check_c1 sanity_check_c3 sanity_check_c4 sanity_check_c5 sanity_check_c7

SSANITY = sanity_check_s0 sanity_check_s1 sanity_check_s2 sanity_check_s3 \
	sanity_check_s4 sanity_check_s5 sanity_check_s6 sanity_check_s7

dcheck_reproducibility.o zcheck_reproducibility.o \
scheck_reproducibility.o ccheck_reproducibility.o \
sanity_check_d.o sanity_check_s.o sanity_check_z.o sanity_check_c.o: %.o :%.c
	@${C_COMPILER} $(CFLAGS)  -c $< -o $@

ifeq (,$N)
N=4096
endif

$(DSANITY): sanity_check_d.o dgenvec.o debug.o dcheck_reproducibility.o
	@$(C_COMPILER) $(CFLAGS) -c $(patsubst %,%.c,$@) -o sanity_check_sub.o
	@$(LDC) sanity_check_d.o dgenvec.o debug.o sanity_check_sub.o dcheck_reproducibility.o -o sanity_check.out $(RBLAS) $(LDFLAGS)
	@./sanity_check.out -n $N
	@-rm sanity_check.out sanity_check_sub.o

$(SSANITY): sanity_check_s.o sgenvec.o debug.o scheck_reproducibility.o
	@$(C_COMPILER) $(CFLAGS) -c $(patsubst %,%.c,$@) -o sanity_check_sub.o
	@$(LDC) sanity_check_s.o sgenvec.o debug.o sanity_check_sub.o scheck_reproducibility.o -o sanity_check.out $(RBLAS) $(LDFLAGS)
	@./sanity_check.out -n $N
	@-rm sanity_check.out sanity_check_sub.o

$(ZSANITY): sanity_check_z.o zgenvec.o dgenvec.o debug.o zcheck_reproducibility.o
	@$(C_COMPILER) $(CFLAGS) -c $(patsubst %,%.c,$@) -o sanity_check_sub.o
	@$(LDC) sanity_check_z.o zgenvec.o dgenvec.o debug.o sanity_check_sub.o zcheck_reproducibility.o -o sanity_check.out $(RBLAS) $(LDFLAGS)
	@./sanity_check.out -n $N
	@-rm sanity_check.out sanity_check_sub.o

$(CSANITY): sanity_check_c.o cgenvec.o sgenvec.o debug.o ccheck_reproducibility.o
	@$(C_COMPILER) $(CFLAGS) -c $(patsubst %,%.c,$@) -o sanity_check_sub.o
	@$(LDC) sanity_check_c.o cgenvec.o sgenvec.o \
		debug.o sanity_check_sub.o ccheck_reproducibility.o \
	-o sanity_check.out $(RBLAS) $(LDFLAGS)
	@./sanity_check.out -n $N
	@-rm sanity_check.out sanity_check_sub.o

check_d_header:
	@echo '+---------+--------------------------------------------------+--------+'
	@echo '| CHECK REPRODUCIBILITY OF RBLAS LIB (DOUBLE)                | STATUS |'
	@echo '+---------+--------------------------------------------------+--------+'

check_s_header:
	@echo '+---------+--------------------------------------------------+--------+'
	@echo '| CHECK REPRODUCIBILITY OF RBLAS LIB (SINGLE)                | STATUS |'
	@echo '+---------+--------------------------------------------------+--------+'

check_z_header:
	@echo '+---------+--------------------------------------------------+--------+'
	@echo '| CHECK REPRODUCIBILITY OF RBLAS LIB (DCOMPLEX)              | STATUS |'
	@echo '+---------+--------------------------------------------------+--------+'

check_c_header:
	@echo '+---------+--------------------------------------------------+--------+'
	@echo '| CHECK REPRODUCIBILITY OF RBLAS LIB (SCOMPLEX)              | STATUS |'
	@echo '+---------+--------------------------------------------------+--------+'

sanity_check_d: check_d_header $(DSANITY)
	@echo '+---------+--------------------------------------------------+--------+'
	@-rm dcheck_reproducibility.o sanity_check_d.o
	
sanity_check_s: check_s_header $(SSANITY)
	@echo '+---------+--------------------------------------------------+--------+'
	@-rm scheck_reproducibility.o sanity_check_s.o

sanity_check_z: check_z_header $(ZSANITY)
	@echo '+---------+--------------------------------------------------+--------+'
	@-rm zcheck_reproducibility.o sanity_check_z.o

sanity_check_c: check_c_header $(CSANITY)
	@echo '+---------+--------------------------------------------------+--------+'
	@-rm ccheck_reproducibility.o sanity_check_c.o

check: sanity_check_d sanity_check_s sanity_check_z sanity_check_c
	@-rm debug.o dgenvec.o sgenvec.o zgenvec.o cgenvec.o

