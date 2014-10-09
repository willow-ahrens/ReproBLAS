common/libtest.a:
	@cd common; make

check: common/libtest.a
	cd checks; make check
