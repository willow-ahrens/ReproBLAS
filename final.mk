# Just a simple example how final.mk can be used for 'install' targets
# You can refer here to any variable defined in the project tree since
# it is included after all rules has been read and processed.
#
# Variables of particular interest:
# INSTALL_BIN_$(dir) - binaries to be installed from directory 'dir'
# INSTALL_LIB_$(dir) - same for libraries
# INSTALL_DOC_$(dir) - and for documentation

BIN_DIR := /tmp/test-ex1/bin
LIB_DIR := /tmp/test-ex1/lib
DOC_DIR := /tmp/test-ex1/share/doc/ex1

INSTALL := install
INSTALL_DATA := install -m 644

install: install-bin install-lib install-doc

install-bin : $(call get_subtree,INSTALL_BIN,$(TOP))
	$(INSTALL) -d $(BIN_DIR)
	$(INSTALL) -t $(BIN_DIR) $^

install-lib : $(call get_subtree,INSTALL_LIB,$(TOP))
	$(INSTALL) -d $(LIB_DIR)
	$(INSTALL) -t $(LIB_DIR) $(filter-out %.a,$^)
	$(INSTALL) -t $(LIB_DIR) -m 644 $(filter %.a,$^)

install-doc: $(call get_subtree,INSTALL_DOC,$(TOP))
	$(INSTALL) -d $(DOC_DIR)
	$(INSTALL_DATA) -t $(DOC_DIR) $^

.PHONY: params replace excise check bench reference tune doc

# Creates the parameter list
params:
	rm -f $(TOP)/src/params.json
	$(foreach SOURCE, $(call get_subtree,COGGED,$(TOP)), $(COG) -D mode=params $(SOURCE);)

# Creates default arguments from current parameter list
default_args:
	$(CALL_PYTHON) $(TOP)/src/gen/default_args.py --params $(TOP)/src/params.json --args $(TOP)/src/default_args.json

# tunes
tune:
	$(CALL_PYTHON) $(TOP)/tune/ReproBLASOpenTuner.py --params $(TOP)/src/params.json --args $(TOP)/src/tuned_args.json --database $(TOP)/tune/ReproBLASOpenTuner.db --trials 50 --no-dups --verbose $(VERBOSE)

# Runs code generators in place
replace:
	$(foreach SOURCE, $(call get_subtree,COGGED,$(TOP)), $(COG) -r $(SOURCE) &&) echo

# Removes generated code from code generators
excise:
	$(foreach SOURCE, $(call get_subtree,COGGED,$(TOP)), $(COG) -r -x $(SOURCE) &&) echo

check:
	$(CALL_PYTHON) $(TOP)/tests/checks/check.py --runmode parallel --verbose $(VERBOSE)

reference:
	rm $(TOP)/tests/checks/data/*
	$(CALL_PYTHON) $(TOP)/tests/checks/reference.py

bench:
	$(CALL_PYTHON) $(TOP)/tests/benchs/bench.py --verbose $(VERBOSE)

doc:
	cd $(TOP); rm -rf doc/*; doxygen config.dox; cd doc/latex;make
