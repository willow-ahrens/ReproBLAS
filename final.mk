# Just a simple example how final.mk can be used for 'install' targets
# You can refer here to any variable defined in the project tree since
# it is included after all rules has been read and processed.
#
# Variables of particular interest:
# INSTALL_BIN_$(dir) - binaries to be installed from directory 'dir'
# INSTALL_LIB_$(dir) - same for libraries
# INSTALL_INC_$(dir) - same for headers
# INSTALL_DOC_$(dir) - and for documentation

INSTALL := install
INSTALL_DATA := install -m 644

.PHONY: check bench acc reference tune dox excise update install install-bin install-lib install-inc install-doc

install: install-bin install-lib install-inc install-doc

install-bin : $(call get_subtree,INSTALL_BIN,$(TOP))
	$(INSTALL) -d $(BIN_DIR)
ifneq ($(strip $(call get_subtree,INSTALL_BIN,$(TOP))),)
	$(INSTALL) $^ $(BIN_DIR)
endif

install-lib : $(call get_subtree,INSTALL_LIB,$(TOP))
	$(INSTALL) -d $(LIB_DIR)
ifneq ($(strip $(filter-out %.a,$(call get_subtree,INSTALL_LIB,$(TOP)))),)
	$(INSTALL) $(filter-out %.a,$^) $(LIB_DIR)
endif
ifneq ($(strip $(filter %.a,$(call get_subtree,INSTALL_LIB,$(TOP)))),)
	$(INSTALL_DATA) $(filter %.a,$^) $(LIB_DIR)
endif

install-inc: $(call get_subtree,INSTALL_INC,$(TOP))
	$(INSTALL) -d $(INC_DIR)
ifneq ($(strip $(call get_subtree,INSTALL_INC,$(TOP))),)
	$(INSTALL_DATA) $^ $(INC_DIR)
endif

install-doc: $(call get_subtree,INSTALL_DOC,$(TOP))
	$(INSTALL) -d $(DOC_DIR)
ifneq ($(strip $(call get_subtree,INSTALL_DOC,$(TOP))),)
	$(INSTALL_DATA) $^ $(DOC_DIR)
endif

# tunes
tune:# update
	$(CALL_PYTHON) $(TOP)/tune/ReproBLASOpenTuner.py --params $(TOP)/src/params.json --args $(TOP)/src/tuned_args.json --database $(TOP)/tune/ReproBLASOpenTuner.db --trials 100 --no-dups --verbose $(VERBOSE) --bail-threshold 7
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

acc:
	$(CALL_PYTHON) $(TOP)/tests/accs/acc.py --runmode parallel --verbose $(VERBOSE)

dox:
	cd $(TOP); rm -rf dox; doxygen config.dox; cd dox/latex; make
	cd $(TOP); cd dox; zip -r reference_website.zip html
	cd $(TOP); rm -f doc/reference_website.zip; rm -f doc/reference_manual.pdf
	cd $(TOP); cp dox/reference_website.zip doc; cp dox/latex/refman.pdf doc/reference_manual.pdf
	cd $(TOP); rm -rf dox;

update: $(GETTER)
	$(GETTER) > $(TOP)/scripts/getter.json
	rm -f $(TOP)/src/params.json
ifeq ($(VERBOSE), true)
	@$(foreach SOURCE, $(call get_subtree,COGGED,$(TOP)), echo "$(COG) $(SOURCE)"; $(COG) -D mode=params $(SOURCE) > $(DEVNULL);)
else
	@$(foreach SOURCE, $(call get_subtree,COGGED,$(TOP)), echo "COG $(SOURCE)"; $(COG) -D mode=params $(SOURCE) > $(DEVNULL);)
endif
	$(CALL_PYTHON) $(TOP)/src/gen/default_args.py --params $(TOP)/src/params.json --args $(TOP)/src/default_args.json
ifeq ($(VERBOSE), true)
	@$(foreach SOURCE, $(call get_subtree,COGGED,$(TOP)), echo "$(COG) $(SOURCE)"; $(COG) -r $(SOURCE);)
else
	@$(foreach SOURCE, $(call get_subtree,COGGED,$(TOP)), echo "COG $(SOURCE)"; $(COG) -r $(SOURCE) > $(DEVNULL);)
endif
