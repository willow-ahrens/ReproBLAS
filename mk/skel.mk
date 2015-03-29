# For the reference here are some automatic variables defined by make.
# There are also their D/F variants e.g. $(<D) - check the manual.
#
# $@ - file name of the target of the rule
# $% - target member name when the target is archive member
# $< - the name of the first dependency
# $? - the names of all dependencies that are newer then the target
# $^ - the names of all dependencies

########################################################################
#                        User defined variables                        #
########################################################################

# VERB_VARS is a list of variables that you'd like to record on per
# directory level.  So if you set it to say AUTOTEST then each each
# directory will have it's own AUTOTEST_$(dir) variable with value taken
# from appropriate Rules.mk
VERB_VARS :=

# OBJ_VARS - like VERB_VARS but all values taken from Rules.mk have
# $(OBJPATH) prepended so instead of saying:
#   INSTALL_$d := $(OBJPATH)/some_target
# you can say simply
#   INSTALL := some_target
OBJ_VARS := INSTALL_BIN INSTALL_LIB PRECIOUS

# DIR_VARS - like VERB_VARS but all values taken from Rules.mk have $(d)
# prepended (unless value is an absolute path) so you can say:
#   INSTALL_DOC := Readme.txt
# instead of:
#   INSTALL_DOC_$d := $(d)/Readme.txt
DIR_VARS := INSTALL_DOC INSTALL_INC COGGED

# NOTE: There is generic macro defined below with which you can get all
# values of given variable from some subtree e.g.:
#   $(call get_subtree,INSTALL,dir)
# will give you value of all INSTALL variables from tree starting at
# 'dir'

########################################################################
#                       Directory specific flags                       #
########################################################################

# You just define in Rules.mk say
# INCLUDES_$(d) := ....
# and this will get expanded properly during compilation (see e.g. COMPILE.c)
# Of course you can still use the target specific variables if you want
# to have special setting for just one target and not the whole
# directory.  See below for definition of @RD variable.
DIR_INCLUDES = $(addprefix -I,$(INCLUDES_$(@RD)))
DIR_CPPFLAGS = $(CPPFLAGS_$(@RD))
DIR_CFLAGS = $(CFLAGS_$(@RD))
DIR_CXXFLAGS = $(CXXFLAGS_$(@RD))
DIR_COGFLAGS = $(COGFLAGS_$(@RD))

########################################################################
#                       Global flags/settings                          #
########################################################################

CFLAGS = $(DIR_CFLAGS)
CXXFLAGS = $(DIR_CXXFLAGS)

# List of includes that all (or at least majority) needs
INCLUDES :=

# Here's an example of settings for preprocessor.  -MMD is to
# automatically build dependency files as a side effect of compilation.
# This has some drawbacks (e.g. when you move/rename a file) but it is
# good enough for me.  You can improve this by using a special script
# that builds the dependency files (one can find examples on the web).
# Note that I'm adding DIR_INCLUDES before INCLUDES so that they have
# precedence.
CPPFLAGS = -MMD -D_REENTRANT -D_POSIX_C_SOURCE=200112L -D__EXTENSIONS__ \
           $(DIR_CPPFLAGS) $(DIR_INCLUDES) $(addprefix -I,$(INCLUDES))

COGFLAGS = $(DIR_COGFLAGS) $(DIR_INCLUDES) $(addprefix -I,$(INCLUDES))

# Linker flags.  The values below will use what you've specified for
# particular target or directory but if you have some flags or libraries
# that should be used for all targets/directories just append them at end.
LDFLAGS = $(LDFLAGS_$(@)) $(addprefix -L,$(LIBDIRS_$(@RD)))

# List of libraries that all targets need (either with specific command
# generated by this makefile system or for which make has built in rules
# since LDLIBS is a variable that implicit make rules are using).
# LDLIBS can be either simple or recursive, but simpler version is
# suggested :).
LDLIBS :=

########################################################################
#                       The end of generic flags                       #
########################################################################

# First we suck in user configuration ...
include $(TOP)/config.mk

# Fill BUILD_ARCH with appropriate value automatically (with some
# cosmetics in case of Cygwin/MinGW :)).
BUILD_ARCH := $(patsubst MINGW32_%,MinGW,$(patsubst CYGWIN_%,Cygwin,$(shell uname -s)))-$(shell uname -m)

# Target platform (where the code will be executed).  I'm not using
# TARGET_ARCH for this since this variable is already used in builtin
# make rules and I want to be able to use those built in rules.  By
# default we are running on the same machine we are building.
HOST_ARCH ?= $(BUILD_ARCH)

# ... host and build specific settings ...
ifneq ($(wildcard $(MK)/config-$(BUILD_ARCH)_$(HOST_ARCH).mk),)
  include $(MK)/config-$(BUILD_ARCH)_$(HOST_ARCH).mk
else
  include $(MK)/config-default.mk
endif

# Parse the configuration variables and turn them into the necessary flags
include $(MK)/config.mk

# ... optional build mode specific flags ...
ifdef BUILD_MODE
  -include $(MK)/build-$(BUILD_MODE).mk
endif


########################################################################
#         A more advanced part - if you change anything below          #
#         you should have at least vague idea how this works :D        #
########################################################################

# I define these for convenience - you can use them in your command for
# updating the target.
DEP_OBJS = $(filter %.o, $^)
DEP_OBJS? = $(filter %.o, $?)
DEP_ARCH = $(filter %.a, $^)
DEP_ARCH? = $(filter %.a, $?)
DEP_LIBS = $(addprefix -L,$(dir $(filter %.$(SOEXT), $^))) $(patsubst lib%.$(SOEXT),-l%,$(notdir $(filter %.$(SOEXT), $^)))

# Kept for backward compatibility - you should stop using these since
# I'm now not dependent on $(OBJDIR)/.fake_file any more
?R = $?
^R = $^

# Targets that match this pattern (make pattern) will use rules defined
# in:
# - def_rules.mk included below (explicit or via `skeleton' macro)
# - built in make rules
# Other targets will have to use _DEPS (and so on) variables which are
# saved in `save_vars' and used in `tgt_rule' (see below).
AUTO_TGTS := %.o

# Where to put the compiled objects.  You can e.g. make it different
# depending on the target platform (e.g. for cross-compilation a good
# choice would be OBJDIR := obj/$(HOST_ARCH)) or debugging being on/off.
#OBJDIR := $(if $(BUILD_MODE),obj/$(BUILD_MODE),obj)
OBJDIR := obj

# Convenience function to convert from a build directory back to the
# "real directory" of a target
define build_to_real_dir
$(if $(strip $(TOP_BUILD_DIR)),$(patsubst $(TOP_BUILD_DIR)%/$(OBJDIR),$(TOP)%,$(1)),$(patsubst %/$(OBJDIR),%,$(1)))
endef

# Convenience function to convert from the "real directory" to the build
# directory
define real_to_build_dir
$(if $(strip $(TOP_BUILD_DIR)),$(TOP_BUILD_DIR)$(subst $(TOP),,$(1))/$(OBJDIR),$(1)/$(OBJDIR))
endef

# By default OBJDIR is relative to the directory of the corresponding Rules.mk
# however you can use TOP_BUILD_DIR to build all objects outside of your
# project tree.  This should be an absolute path.  Note that it can be
# also inside your project like example below.
TOP_BUILD_DIR := $(if $(BUILD_MODE),$(TOP)/build/$(BUILD_MODE)/$(HOST_ARCH),$(TOP)/build/$(HOST_ARCH))
OBJPATH = $(call real_to_build_dir,$(d))
CLEAN_DIR = $(call real_to_build_dir,$(subst clean_dir_,,$@))
DIST_CLEAN_DIR = $(patsubst %/$(OBJDIR),%/$(firstword $(subst /, ,$(OBJDIR))),\
				 $(call real_to_build_dir,$(subst dist_clean_,,$@)))

# This variable contains a list of subdirectories where to look for
# sources.  That is if you have some/dir/Rules.mk where you name object
# say client.o this object will be created in some/dir/$(OBJDIR)/ and
# corresponding source file will be searched in some/dir and in
# some/dir/{x,y,z,...} where "x y z ..." is value of this variable.
SRCS_VPATH := src

# Target "real directory" - this is used above already and is most
# reliable way to refer to "per directory flags".  In theory one could
# use automatic variable already defined by make "<D" but this will not
# work well when somebody uses SRCS_VPATH variable.
@RD = $(call build_to_real_dir,$(@D))

# These are commands that are used to update the target.  If you have
# a target that make handles with built in rules just add its pattern to
# the AUTO_TGTS below.  Otherwise you have to supply the command and you
# can either do it explicitly with _CMD variable or based on the
# target's suffix and corresponding MAKECMD variable.  For example %.a
# are # updated by MAKECMD.a (exemplary setting below).  If the target
# is not filtered out by AUTO_TGTS and there's neither _CMD nor suffix
# specific command to build the target DEFAULT_MAKECMD is used.

# Creating archives gets more complicated if some dependencies are themselves
# archives. In this case, the contents have to be extracted and archived again
# in a larger archive.
EXTRACT_DIR = $@_extract
MAKECMD.a = $(call echo_cmd,AR $@) \
	$(if $(DEP_OBJS?), \
		$(AR) $(ARFLAGS) $@ $(DEP_OBJS?) \
		&&) \
	$(if $(DEP_ARCH?), \
		mkdir -p $(EXTRACT_DIR) \
		&& cd $(EXTRACT_DIR) \
		$(foreach lib,$(DEP_ARCH?),&& $(AR) xo $(lib)) \
		&& cd - \
		&& $(AR) $(ARFLAGS) $@ $(EXTRACT_DIR)/*.o \
		&& rm -rf $(EXTRACT_DIR) \
		&&) \
	$(RANLIB) $@

MAKECMD.$(SOEXT) = $(LINK.cc) $(DEP_OBJS) $(DEP_ARCH) $(DEP_LIBS) $(LIBS_$(@)) $(LDLIBS) -shared -o $@
DEFAULT_MAKECMD = $(LINK.cc) $(DEP_OBJS) $(DEP_ARCH) $(DEP_LIBS) $(LIBS_$(@)) $(LDLIBS) -o $@

########################################################################
# Below is a "Blood sugar sex^H^H^Hmake magik" :) - don't touch it     #
# unless you know what you are doing.                                  #
########################################################################

# This can be useful.  E.g. if you want to set INCLUDES_$(d) for given
# $(d) to the same value as includes for its parent directory plus some
# add ons then: INCLUDES_$(d) := $(INCLUDES_$(parent_dir)) ...
parent_dir = $(patsubst %/,%,$(dir $(d)))

define include_subdir_rules
dir_stack := $(d) $(dir_stack)
d := $(d)/$(1)
$$(eval $$(value HEADER))
include $(addsuffix /Rules.mk,$$(d))
$$(eval $$(value FOOTER))
d := $$(firstword $$(dir_stack))
dir_stack := $$(wordlist 2,$$(words $$(dir_stack)),$$(dir_stack))
endef

define save_vars
DEPS_$(1)$(2) = $(value $(2)_DEPS)
LIBS_$(1)$(2) = $(value $(2)_LIBS)
LDFLAGS_$(1)$(2) = $(value $(2)_LDFLAGS)
CMD_$(1)$(2) = $(value $(2)_CMD)
$(2)_DEPS =
$(2)_LIBS =
$(2)_LDFLAGS =
$(2)_CMD =
endef

define tgt_rule
abs_deps := $$(foreach dep,$$(DEPS_$(1)),$$(if $$(or $$(filter /%,$$(dep)),$$(filter $$$$%,$$(dep))),$$(dep),$$(addprefix $(OBJPATH)/,$$(dep))))
-include $$(addsuffix .d,$$(basename $$(abs_deps)))
$(1): $$(abs_deps) $(if $(findstring $(OBJDIR),$(1)),| $(OBJPATH),)
	$$(or $$(CMD_$(1)),$$(MAKECMD$$(suffix $$@)),$$(DEFAULT_MAKECMD))
endef

# subtree_tgts is now just a special case of a more general get_subtree
# macro since $(call get_subtree,TARGETS,dir) has the same effect but
# I'm keeping it for backward compatibility
define subtree_tgts
$(TARGETS_$(1)) $(foreach sd,$(SUBDIRS_$(1)),$(call subtree_tgts,$(sd)))
endef

define get_subtree
$($(1)_$(2)) $(foreach sd,$(SUBDIRS_$(2)),$(call get_subtree,$(1),$(sd)))
endef

# if we are using out of project build tree then there is no need to
# have dist_clean on per directory level and the one below is enough
ifneq ($(strip $(TOP_BUILD_DIR)),)
dist_clean :
	rm -rf $(TOP_BUILD_DIR)
endif

# Suck in the default rules
include $(MK)/def_rules.mk
