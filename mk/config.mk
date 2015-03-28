# Fill BUILD_ARCH with appropriate value automatically (with some
# cosmetics in case of Cygwin/MinGW :)).
BUILD_ARCH := $(patsubst MINGW32_%,MinGW,$(patsubst CYGWIN_%,Cygwin,$(shell uname -s)))-$(shell uname -m)

# Target platform (where the code will be executed).  I'm not using
# TARGET_ARCH for this since this variable is already used in builtin
# make rules and I want to be able to use those built in rules.  By
# default we are running on the same machine we are building.
HOST_ARCH := $(BUILD_ARCH)

# Build mode e.g. debug, profile, release.  Build specific mode flags
# can be entered in $(MK)/build-$(BUILD_MODE).mk file e.g. for debug
# following seems to be a reasonable contents
#CFLAGS   += -ggdb
#CXXFLAGS += -ggdb
#CPPFLAGS += -DDEBUG
#LDFLAGS  += -ggdb
# If you don't plan on having different build modes then just comment
# below or set it to empty.
BUILD_MODE := release

# To have some per directory setting automatically propagated to all
# subdirs then uncomment below.  That way you can have all project
# compiled with "global" settings and easily switch flags for some
# subtree just by setting per directory settings at the top dir of the
# subtree.  You may of course overwrite inherited values and you can
# turn inheritance in some part by just clearing INHERIT_DIR_VARS_$(d)
# This is a global inheritance flag - you might want to turn it on only
# in some directory (just set INHERIT_DIR_VARS_$(d) there).
INHERIT_DIR_VARS := CPPFLAGS INCLUDES CFLAGS CXXFLAGS

# Default optimization flags.
OPTFLAGS := -O3

