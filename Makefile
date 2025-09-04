#
# This file is part of GRAMPC - (https://github.com/grampc/grampc)
#
# GRAMPC -- A software framework for embedded nonlinear model predictive
# control using a gradient-based augmented Lagrangian approach
#
# Copyright 2014-2025 by Knut Graichen, Andreas Voelz, Thore Wietzke,
# Tobias Englert (<v2.3), Felix Mesmer (<v2.3), Soenke Rhein (<v2.3),
# Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
# All rights reserved.
#
# GRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
#


#
# compiler
#
COMPILER = gcc
AR       = ar rcs
RM       = rm -f
MKDIR    = mkdir -p
OBJEXT   = o
LIBEXT   = a
LIBNAME  = grampc


#
# path of header and source files
#
HEADER_PATH = ./include
SOURCE_PATH = ./src
LIBS_PATH   = ./libs


#
# flags and settings
#
HEADER = -I${HEADER_PATH}
CFLAGS = -c -O3 -Wall
LIBS   = 

#
# warning text to notify users of change in probfct signatures
#
$(warning Probfct signatures changed with GRAMPC 2.3. The ordering of dfd*_vec and dHdxdt changed to (out, t, x, u, p, vec, param, userparam). Make sure, your implementation uses the correct new signatures.)


#
# compile fixed-size model
#
ifeq ($(FIXEDSIZE), 1)
	LIBS_PATH = ${PROBLEM_PATH}
	LIBNAME := $(LIBNAME)_fixedsize
	HEADER += -I${PROBLEM_PATH}
	CFLAGS += -DFIXEDSIZE
endif


#
# object files and lib
#		
GRAMPC_OBJS = \
        $(SOURCE_PATH)/discrete.$(OBJEXT) \
		$(SOURCE_PATH)/ruku45.$(OBJEXT) \
		$(SOURCE_PATH)/rodas.$(OBJEXT) \
		$(SOURCE_PATH)/trapezoidal.$(OBJEXT) \
		$(SOURCE_PATH)/simpson.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_alloc.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_erk.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_fixedsize.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_init.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_mess.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_run.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_configreader.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_setopt.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_setparam.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_util.$(OBJEXT) \
		$(SOURCE_PATH)/timing.$(OBJEXT)
GRAMPC_LIB = \
		$(LIBS_PATH)/lib$(LIBNAME).$(LIBEXT)


#
# targets
#
all: $(GRAMPC_LIB)

$(GRAMPC_LIB): $(GRAMPC_OBJS)
			@echo 'Building library: $@'
			$(MKDIR) $(LIBS_PATH)
			$(AR) "$@" $(GRAMPC_OBJS)
			$(RM) $(GRAMPC_OBJS)
			@echo 'Finished building: $@'
			@echo ''

$(GRAMPC_OBJS): %.$(OBJEXT): %.c
			@echo 'Building file: $@'
			$(COMPILER) -o "$@" "$<" $(CFLAGS) $(HEADER) $(LIBS) 
			@echo 'Finished building: $@'
			@echo ''

clean: 
			@echo 'Cleaning up built files:'
			$(RM) $(GRAMPC_LIB) $(GRAMPC_OBJS)
			@echo ''

#
# end of make
#		
