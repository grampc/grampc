#
# This file is part of GRAMPC - (https://sourceforge.net/projects/grampc/)
#
# GRAMPC -- A software framework for embedded nonlinear model predictive
# control using a gradient-based augmented Lagrangian approach
#
# Copyright 2014-2019 by Tobias Englert, Knut Graichen, Felix Mesmer,
# Soenke Rhein, Andreas Voelz, Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
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
		$(SOURCE_PATH)/euler1.$(OBJEXT) \
		$(SOURCE_PATH)/eulermod2.$(OBJEXT) \
		$(SOURCE_PATH)/heun2.$(OBJEXT) \
		$(SOURCE_PATH)/ruku45.$(OBJEXT) \
		$(SOURCE_PATH)/rodas.$(OBJEXT) \
		$(SOURCE_PATH)/trapezodial.$(OBJEXT) \
		$(SOURCE_PATH)/simpson.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_alloc.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_fixedsize.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_init.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_mess.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_run.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_setopt.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_setparam.$(OBJEXT) \
		$(SOURCE_PATH)/grampc_util.$(OBJEXT) 
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
