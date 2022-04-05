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
# input and output files
#
SRCS = main_CRANE_3D.c \
		probfct_CRANE_3D.c
OUTPUT = startMPC


#
# compiler
#
COMPILER    = gcc
RM          = rm -f
MKDIR       = mkdir -p
OBJEXT      = o
LIBEXT      = a
lib         = -l
PATH_PREFIX = ../..


#
# path of header and lib files
#
HEADER_PATH = $(PATH_PREFIX)/include
LIBS_PATH   = $(PATH_PREFIX)/libs


#
# flags and settings
#
HEADER = -I$(HEADER_PATH) -I$(CURDIR)
CFLAGS = -O3 -Wall
LIBS   = -L$(LIBS_PATH) \
		 $(lib)grampc$(a) \
		 $(lib)m


OBJS = $(SRCS:.c=.$(OBJEXT))


#
# compile fixed-size model
#
ifeq ($(FIXEDSIZE), 1)
	OUTPUT := $(OUTPUT)_fixedsize
	CFLAGS += -DFIXEDSIZE
	LIBS = -L./ \
		   $(lib)grampc_fixedsize \
		   $(lib)m
	export PROBLEM_PATH=$(CURDIR)
endif


#
# targets
#
all: grampc_all $(OUTPUT)

$(OUTPUT): $(OBJS)
	@echo 'Building target: $@'
	$(MKDIR) ./res
	$(COMPILER) -o "$@" $(OBJS) $(HEADER) $(LIBS)
	$(RM) $(OBJS)
	@echo 'Finished building: $@'
	@echo ''

$(OBJS): %.$(OBJEXT): %.c
	@echo 'Building file: $@'
	$(COMPILER) -o "$@" "$<" $(CFLAGS) $(HEADER) -c
	@echo ''

clean: grampc_clean
	@echo 'Cleaning up built files:'
	$(RM) $(OUTPUT) $(OBJS)
	@echo ''


#
# additional targets in fixed-size mode
#
grampc_all:
ifeq ($(FIXEDSIZE), 1)
	@echo 'Building GRAMPC'
	$(MAKE) -C $(PATH_PREFIX) all
	@echo ''
endif

grampc_clean:
ifeq ($(FIXEDSIZE), 1)
	@echo 'Cleaning GRAMPC'
	$(MAKE) -C $(PATH_PREFIX) clean
	@echo ''
endif

#
# end of make
#
