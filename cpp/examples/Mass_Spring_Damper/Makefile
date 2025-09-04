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
# input and output files
#
SRCS = main_MSD.cpp problem_description_MSD.cpp
OUTPUT = startMPC


#
# compiler
#
CPP_COMPILER = g++
RM           = rm -f
MKDIR        = mkdir -p
OBJEXT       = o
LIBEXT       = a
lib          = -l
PATH_PREFIX  = ../..


#
# path of header and lib files
#
HEADER_PATH     = $(PATH_PREFIX)/../include
CPP_HEADER_PATH = $(PATH_PREFIX)/include
CPP_LIBS_PATH   = $(PATH_PREFIX)/libs


#
# flags and settings
#
CPP_HEADER = -I${HEADER_PATH} -I${CPP_HEADER_PATH}
CPP_CFLAGS = -c -O3 -Wall
CPP_LIBS   = -L$(CPP_LIBS_PATH) \
		     $(lib)grampc$(a) \
		     $(lib)m

CPP_OBJS = $(SRCS:.cpp=.$(OBJEXT))

#
# targets
#
all: $(OUTPUT)

$(OUTPUT): $(CPP_OBJS)
	@echo 'Building target: $@'
#	$(MKDIR) ./res
	$(CPP_COMPILER) -o "$@" $(CPP_OBJS) $(CPP_HEADER) $(CPP_LIBS)
	$(RM) $(CPP_OBJS)
	@echo 'Finished building: $@'
	@echo ''

$(CPP_OBJS): %.$(OBJEXT): %.cpp
	@echo 'Building file: $@'
	$(CPP_COMPILER) -o "$@" "$<" $(CPP_CFLAGS) $(CPP_HEADER)
	@echo ''

clean:
	@echo 'Cleaning up built files:'
	$(RM) $(OUTPUT) $(CPP_OBJS)
	@echo ''

#
# end of make
#
