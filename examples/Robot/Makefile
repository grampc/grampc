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
# input and output files
#
SRCSOCP = main_OCP_ROBOT.c 
SRCPROB = probfct_ROBOT.c	
		
OUTPUTOCP = startOCP

PATHRESOCP = ./resOCP


#
# compiler
#
COMPILER    = gcc
RM          = rm -f
MKDIR       = mkdir -p
OBJEXT      = o
LIBEXT      = a
lib         = -l
a           = 
PATH_PREFIX = ../..

#
# path of header and libs files
#
HEADER_PATH = $(PATH_PREFIX)/include
LIBS_PATH   = $(PATH_PREFIX)/libs


#
# flags and settings
#
HEADER = -I$(HEADER_PATH)
CFLAGS = -O3 -Wall
LIBS   = -L$(LIBS_PATH) \
		 $(lib)m \
	 	 $(lib)grampc$(a)

OBJSOCP = $(SRCSOCP:.c=.$(OBJEXT))
OBJPROB = $(SRCPROB:.c=.$(OBJEXT))


#
# targets
#
all: $(OUTPUTOCP)
OCP: $(OUTPUTOCP)

$(OUTPUTOCP): $(OBJSOCP) $(OBJPROB)
	@echo 'Building target: $@'
	$(MKDIR) $(PATHRESOCP)
	$(COMPILER) -o "$@" $(OBJSOCP) $(OBJPROB) $(HEADER) $(LIBS)
	@echo 'Finished building: $@'
	@echo '' 	

$(OBJSOCP): %.$(OBJEXT): %.c
	@echo 'Building file: $@'
	$(COMPILER) -o "$@" "$<" $(CFLAGS) $(HEADER) -c
	@echo '' 

$(OBJPROB): %.$(OBJEXT): %.c
	@echo 'Building file: $@'
	$(COMPILER) -o "$@" "$<" $(CFLAGS) $(HEADER) -c
	@echo '' 
	
clean:
	@echo 'Cleaning up built files:'
	$(RM) $(OUTPUTOCP) $(OBJSOCP)
	$(RM) $(OBJPROB)
	$(RM) $(PATHRESOCP)/*.txt
	@echo ''

#
# end of make
#
