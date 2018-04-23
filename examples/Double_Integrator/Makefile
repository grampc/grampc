#
# This file is part of GRAMPC - (https://sourceforge.net/projects/grampc/)
#
# GRAMPC -- A software framework for embedded nonlinear model predictive
# control using a gradient-based augmented Lagrangian approach
#
# Copyright (C) 2014-2018 by Tobias Englert, Knut Graichen, Felix Mesmer,
# Soenke Rhein, Andreas Voelz, Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
# Developed at the Institute of Measurement, Control, and Microtechnology,
# Ulm University. All rights reserved.
#
# GRAMPC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as 
# published by the Free Software Foundation, either version 3 of 
# the License, or (at your option) any later version.
#
# GRAMPC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public 
# License along with GRAMPC. If not, see <http://www.gnu.org/licenses/>.
#


#
# input and output files
#
SRCSMPC = main_DOUBLE_INTEGRATOR.c
SRCSOCP = main_OCP_DOUBLE_INTEGRATOR.c 
SRCPROB = probfct_DOUBLE_INTEGRATOR.c	
SRCPROBOCP = probfct_DOUBLE_INTEGRATOR_OCP.c
		
OUTPUTMPC = startMPC
OUTPUTOCP = startOCP

PATHRESMPC = ./res
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

OBJSMPC = $(SRCSMPC:.c=.$(OBJEXT))
OBJSOCP = $(SRCSOCP:.c=.$(OBJEXT))
OBJPROB = $(SRCPROB:.c=.$(OBJEXT))
OBJPROBOCP = $(SRCPROBOCP:.c=.$(OBJEXT))


#
# targets
#
all: $(OUTPUTMPC)  $(OUTPUTOCP)
MPC: $(OUTPUTMPC)
OCP: $(OUTPUTOCP)

$(OUTPUTMPC): $(OBJSMPC) $(OBJPROB)
	@echo 'Building target: $@'
	$(MKDIR) $(PATHRESMPC)
	$(COMPILER) -o "$@" $(OBJSMPC) $(OBJPROB) $(HEADER) $(LIBS)
	@echo 'Finished building: $@'
	@echo '' 

$(OUTPUTOCP): $(OBJSOCP) $(OBJPROBOCP)
	@echo 'Building target: $@'
	$(MKDIR) $(PATHRESOCP)
	$(COMPILER) -o "$@" $(OBJSOCP) $(OBJPROBOCP) $(HEADER) $(LIBS)
	@echo 'Finished building: $@'
	@echo '' 	
	
$(OBJSMPC): %.$(OBJEXT): %.c
	@echo 'Building file: $@'
	$(COMPILER) -o "$@" "$<" $(CFLAGS) $(HEADER) -c
	@echo '' 

$(OBJSOCP): %.$(OBJEXT): %.c
	@echo 'Building file: $@'
	$(COMPILER) -o "$@" "$<" $(CFLAGS) $(HEADER) -c
	@echo '' 

$(OBJPROB): %.$(OBJEXT): %.c
	@echo 'Building file: $@'
	$(COMPILER) -o "$@" "$<" $(CFLAGS) $(HEADER) -c
	@echo '' 

$(OBJPROBOCP): %.$(OBJEXT): %.c
	@echo 'Building file: $@'
	$(COMPILER) -o "$@" "$<" $(CFLAGS) $(HEADER) -c
	@echo ''
	
clean:
	@echo 'Cleaning up built files:'
	$(RM) $(OUTPUTMPC) $(OBJSMPC)
	$(RM) $(OUTPUTOCP) $(OBJSOCP)
	$(RM) $(OBJPROB)
	$(RM) $(OBJPROBOCP)
	$(RM) $(PATHRESMPC)/*.txt
	$(RM) $(PATHRESOCP)/*.txt
	@echo ''

#
# end of make
#
