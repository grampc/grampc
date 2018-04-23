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
SRCS_NLC_4 = main_NLChain_4.c
SRCS_NLC_6 = main_NLChain_6.c
SRCS_NLC_8 = main_NLChain_8.c
SRCS_NLC_10 = main_NLChain_10.c
SRCS_NLC_12 = main_NLChain_12.c 
SRCPROB_NLC_4 = probfct_NLChain_4.c	
SRCPROB_NLC_6 = probfct_NLChain_6.c
SRCPROB_NLC_8 = probfct_NLChain_8.c
SRCPROB_NLC_10 = probfct_NLChain_10.c
SRCPROB_NLC_12 = probfct_NLChain_12.c
		
OUTPUT_NLC_4 = startMPC_NLC_4
OUTPUT_NLC_6 = startMPC_NLC_6
OUTPUT_NLC_8 = startMPC_NLC_8
OUTPUT_NLC_10 = startMPC_NLC_10
OUTPUT_NLC_12 = startMPC_NLC_12

PATHRESNLC_4 = ./resNLC_4
PATHRESNLC_6 = ./resNLC_6
PATHRESNLC_8 = ./resNLC_8
PATHRESNLC_10 = ./resNLC_10
PATHRESNLC_12 = ./resNLC_12



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

OBJS_NLC_4 = $(SRCS_NLC_4:.c=.$(OBJEXT))
OBJS_NLC_6 = $(SRCS_NLC_6:.c=.$(OBJEXT))
OBJS_NLC_8 = $(SRCS_NLC_8:.c=.$(OBJEXT))
OBJS_NLC_10 = $(SRCS_NLC_10:.c=.$(OBJEXT))
OBJS_NLC_12 = $(SRCS_NLC_12:.c=.$(OBJEXT))

OBJPROB_NLC_4 = $(SRCPROB_NLC_4:.c=.$(OBJEXT))
OBJPROB_NLC_6 = $(SRCPROB_NLC_6:.c=.$(OBJEXT))
OBJPROB_NLC_8 = $(SRCPROB_NLC_8:.c=.$(OBJEXT))
OBJPROB_NLC_10 = $(SRCPROB_NLC_10:.c=.$(OBJEXT))
OBJPROB_NLC_12 = $(SRCPROB_NLC_12:.c=.$(OBJEXT))


#
# targets
#
all: $(OUTPUT_NLC_4)  $(OUTPUT_NLC_6) $(OUTPUT_NLC_8)  $(OUTPUT_NLC_10) $(OUTPUT_NLC_12)
NLC_4: $(OUTPUT_NLC_4)
NLC_6: $(OUTPUT_NLC_6)
NLC_8: $(OUTPUT_NLC_8)
NLC_10: $(OUTPUT_NLC_10)
NLC_12: $(OUTPUT_NLC_12)


$(OUTPUT_NLC_4): $(OBJS_NLC_4) $(OBJPROB_NLC_4)
	@echo 'Building target: $@'
	$(MKDIR) $(PATHRESNLC_4)
	$(COMPILER) -o "$@" $(OBJS_NLC_4) $(OBJPROB_NLC_4) $(HEADER) $(LIBS)
	@echo 'Finished building: $@'
	@echo '' 
	
$(OUTPUT_NLC_6): $(OBJS_NLC_6) $(OBJPROB_NLC_6)
	@echo 'Building target: $@'
	$(MKDIR) $(PATHRESNLC_6)
	$(COMPILER) -o "$@" $(OBJS_NLC_6) $(OBJPROB_NLC_6) $(HEADER) $(LIBS)
	@echo 'Finished building: $@'
	@echo '' 

$(OUTPUT_NLC_8): $(OBJS_NLC_8) $(OBJPROB_NLC_8)
	@echo 'Building target: $@'
	$(MKDIR) $(PATHRESNLC_8)
	$(COMPILER) -o "$@" $(OBJS_NLC_8) $(OBJPROB_NLC_8) $(HEADER) $(LIBS)
	@echo 'Finished building: $@'
	@echo '' 	
	
$(OUTPUT_NLC_10): $(OBJS_NLC_10) $(OBJPROB_NLC_10)
	@echo 'Building target: $@'
	$(MKDIR) $(PATHRESNLC_10)
	$(COMPILER) -o "$@" $(OBJS_NLC_10) $(OBJPROB_NLC_10) $(HEADER) $(LIBS)
	@echo 'Finished building: $@'
	@echo '' 	
	
$(OUTPUT_NLC_12): $(OBJS_NLC_12) $(OBJPROB_NLC_12)
	@echo 'Building target: $@'
	$(MKDIR) $(PATHRESNLC_12)
	$(COMPILER) -o "$@" $(OBJS_NLC_12) $(OBJPROB_NLC_12) $(HEADER) $(LIBS)
	@echo 'Finished building: $@'
	@echo '' 	
	
$(OBJS_NLC_4): %.$(OBJEXT): %.c
	@echo 'Building file: $@'
	$(COMPILER) -o "$@" "$<" $(CFLAGS) $(HEADER) -c
	@echo '' 

$(OBJS_NLC_6): %.$(OBJEXT): %.c
	@echo 'Building file: $@'
	$(COMPILER) -o "$@" "$<" $(CFLAGS) $(HEADER) -c
	@echo '' 

$(OBJS_NLC_8): %.$(OBJEXT): %.c
	@echo 'Building file: $@'
	$(COMPILER) -o "$@" "$<" $(CFLAGS) $(HEADER) -c
	@echo '' 
	
$(OBJS_NLC_10): %.$(OBJEXT): %.c
	@echo 'Building file: $@'
	$(COMPILER) -o "$@" "$<" $(CFLAGS) $(HEADER) -c
	@echo '' 
	
$(OBJS_NLC_12): %.$(OBJEXT): %.c
	@echo 'Building file: $@'
	$(COMPILER) -o "$@" "$<" $(CFLAGS) $(HEADER) -c
	@echo '' 
	
$(OBJPROB_NLC_4): %.$(OBJEXT): %.c
	@echo 'Building file: $@'
	$(COMPILER) -o "$@" "$<" $(CFLAGS) $(HEADER) -c
	@echo '' 
	
$(OBJPROB_NLC_6): %.$(OBJEXT): %.c
	@echo 'Building file: $@'
	$(COMPILER) -o "$@" "$<" $(CFLAGS) $(HEADER) -c
	@echo ''

$(OBJPROB_NLC_8): %.$(OBJEXT): %.c
	@echo 'Building file: $@'
	$(COMPILER) -o "$@" "$<" $(CFLAGS) $(HEADER) -c
	@echo ''

$(OBJPROB_NLC_10): %.$(OBJEXT): %.c
	@echo 'Building file: $@'
	$(COMPILER) -o "$@" "$<" $(CFLAGS) $(HEADER) -c
	@echo ''

$(OBJPROB_NLC_12): %.$(OBJEXT): %.c
	@echo 'Building file: $@'
	$(COMPILER) -o "$@" "$<" $(CFLAGS) $(HEADER) -c
	@echo ''	
	
clean:
	@echo 'Cleaning up built files:'
	$(RM) $(OUTPUT_NLC_4) $(OBJS_NLC_4)
	$(RM) $(OUTPUT_NLC_6) $(OBJS_NLC_6)
	$(RM) $(OUTPUT_NLC_8) $(OBJS_NLC_8)
	$(RM) $(OUTPUT_NLC_10) $(OBJS_NLC_10)
	$(RM) $(OUTPUT_NLC_12) $(OBJS_NLC_12)
	
	$(RM) $(OBJPROB_NLC_4)
	$(RM) $(OBJPROB_NLC_6)
	$(RM) $(OBJPROB_NLC_8)
	$(RM) $(OBJPROB_NLC_10)
	$(RM) $(OBJPROB_NLC_12)
	
	$(RM) $(PATHRESNLC_4)/*.txt
	$(RM) $(PATHRESNLC_6)/*.txt
	$(RM) $(PATHRESNLC_8)/*.txt
	$(RM) $(PATHRESNLC_10)/*.txt
	$(RM) $(PATHRESNLC_12)/*.txt
	@echo ''

#
# end of make
#
