/*
 *
 * This file is part of GRAMPC.
 *
 * GRAMPC - a gradient-based MPC software for real-time applications
 *
 * Copyright (C) 2014 by Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Developed at the Institute of Measurement, Control, and
 * Microtechnology, University of Ulm. All rights reserved.
 *
 * GRAMPC is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as 
 * published by the Free Software Foundation, either version 3 of 
 * the License, or (at your option) any later version.
 *
 * GRAMPC is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public 
 * License along with GRAMPC. If not, see <http://www.gnu.org/licenses/>.
 *
 */


/*
 *
 * File: grampc_mess.c
 * Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Date: February 2014
 * Version: v1.0
 *
 * Error printing file for GRAMPC.
 *
 */


#include "grampc_mess.h"

void grampc_error(char *mess)
{
  printError(mess);
}

void grampc_error_addstring(char *mess, char *addstring)
{
  printErrorAddString(mess,addstring);
}

void grampc_error_optname(char *optname)
{
  printErrorAddString(INVALID_OPTION_NAME,optname);
}

void grampc_error_optvalue(char *optname)
{
  printErrorAddString(INVALID_OPTION_VALUE,optname);
}

void grampc_error_paramname(char *paramname)
{
  printErrorAddString(INVALID_PARAM_NAME,paramname);
}

void grampc_error_paramvalue(char *paramname)
{
  printErrorAddString(INVALID_PARAM_VALUE,paramname);
}