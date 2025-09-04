/* This file is part of GRAMPC - (https://github.com/grampc/grampc)
 *
 * GRAMPC -- A software framework for embedded nonlinear model predictive
 * control using a gradient-based augmented Lagrangian approach
 *
 * Copyright 2014-2025 by Knut Graichen, Andreas Voelz, Thore Wietzke,
 * Tobias Englert (<v2.3), Felix Mesmer (<v2.3), Soenke Rhein (<v2.3),
 * Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
 * All rights reserved.
 *
 * GRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
 *
 */


#ifndef GRAMPC_CONFIGREADER_H_
#define GRAMPC_CONFIGREADER_H_

#include "grampc_init.h"
#include "grampc_mess.h"
#include "grampc_util.h"
#include "grampc_setopt.h"
#include "grampc_setparam.h"

#include <stdlib.h>
#include <stdio.h>

#define CONFIG_READER "Config Reader"

void grampc_ltrim(typeChar *str, typeInt maxLineLength, typeInt type);
void grampc_rtrim(typeChar *str, typeInt maxLineLength);
void grampc_convert_to_vector(typeRNum *out, typeChar *configName, typeChar *configValue, typeInt size);

void grampc_get_config_from_file(const typeGRAMPC *grampc, const typeChar *fileName);

#if USE_typeRNum == USE_FLOAT
#define STRTOtypeRNum(str, end) strtof(str, end)
#else
#define STRTOtypeRNum(str, end) strtod(str, end)
#endif

#define STRTOtypeInt(str, end) (typeInt)strtol(str, end, 10)

#endif /* GRAMPC_CONFIGREADER_H_ */
