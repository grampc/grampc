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

#include "grampc_configreader.h"


/** Trim leading white-spaces and braces for vectors **/
void grampc_ltrim(typeChar *str, typeInt maxLineLength, typeInt type)
{
    typeInt index, i, j;
    index = 0;

    if (str[0] != '\n') {
        if (type == 0) { /* trim whitespaces for scalar fields */
            /* Find last index of whitespace */
            while( (str[index] == ' ' || str[index] == '\t' || str[index] == '\r' || str[index] == '\n' || str[index] == '\v' || str[index] == '\f') && index < maxLineLength - 1)
            {
                index++;
            }

            if(index != 0)
            {
                /* Shift all trailing characters to its left */
                i = 0;
                while (str[i + index] != '\0')
                {
                    str[i] = str[i + index];
                    i++;
                }
                j = i;
                while (j < i + index) {
                    str[j] = '\0'; // Make sure that string is NULL terminated
                    j++;
                }
            }
        }
        else { 
            /* Find last index of whitespace character and braces for vectors */
            while( (str[index] == ' ' || str[index] == '\t' || str[index] == '\r' || str[index] == '\n' || str[index] == '\v' || str[index] == '\f' || str[index] == '{' || str[index] == '[') && index < maxLineLength -1)
            {
                index++;
            }

            if(index != 0)
            {
                /* Shift all trailing characters to its left */
                i = 0;
                while(str[i + index] != '\0')
                {
                    str[i] = str[i + index];
                    i++;
                }
                j = i;
                while (j < i + index) {
                    str[j] = '\0'; // Make sure that string is NULL terminated
                    j++;
                }
            }
        }
    }
}


/** Trim trailing white-spaces and braces for vectors **/
void grampc_rtrim(typeChar *str, typeInt maxLineLength)
{
    typeInt index;
    index = maxLineLength-1;

    /* Find first index of whitespace or braces for vectors character from the back */
    while(str[index] == '\0' || str[index] == ' ' || str[index] == '\t' || str[index] == '\r' || str[index] == '\n' || str[index] == '\v' || str[index] == '\f'|| str[index] == '}' || str[index] == ']')
    {
        str[index] = '\0'; // Make sure that string is NULL terminated
        index--;
    }
}


/** Convert char pointer to array **/
void grampc_convert_to_vector(typeRNum *out, typeChar *configName, typeChar *configValue, typeInt size)
{
    typeChar *str = configValue;
    typeChar *end;

    typeInt i = 0;
    typeRNum f;

    if (out == NULL)
    {
        grampc_error_addstring("Can't allocate memory for string to vector conversion.\n", CONFIG_READER);
    }

    for (f = STRTOtypeRNum(str, &end); str != end; f = STRTOtypeRNum(str, &end))
    {
        /* Check for right vector sizes (input vec too long) */
        if (i == size) {
            grampc_error_addstring(INVALID_NO_ELEMENTS, configName);
        }

        /* Fill current array element */
        out[i] = f;
        i = i+1;

        /* Go to next element */
        str = end;

        /* Skip one index if vector is devided by commata or semicolons */
        if (str[0] == ',' || str[0] == ';')
            str = str + 1;
    }

    /* Check for right vector sizes (input vec too short) */
    if (i != size) {
        grampc_error_addstring(INVALID_NO_ELEMENTS, configName);
    }
}


/** Convert char pointer to integer array **/
void grampc_convert_to_int_vector(typeInt *out, typeChar *configName, typeChar *configValue, typeInt size)
{
    typeChar *str = configValue;
    typeChar *end;

    typeInt i = 0;
    typeInt f;

    if (out == NULL)
    {
        grampc_error_addstring("Can't allocate memory for string to vector conversion.\n", CONFIG_READER);
    }

    for (f = STRTOtypeInt(str, &end); str != end; f = STRTOtypeInt(str, &end))
    {
        /* Check for right vector sizes (input vec too long) */
        if (i == size) {
            grampc_error_addstring(INVALID_NO_ELEMENTS, configName);
        }

        /* Fill current array element */
        out[i] = f;
        i = i+1;

        /* Go to next element */
        str = end;

        /* Skip one index if vector is devided by commata or semicolons */
        if (str[0] == ',' || str[0] == ';')
            str = str + 1;
    }

    /* Check for right vector sizes (input vec too short) */
    if (i != size) {
        grampc_error_addstring(INVALID_NO_ELEMENTS, configName);
    }
}


/** Read configuration file and set options and parameter **/
void grampc_get_config_from_file(const typeGRAMPC *grampc, const typeChar *fileName)
{
    /* Open configuration file */
    FILE *configFile = fopen(fileName, "rt");

    if (configFile == NULL)
    {
        grampc_error_addstring("Can't open config file.\n", CONFIG_READER);
    }

    /* Allocate memory for reading from file */
    typeChar *line, *currentSectionName, *configName, *configValue;
    typeChar errorString[60];
    typeInt maxLineLength = MAX(grampc->param->Nx, MAX(grampc->param->Nu, MAX(grampc->param->Np, grampc->param->Nc))) * 10 + 100;
    line = (typeChar*)calloc(maxLineLength, sizeof(*line));
    currentSectionName = (typeChar*)calloc(maxLineLength, sizeof(*currentSectionName));
    configName = (typeChar*)calloc(maxLineLength, sizeof(*configName));
    configValue = (typeChar*)calloc(maxLineLength, sizeof(*configValue));

    if (line == NULL || currentSectionName == NULL || configName == NULL || configValue == NULL)
    {
        grampc_error_addstring("Can't allocate enough memory.\n", CONFIG_READER);
    }

    typeChar *tmp = NULL; // temporary variable needed for str to floating point conversion
    typeInt i;
    typeInt line_number = 0;

    /** Loop over all lines in the configuration file **/
    while (fgets(line, maxLineLength, configFile) != NULL) {

        /* Check if the line buffer can't read the whole line */
        line_number += 1;
        if (strlen(line) == maxLineLength - 1)
        {
            snprintf(errorString, 60, "Line length at line number %d too long.\n", line_number);
            grampc_error_addstring(errorString, CONFIG_READER);
        }
        
        /* Trim leading whitspaces */
        grampc_ltrim(line, maxLineLength, 0);

        /* Skip comments and empty lines */
        if (line[0] == '\n' || line[0] == ';' || line[0] == '#' || line[0] == 0) {
            /* Reset line pointer */
            memset(line, 0 , maxLineLength);
        }

        /* New section */
        else if (line[0] == '[') {
            /* Clear section name pointer */
            memset(currentSectionName, 0 , strlen(currentSectionName));

            /* Fill pointer elements with new section name */
            i = 1;
            while (i < maxLineLength - 1 && !(line[i] == ']')) {
                currentSectionName[i-1] = line[i];
                i++;
            }

            /* Trim whitespace */
            grampc_rtrim(currentSectionName, maxLineLength);
            grampc_ltrim(currentSectionName, maxLineLength, 0);

            /* Reset line pointer */
            memset(line, 0 , maxLineLength);
        }

        /* Not a comment/empty line or a section line, must be a name [=] value pair */
        else {
            /* Clear config name and value pointers */
            memset(configName, 0 , strlen(configName));
            memset(configValue, 0 , maxLineLength);

            /* Find position of = */
            typeChar *end = strchr(line, (typeInt)'=');
            if(end == NULL) {
                snprintf(errorString, 60, "Line %d contains no valid Name=Value pair.\n", line_number);
                grampc_error_addstring(errorString, CONFIG_READER);
            }
            typeInt pos = (typeInt)(end-line);

            /* Fill new name and value elements */
            for (i = 0; i < pos; i++)
            {
                configName[i] = line[i];
            }
            for (i = pos+1; i < maxLineLength; i++)
            {
                configValue[i-pos-1] = line[i];
            }

            /* Reset line pointer */
            memset(line, 0 , maxLineLength);

            /* Trim whitespace and braces for vectors */
            grampc_rtrim(configName, maxLineLength);
            grampc_rtrim(configValue, maxLineLength);
            grampc_ltrim(configValue, maxLineLength, 1);

            /** Set parameter and options **/
            /** MPC parameter **/
            if (!strcmp(currentSectionName, "GRAMPC parameter")) {

                /** MPC parameter of type real **/
                if (!strcmp(configName, "Thor") || !strcmp(configName, "Tmax") || !strcmp(configName, "Tmin") || !strcmp(configName, "dt") || !strcmp(configName, "t0")) {
                    /* Convert from char* to double and set values */
                    grampc_setparam_real(grampc, configName, STRTOtypeRNum(configValue, &tmp));
                }

                /** MPC parameter of type real vector **/
                else if (!strcmp(configName, "x0") || !strcmp(configName, "xdes") || !strcmp(configName, "u0") || !strcmp(configName, "udes") || !strcmp(configName, "umax") || !strcmp(configName, "umin") || !strcmp(configName, "p0") || !strcmp(configName, "pmax") || !strcmp(configName, "pmin")) {
                    /* Set size of parameter vector */
                    typeRNum *value;
                    if (!strcmp(configName, "x0") || !strcmp(configName, "xdes")) {
                        value = (typeRNum*)calloc(grampc->param->Nx, sizeof(*value));
                        /* Convert char pointer to vector */
                        grampc_convert_to_vector(value, configName, configValue, grampc->param->Nx);
                    }
                    else if (!strcmp(configName, "u0") || !strcmp(configName, "udes") || !strcmp(configName, "umax") || !strcmp(configName, "umin")) {
                        value = (typeRNum*)calloc(grampc->param->Nu, sizeof(*value));
                        /* Convert char pointer to vector */
                        grampc_convert_to_vector(value, configName, configValue, grampc->param->Nu);
                    }
                    else {
                        value = (typeRNum*)calloc(grampc->param->Np, sizeof(*value));
                        /* Convert char pointer to vector */
                        grampc_convert_to_vector(value, configName, configValue, grampc->param->Np);
                    }

                    /* Set parameter */
                    grampc_setparam_real_vector(grampc, configName, value);
                    free(value);
                }

                /** Undefined parameter name **/
                else {
                    grampc_error_addstring(INVALID_PARAM_NAME, configName);
                }
            }

            /** MPC options **/
            else if (!strcmp(currentSectionName, "GRAMPC options")) {

                /** MPC options of type real **/
                if (!strcmp(configName, "IntegratorRelTol") || !strcmp(configName, "IntegratorAbsTol") || !strcmp(configName, "IntegratorMinStepSize") || !strcmp(configName, "LineSearchMax") || !strcmp(configName, "LineSearchMin") || !strcmp(configName, "LineSearchInit") || !strcmp(configName, "LineSearchIntervalFactor") || !strcmp(configName, "LineSearchAdaptFactor") || !strcmp(configName, "LineSearchIntervalTol") || !strcmp(configName, "OptimParamLineSearchFactor") || !strcmp(configName, "OptimTimeLineSearchFactor") || !strcmp(configName, "TScale") || !strcmp(configName, "TOffset") || !strcmp(configName, "JScale") || !strcmp(configName, "MultiplierMax") || !strcmp(configName, "MultiplierDampingFactor") || !strcmp(configName, "PenaltyMax") || !strcmp(configName, "PenaltyMin") || !strcmp(configName, "PenaltyIncreaseFactor") || !strcmp(configName, "PenaltyDecreaseFactor") || !strcmp(configName, "PenaltyIncreaseThreshold") || !strcmp(configName, "AugLagUpdateGradientRelTol") || !strcmp(configName, "ConvergenceGradientRelTol")) {
                    /* Convert from char* to double and set values */
                    grampc_setopt_real(grampc, configName, STRTOtypeRNum(configValue, &tmp));
                }

                /** MPC options of type integer **/
                else if (!strcmp(configName, "MaxGradIter") || !strcmp(configName, "MaxMultIter") || !strcmp(configName, "Nhor") || !strcmp(configName, "IntegratorMaxSteps")) {
                    /* Convert from char* to integer and set values */
                    grampc_setopt_int(grampc, configName, atoi(configValue));
                }

                /** MPC options of type string **/
                else if (!strcmp(configName, "ShiftControl") || !strcmp(configName, "ScaleProblem") || !strcmp(configName, "TimeDiscretization") || !strcmp(configName, "IntegratorCost") || !strcmp(configName, "Integrator") || !strcmp(configName, "GradEvalType") || !strcmp(configName, "LineSearchType") || !strcmp(configName, "LineSearchExpAutoFallback") || !strcmp(configName, "OptimControl") || !strcmp(configName, "OptimParam") || !strcmp(configName, "OptimTime") || !strcmp(configName, "IntegralCost") || !strcmp(configName, "TerminalCost") || !strcmp(configName, "EqualityConstraints") || !strcmp(configName, "InequalityConstraints") || !strcmp(configName, "TerminalEqualityConstraints") || !strcmp(configName, "TerminalInequalityConstraints") || !strcmp(configName, "ConstraintsHandling") || !strcmp(configName, "ConvergenceCheck")) {
                    grampc_setopt_string(grampc, configName, configValue);
                }

                /** MPC options of type real vector **/
                else if (!strcmp(configName, "xScale") || !strcmp(configName, "xOffset") || !strcmp(configName, "uScale") || !strcmp(configName, "uOffset") || !strcmp(configName, "pScale") || !strcmp(configName, "pOffset") || !strcmp(configName, "cScale") || !strcmp(configName, "ConstraintsAbsTol")) {
                    /* Set size of parameter vector */
                    typeRNum *value;
                    if (!strcmp(configName, "xScale") || !strcmp(configName, "xOffset")) {
                        value = (typeRNum*)calloc(grampc->param->Nx, sizeof(*value));
                        /* Convert char pointer to vector */
                        grampc_convert_to_vector(value, configName, configValue,grampc->param->Nx);
                    }
                    else if (!strcmp(configName, "uScale") || !strcmp(configName, "uOffset")) {
                        value = (typeRNum*)calloc(grampc->param->Nu, sizeof(*value));
                        /* Convert char pointer to vector */
                        grampc_convert_to_vector(value, configName, configValue,grampc->param->Nu);
                    }
                    else if (!strcmp(configName, "cScale") || !strcmp(configName, "ConstraintsAbsTol")) {
                        value = (typeRNum*)calloc(grampc->param->Nc, sizeof(*value));
                        /* Convert char pointer to vector */
                        grampc_convert_to_vector(value, configName, configValue, grampc->param->Nc);
                    }
                    else {
                        value = (typeRNum*)calloc(grampc->param->Np, sizeof(*value));
                        /* Convert char pointer to vector */
                        grampc_convert_to_vector(value, configName, configValue, grampc->param->Np);
                    }

                    /* Set parameter */
                    grampc_setopt_real_vector(grampc, configName, value);
                    free(value);
                }

                /** MPC options of type integer vector **/
                else if (!strcmp(configName, "FlagsRodas")) {
                    /* Set size of parameter vector */
                    typeInt *value;
                    value = (typeInt*)calloc(8, sizeof(*value));

                    /* Convert char pointer to vector */
                    grampc_convert_to_int_vector(value, configName, configValue, 8);

                    /* Set parameter */
                    grampc_setopt_int_vector(grampc, configName, value);
                    free(value);
                }

                /** Undefined option name **/
                else {
                    grampc_error_addstring(INVALID_OPTION_NAME, configName);
                }
            }

            /** Undefined section name **/
            else {
                grampc_error_addstring(INVALID_SECTION_NAME, currentSectionName);
            }
        }
    }
    /** End loop over all lines in the configuration file **/

    /* Free char pointer memory */
    free(line);
    free(currentSectionName);
    free(configName);
    free(configValue);

    /* Close configuration file */
    fclose(configFile);
}
