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
 * File: main.c
 * Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
 * Date: February 2014
 * Version: v1.0
 *
 * Mainfile for example Crane_2D.
 *
 */


#include "grampc.h"
#include <time.h>

#define NX	6
#define NU	2

int main()
{
  typeGRAMPC *grampc;
  typeInt iMPC, MaxIter, i;
  typeRNum Tsim  = 8.0;
  typeRNum t     = 0.0;
  typeRNum Thor  = 1.5;
  typeRNum dt    = 0.001;
  typeInt Nhor   = 40;
  typeInt NpCost = 14;

  typeRNum x0[NX]   = {-2.0,0.0,2.0,0.0,0.0,0.0};
  typeRNum u0[NU]   = {0.0,0.0};
  typeRNum xdes[NX] = {2.0,0.0,0.4,0.0,0.0,0.0};
  typeRNum udes[NU] = {0.0,0.0};

  typeRNum xScale[NX]  = {1.0,1.0,1.0,1.0,1.0,1.0};
  typeRNum xOffset[NX] = {0.0,0.0,0.0,0.0,0.0,0.0};
  typeRNum uScale[NU]  = {1.0,1.0};
  typeRNum uOffset[NU] = {0.0,0.0};

  typeRNum umax[NU] = {2.0,2.0};
  typeRNum umin[NU] = {-2.0,-2.0};

  typeRNum pCost[14] = {1.0,1.0,1.0,1.0,1.0,1.0,
			1.0,1.0,1.0,1.0,1.0,1.0,
			0.01,0.01};

  FILE *file_x;
  FILE *file_u;
  FILE *file_J;
  FILE *file_t;

  clock_t tic, toc;
  typeRNum *CPUtimeVec;
  typeRNum CPUtime = 0;

  file_x = fopen("xvec.txt","w");
  file_u = fopen("uvec.txt","w");
  file_J = fopen("Jvec.txt","w");
  file_t = fopen("tvec.txt","w");

  grampc_init(&grampc);
  /* grampc_setopt_int(grampc,"MaxIter",2); */
  /* grampc_setopt_string(grampc,"ShiftControl","off"); */
  /* grampc_setopt_string(grampc,"ScaleProblem","on"); */
  /* grampc_setopt_string(grampc,"CostIntegrator","simpson"); */
  grampc_setopt_string(grampc,"Integrator","euler");
  /* grampc_setopt_real(grampc,"IntegratorRelTol",1e-3); */
  /* grampc_setopt_real(grampc,"IntegratorAbsTol",1e-1); */
  grampc_setopt_string(grampc,"LineSearchType","explicit2");
  /* grampc_setopt_real(grampc,"LineSearchMax",2.0); */
  /* grampc_setopt_real(grampc,"LineSearchMin",1e-3); */
  /* grampc_setopt_real(grampc,"LineSearchInit",0.1); */
  /* grampc_setopt_real(grampc,"LineSearchIntervalFactor",0.5); */
  /* grampc_setopt_real(grampc,"LineSearchAdaptFactor",2.0); */
  /* grampc_setopt_real(grampc,"LineSearchIntervalTol",1e-3); */
  /* grampc_setopt_string(grampc,"FinalCost","on"); */
  /* grampc_setopt_string(grampc,"IntegralCost","on"); */
  /* grampc_setopt_string(grampc,"JacobianX","sysjacx"); */
  /* grampc_setopt_string(grampc,"JacobianU","sysjacu"); */
  grampc_printopt(grampc);

  grampc_setparam_real(grampc,"Thor",Thor);
  grampc_setparam_real(grampc,"dt",dt);
  /* grampc_setparam_vector(grampc,"xScale",xScale); */
  /* grampc_setparam_vector(grampc,"xOffset",xOffset); */
  /* grampc_setparam_vector(grampc,"uScale",uScale); */
  /* grampc_setparam_vector(grampc,"uOffset",uOffset); */
  grampc_setparam_int(grampc,"Nhor",Nhor);
  grampc_setparam_vector(grampc,"umax",umax);
  grampc_setparam_vector(grampc,"umin",umin);
  grampc_setparam_int(grampc,"NpCost",NpCost);
  grampc_setparam_vector(grampc,"pCost",pCost);
  grampc_setparam_vector(grampc,"xk",x0);
  grampc_setparam_vector(grampc,"u0",u0);
  grampc_setparam_vector(grampc,"xdes",xdes);
  grampc_setparam_vector(grampc,"udes",udes);
  grampc_printparam(grampc);

  MaxIter = (int)(Tsim/grampc->param->dt);

  CPUtimeVec = calloc(MaxIter+1,sizeof(*CPUtimeVec));

  printf("MPC running ...\n");
    for (iMPC = 1; iMPC <= MaxIter+1; iMPC++) {
    t = t + grampc->param->dt;
    tic = clock();
    grampc_run(grampc);
    grampc_setparam_vector(grampc,"xk",grampc->sol->xnext);
    toc = clock();
    CPUtimeVec[iMPC-1] = (typeRNum)((toc-tic)*1000/CLOCKS_PER_SEC);
    for (i = 0; i <= grampc->param->Nx-1; i++) {
      if (i <= grampc->param->Nx-2) {
        fprintf(file_x,"%.5f ,",grampc->sol->xnext[i]);
      }
    }
    for (i = 0; i <= grampc->param->Nu-2; i++) {
      fprintf(file_u,"%.5f ,",grampc->sol->unext[i]);
    }
    fprintf(file_x,"%.5f;\n",grampc->sol->xnext[grampc->param->Nx-1]);
    fprintf(file_u,"%.5f;\n",grampc->sol->unext[grampc->param->Nu-1]);
    fprintf(file_J,"%.5f;\n",grampc->sol->J[0]);
    fprintf(file_t,"%.5f;\n",t);
  }

  fclose(file_x);
  fclose(file_u);
  fclose(file_J);
  fclose(file_t);

  for (i = 0; i <= MaxIter; i++) {
    CPUtime = CPUtime + CPUtimeVec[i]/(MaxIter+1);
  }

  grampc_free(&grampc);
  free(CPUtimeVec);

  printf("MPC finished. Average computation time: %.3f ms.\n",CPUtime);

  return 0;
}
