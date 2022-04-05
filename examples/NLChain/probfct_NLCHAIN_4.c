/* This file is part of GRAMPC - (https://sourceforge.net/projects/grampc/)
 *
 * GRAMPC -- A software framework for embedded nonlinear model predictive
 * control using a gradient-based augmented Lagrangian approach
 *
 * Copyright 2014-2019 by Tobias Englert, Knut Graichen, Felix Mesmer,
 * Soenke Rhein, Andreas Voelz, Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
 * All rights reserved.
 *
 * GRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
 *
 *
 *
 *
 *
 *
 *
 *
 * This probfct-file describes the nonlinear chain problem with four masses from
 * Kirches, C., Wirsching, L., Bock, H., Schloder, J.: Efficient Direct multiple
 * Shooting for Nonlinear Model Predictive Control on Long Horizons. Journal of
 * Process Control 22(3), 540-550 (2012)
 *
 *                                           _T
 *		                            /
 *      min    J(u,p,T;x0) = V(T,x(T),p) + / l(t,x(t),u(t),p) dt
 *   u(.),p,T                            _/
 *                                      0
 *             .
 *      s.t.   x(t) = f(t0+t,x(t),u(t),p), x(0) = x0
 *             u_min <= u(t) <= u_max
 *
 */

#include "probfct.h"

#define POW2(a) ((a)*(a))

 /** OCP dimensions: states (Nx), controls (Nu), parameters (Np), equalities (Ng), 
    inequalities (Nh), terminal equalities (NgT), terminal inequalities (NhT) **/
void ocp_dim(typeInt *Nx, typeInt *Nu, typeInt *Np, typeInt *Ng, typeInt *Nh, typeInt *NgT, typeInt *NhT, typeUSERPARAM *userparam)
{
	*Nx = 21;
	*Nu = 3;
	*Np = 0;
	*Nh = 0;
	*Ng = 0;
	*NgT = 0;
	*NhT = 0;
}

/** System function f(t,x,u,p,userparam) 
    ------------------------------------ **/
void ffct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	typeRNum yaux100, yaux101, yaux102, yaux103, yaux104, yaux105, yaux106, yaux107, yaux13, yaux14;
	typeRNum yaux15, yaux16, yaux17, yaux18, yaux19, yaux20, yaux21, yaux22, yaux24, yaux25;
	typeRNum yaux26, yaux27, yaux28, yaux29, yaux30, yaux31, yaux32, yaux33, yaux34, yaux35;
	typeRNum yaux36, yaux37, yaux38, yaux39, yaux40, yaux41, yaux46, yaux47, yaux52, yaux53;
	typeRNum yaux59, yaux60, yaux61, yaux62, yaux63, yaux64, yaux65, yaux66, yaux67, yaux68;
	typeRNum yaux69, yaux70, yaux71, yaux72, yaux73, yaux74, yaux75, yaux80, yaux85, yaux91;
	typeRNum yaux93, yaux94, yaux95, yaux96, yaux97, yaux98, yaux99;

	yaux13 = x[0];
	yaux25 = x[3];
	yaux15 = x[1];
	yaux17 = x[2];
	yaux14 = yaux13 * yaux13;
	yaux16 = yaux15 * yaux15;
	yaux18 = yaux17 * yaux17;
	yaux19 = yaux14 + yaux16 + yaux18;
	yaux20 = 1 / SQRT(yaux19);
	yaux21 = -0.55*yaux20;
	yaux22 = 4 + yaux21;
	yaux30 = x[4];
	yaux27 = -yaux25;
	yaux28 = yaux13 + yaux27;
	yaux29 = yaux28 * yaux28;
	yaux31 = -yaux30;
	yaux32 = yaux15 + yaux31;
	yaux33 = yaux32 * yaux32;
	yaux34 = x[5];
	yaux35 = -yaux34;
	yaux36 = yaux17 + yaux35;
	yaux37 = yaux36 * yaux36;
	yaux38 = yaux29 + yaux33 + yaux37;
	yaux39 = 1 / SQRT(yaux38);
	yaux40 = -0.55*yaux39;
	yaux41 = 4 + yaux40;
	yaux24 = -yaux13;
	yaux26 = yaux24 + yaux25;
	yaux59 = x[6];
	yaux46 = -yaux15;
	yaux47 = yaux30 + yaux46;
	yaux64 = x[7];
	yaux61 = -yaux59;
	yaux62 = yaux25 + yaux61;
	yaux63 = yaux62 * yaux62;
	yaux65 = -yaux64;
	yaux66 = yaux30 + yaux65;
	yaux67 = yaux66 * yaux66;
	yaux68 = x[8];
	yaux69 = -yaux68;
	yaux70 = yaux34 + yaux69;
	yaux71 = yaux70 * yaux70;
	yaux72 = yaux63 + yaux67 + yaux71;
	yaux73 = 1 / SQRT(yaux72);
	yaux74 = -0.55*yaux73;
	yaux75 = 4 + yaux74;
	yaux52 = -yaux17;
	yaux53 = yaux34 + yaux52;
	yaux60 = yaux27 + yaux59;
	yaux91 = x[9];
	yaux80 = yaux31 + yaux64;
	yaux96 = x[10];
	yaux93 = -yaux91;
	yaux94 = yaux59 + yaux93;
	yaux95 = yaux94 * yaux94;
	yaux97 = -yaux96;
	yaux98 = yaux64 + yaux97;
	yaux99 = yaux98 * yaux98;
	yaux100 = x[11];
	yaux101 = -yaux100;
	yaux102 = yaux101 + yaux68;
	yaux103 = yaux102 * yaux102;
	yaux104 = yaux103 + yaux95 + yaux99;
	yaux105 = 1 / SQRT(yaux104);
	yaux106 = -0.55*yaux105;
	yaux107 = 4 + yaux106;
	yaux85 = yaux35 + yaux68;

	out[0] = x[12];
	out[1] = x[13];
	out[2] = x[14];
	out[3] = x[15];
	out[4] = x[16];
	out[5] = x[17];
	out[6] = x[18];
	out[7] = x[19];
	out[8] = x[20];
	out[9] = u[0];
	out[10] = u[1];
	out[11] = u[2];
	out[12] = 8.88888888888889*(-0.1*yaux13*yaux22 + 0.1*yaux26*yaux41);
	out[13] = 8.88888888888889*(-0.1*yaux15*yaux22 + 0.1*yaux41*yaux47);
	out[14] = -9.81 + 8.88888888888889*(-0.1*yaux17*yaux22 + 0.1*yaux41*yaux53);
	out[15] = 8.88888888888889*(-0.1*yaux26*yaux41 + 0.1*yaux60*yaux75);
	out[16] = 8.88888888888889*(-0.1*yaux41*yaux47 + 0.1*yaux75*yaux80);
	out[17] = -9.81 + 8.88888888888889*(-0.1*yaux41*yaux53 + 0.1*yaux75*yaux85);
	out[18] = 8.88888888888889*(-0.1*yaux60*yaux75 + 0.1*yaux107*(yaux61 + yaux91));
	out[19] = 8.88888888888889*(-0.1*yaux75*yaux80 + 0.1*yaux107*(yaux65 + yaux96));
	out[20] = -9.81 + 8.88888888888889*(0.1*yaux107*(yaux100 + yaux69) - 0.1*yaux75*yaux85);
}
/** Jacobian df/dx multiplied by vector vec, i.e. (df/dx)^T*vec or vec^T*(df/dx) **/
void dfdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	typeRNum yaux1, yaux10, yaux100, yaux101, yaux102, yaux103, yaux104, yaux105, yaux106, yaux107;
	typeRNum yaux108, yaux109, yaux11, yaux110, yaux111, yaux112, yaux113, yaux114, yaux115, yaux116;
	typeRNum yaux117, yaux119, yaux12, yaux120, yaux123, yaux124, yaux127, yaux129, yaux13, yaux130;
	typeRNum yaux132, yaux139, yaux14, yaux142, yaux143, yaux145, yaux15, yaux150, yaux157, yaux16;
	typeRNum yaux160, yaux165, yaux166, yaux168, yaux17, yaux175, yaux176, yaux177, yaux178, yaux179;
	typeRNum yaux18, yaux180, yaux181, yaux182, yaux183, yaux184, yaux185, yaux186, yaux187, yaux188;
	typeRNum yaux189, yaux19, yaux190, yaux191, yaux192, yaux193, yaux197, yaux2, yaux20, yaux201;
	typeRNum yaux21, yaux22, yaux23, yaux234, yaux235, yaux236, yaux237, yaux24, yaux25, yaux26;
	typeRNum yaux28, yaux29, yaux3, yaux30, yaux31, yaux32, yaux33, yaux35, yaux36, yaux37;
	typeRNum yaux38, yaux4, yaux41, yaux42, yaux43, yaux44, yaux45, yaux48, yaux5, yaux50;
	typeRNum yaux51, yaux52, yaux54, yaux55, yaux56, yaux6, yaux60, yaux63, yaux64, yaux67;
	typeRNum yaux7, yaux72, yaux73, yaux77, yaux8, yaux80, yaux85, yaux86, yaux89, yaux9;
	typeRNum yaux96, yaux97, yaux98, yaux99;

	yaux2 = x[0];
	yaux3 = x[3];
	yaux4 = -yaux3;
	yaux5 = yaux2 + yaux4;
	yaux6 = yaux5 * yaux5;
	yaux7 = x[1];
	yaux8 = x[4];
	yaux9 = -yaux8;
	yaux10 = yaux7 + yaux9;
	yaux11 = yaux10 * yaux10;
	yaux12 = x[2];
	yaux13 = x[5];
	yaux14 = -yaux13;
	yaux15 = yaux12 + yaux14;
	yaux16 = yaux15 * yaux15;
	yaux17 = yaux11 + yaux16 + yaux6;
	yaux29 = yaux2 * yaux2;
	yaux30 = yaux7 * yaux7;
	yaux31 = yaux12 * yaux12;
	yaux32 = yaux29 + yaux30 + yaux31;
	yaux18 = 1 / SQRT(yaux17);
	yaux19 = -0.55*yaux18;
	yaux20 = 4 + yaux19;
	yaux22 = -yaux2;
	yaux23 = yaux22 + yaux3;
	yaux24 = POW(yaux17, -1.5);
	yaux33 = POW(yaux32, -1.5);
	yaux43 = -yaux7;
	yaux44 = yaux43 + yaux8;
	yaux51 = -yaux12;
	yaux52 = yaux13 + yaux51;
	yaux28 = vec[12];
	yaux42 = -0.05500000000000001*yaux2*yaux33*yaux7;
	yaux48 = vec[16];
	yaux21 = 0.1*yaux20;
	yaux41 = vec[13];
	yaux35 = 1 / SQRT(yaux32);
	yaux36 = 0.05500000000000001*yaux35;
	yaux37 = -0.1*yaux20;
	yaux1 = vec[15];
	yaux50 = vec[17];
	yaux54 = vec[14];
	yaux55 = -0.05500000000000001*yaux12*yaux2*yaux33;
	yaux72 = -0.05500000000000001*yaux12*yaux33*yaux7;
	yaux25 = -0.05500000000000001*yaux23*yaux24*yaux5;
	yaux26 = yaux21 + yaux25;
	yaux97 = x[6];
	yaux98 = -yaux97;
	yaux99 = yaux3 + yaux98;
	yaux100 = yaux99 * yaux99;
	yaux101 = x[7];
	yaux102 = -yaux101;
	yaux103 = yaux102 + yaux8;
	yaux104 = yaux103 * yaux103;
	yaux105 = x[8];
	yaux106 = -yaux105;
	yaux107 = yaux106 + yaux13;
	yaux108 = yaux107 * yaux107;
	yaux109 = yaux100 + yaux104 + yaux108;
	yaux38 = 0.05500000000000001*yaux23*yaux24*yaux5;
	yaux110 = 1 / SQRT(yaux109);
	yaux111 = -0.55*yaux110;
	yaux112 = 4 + yaux111;
	yaux114 = yaux4 + yaux97;
	yaux115 = POW(yaux109, -1.5);
	yaux45 = 0.05500000000000001*yaux24*yaux44*yaux5;
	yaux123 = yaux101 + yaux9;
	yaux56 = 0.05500000000000001*yaux24*yaux5*yaux52;
	yaux130 = yaux105 + yaux14;
	yaux63 = -0.05500000000000001*yaux10*yaux24*yaux44;
	yaux64 = yaux21 + yaux63;
	yaux60 = 0.05500000000000001*yaux10*yaux23*yaux24;
	yaux127 = vec[19];
	yaux113 = 0.1*yaux112;
	yaux67 = 0.05500000000000001*yaux10*yaux24*yaux44;
	yaux119 = -0.1*yaux112;
	yaux96 = vec[18];
	yaux129 = vec[20];
	yaux73 = 0.05500000000000001*yaux10*yaux24*yaux52;
	yaux85 = -0.05500000000000001*yaux15*yaux24*yaux52;
	yaux86 = yaux21 + yaux85;
	yaux77 = 0.05500000000000001*yaux15*yaux23*yaux24;
	yaux80 = 0.05500000000000001*yaux15*yaux24*yaux44;
	yaux89 = 0.05500000000000001*yaux15*yaux24*yaux52;
	yaux116 = -0.05500000000000001*yaux114*yaux115*yaux99;
	yaux117 = yaux113 + yaux116;
	yaux120 = 0.05500000000000001*yaux114*yaux115*yaux99;
	yaux175 = x[9];
	yaux176 = -yaux175;
	yaux177 = yaux176 + yaux97;
	yaux178 = yaux177 * yaux177;
	yaux179 = x[10];
	yaux180 = -yaux179;
	yaux181 = yaux101 + yaux180;
	yaux182 = yaux181 * yaux181;
	yaux183 = x[11];
	yaux184 = -yaux183;
	yaux185 = yaux105 + yaux184;
	yaux186 = yaux185 * yaux185;
	yaux187 = yaux178 + yaux182 + yaux186;
	yaux124 = 0.05500000000000001*yaux115*yaux123*yaux99;
	yaux193 = POW(yaux187, -1.5);
	yaux132 = 0.05500000000000001*yaux115*yaux130*yaux99;
	yaux142 = -0.05500000000000001*yaux103*yaux115*yaux123;
	yaux143 = yaux113 + yaux142;
	yaux139 = 0.05500000000000001*yaux103*yaux114*yaux115;
	yaux192 = yaux175 + yaux98;
	yaux145 = 0.05500000000000001*yaux103*yaux115*yaux123;
	yaux188 = 1 / SQRT(yaux187);
	yaux189 = -0.55*yaux188;
	yaux190 = 4 + yaux189;
	yaux191 = -0.1*yaux190;
	yaux197 = yaux102 + yaux179;
	yaux150 = 0.05500000000000001*yaux103*yaux115*yaux130;
	yaux201 = yaux106 + yaux183;
	yaux165 = -0.05500000000000001*yaux107*yaux115*yaux130;
	yaux166 = yaux113 + yaux165;
	yaux157 = 0.05500000000000001*yaux107*yaux114*yaux115;
	yaux160 = 0.05500000000000001*yaux107*yaux115*yaux123;
	yaux168 = 0.05500000000000001*yaux107*yaux115*yaux130;
	yaux234 = SQRT(yaux187);
	yaux235 = -yaux234;
	yaux236 = 0.1375 + yaux235;
	yaux237 = -7.2727272727272725*yaux187*yaux236;

	out[0] = 8.88888888888889*yaux1*yaux26 + 8.88888888888889*yaux28*(-0.4 - 0.05500000000000001*yaux29*yaux33 + yaux36 + yaux37 + yaux38) + 8.88888888888889*yaux41*(yaux42 + yaux45) - 0.488888888888889*yaux24*yaux44*yaux48*yaux5 - 0.488888888888889*yaux24*yaux5*yaux50*yaux52 + 8.88888888888889*yaux54*(yaux55 + yaux56);
	out[1] = -0.488888888888889*yaux1*yaux10*yaux23*yaux24 - 0.488888888888889*yaux10*yaux24*yaux50*yaux52 + 8.88888888888889*yaux28*(yaux42 + yaux60) + 8.88888888888889*yaux48*yaux64 + 8.88888888888889*yaux41*(-0.4 - 0.05500000000000001*yaux30*yaux33 + yaux36 + yaux37 + yaux67) + 8.88888888888889*yaux54*(yaux72 + yaux73);
	out[2] = -0.488888888888889*yaux1*yaux15*yaux23*yaux24 - 0.488888888888889*yaux15*yaux24*yaux44*yaux48 + 8.88888888888889*yaux28*(yaux55 + yaux77) + 8.88888888888889*yaux41*(yaux72 + yaux80) + 8.88888888888889*yaux50*yaux86 + 8.88888888888889*yaux54*(-0.4 - 0.05500000000000001*yaux31*yaux33 + yaux36 + yaux37 + yaux89);
	out[3] = 8.88888888888889*yaux26*yaux28 + 8.88888888888889*yaux1*(yaux119 + yaux120 + yaux37 + yaux38) + 8.88888888888889*(yaux124 + yaux45)*yaux48 - 0.488888888888889*yaux24*yaux41*yaux44*yaux5 - 0.488888888888889*yaux24*yaux5*yaux52*yaux54 + 8.88888888888889*yaux50*(yaux132 + yaux56) + 8.88888888888889*yaux117*yaux96 - 0.488888888888889*yaux115*yaux123*yaux127*yaux99 - 0.488888888888889*yaux115*yaux129*yaux130*yaux99;
	out[4] = -0.488888888888889*yaux103*yaux115*yaux129*yaux130 + 8.88888888888889*yaux127*yaux143 - 0.488888888888889*yaux10*yaux23*yaux24*yaux28 - 0.488888888888889*yaux10*yaux24*yaux52*yaux54 + 8.88888888888889*yaux1*(yaux139 + yaux60) + 8.88888888888889*yaux41*yaux64 + 8.88888888888889*yaux48*(yaux119 + yaux145 + yaux37 + yaux67) + 8.88888888888889*yaux50*(yaux150 + yaux73) - 0.488888888888889*yaux103*yaux114*yaux115*yaux96;
	out[5] = -0.488888888888889*yaux107*yaux115*yaux123*yaux127 + 8.88888888888889*yaux129*yaux166 - 0.488888888888889*yaux15*yaux23*yaux24*yaux28 - 0.488888888888889*yaux15*yaux24*yaux41*yaux44 + 8.88888888888889*yaux1*(yaux157 + yaux77) + 8.88888888888889*yaux48*(yaux160 + yaux80) + 8.88888888888889*yaux54*yaux86 + 8.88888888888889*yaux50*(yaux119 + yaux168 + yaux37 + yaux89) - 0.488888888888889*yaux107*yaux114*yaux115*yaux96;
	out[6] = 8.88888888888889*yaux1*yaux117 + 8.88888888888889*yaux127*(yaux124 + 0.05500000000000001*yaux177*yaux193*yaux197) + 8.88888888888889*yaux129*(yaux132 + 0.05500000000000001*yaux177*yaux193*yaux201) + 8.88888888888889*(yaux119 + yaux120 + yaux191 + 0.05500000000000001*yaux177*yaux192*yaux193)*yaux96 - 0.488888888888889*yaux115*yaux123*yaux48*yaux99 - 0.488888888888889*yaux115*yaux130*yaux50*yaux99;
	out[7] = -0.488888888888889*yaux1*yaux103*yaux114*yaux115 + 8.88888888888889*yaux127*(yaux119 + yaux145 + yaux191 + 0.05500000000000001*yaux181*yaux193*yaux197) + 8.88888888888889*yaux129*(yaux150 + 0.05500000000000001*yaux181*yaux193*yaux201) + 8.88888888888889*yaux143*yaux48 - 0.488888888888889*yaux103*yaux115*yaux130*yaux50 + 8.88888888888889*(yaux139 + 0.05500000000000001*yaux181*yaux192*yaux193)*yaux96;
	out[8] = -0.488888888888889*yaux1*yaux107*yaux114*yaux115 + 8.88888888888889*yaux127*(yaux160 + 0.05500000000000001*yaux185*yaux193*yaux197) + 8.88888888888889*yaux129*(yaux119 + yaux168 + yaux191 + 0.05500000000000001*yaux185*yaux193*yaux201) - 0.488888888888889*yaux107*yaux115*yaux123*yaux48 + 8.88888888888889*yaux166*yaux50 + 8.88888888888889*(yaux157 + 0.05500000000000001*yaux185*yaux192*yaux193)*yaux96;
	out[9] = 0.488888888888889*yaux193*(1 * yaux127*yaux177*yaux181 + 1 * yaux129*yaux177*yaux185 + 1 * (1 * yaux178 + yaux237)*yaux96);
	out[10] = 0.488888888888889*yaux193*(1 * yaux129*yaux181*yaux185 + 1 * yaux127*(1 * yaux182 + yaux237) + 1 * yaux177*yaux181*yaux96);
	out[11] = -3.555555555555556*yaux193*(-0.1375*yaux127*yaux181*yaux185 + 1 * yaux129*(-0.1375*yaux186 + 1 * yaux187*yaux236) - 0.1375*yaux177*yaux185*yaux96);
	out[12] = vec[0];
	out[13] = vec[1];
	out[14] = vec[2];
	out[15] = vec[3];
	out[16] = vec[4];
	out[17] = vec[5];
	out[18] = vec[6];
	out[19] = vec[7];
	out[20] = vec[8];
}
/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	out[0] = vec[9];
	out[1] = vec[10];
	out[2] = vec[11];
}
/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
void dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}


/** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam) 
    -------------------------------------------------- **/
void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	ctypeRNum* param = (ctypeRNum*)userparam;
	typeRNum yaux11;

	yaux11 = x[9];

	out[0] = param[2] * (POW2(u[0]) + POW2(u[1]) + POW2(u[2]))
		+ param[1] * (56.25 - 15 * yaux11 + yaux11 * yaux11 + POW2(x[10]) + POW2(x[11]))
		+ param[0] * (POW2(x[12]) + POW2(x[13]) + POW2(x[14]) + POW2(x[15]) + POW2(x[16])
			+ POW2(x[17]) + POW2(x[18]) + POW2(x[19]) + POW2(x[20]));
}
/** Gradient dl/dx **/
void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	ctypeRNum* param = (ctypeRNum*)userparam;
	typeRNum yaux1, yaux9;

	yaux1 = param[1];
	yaux9 = param[0];

	out[0] = 0;
	out[1] = 0;
	out[2] = 0;
	out[3] = 0;
	out[4] = 0;
	out[5] = 0;
	out[6] = 0;
	out[7] = 0;
	out[8] = 0;
	out[9] = 2 * yaux1*(-7.5 + x[9]);
	out[10] = 2 * yaux1*x[10];
	out[11] = 2 * yaux1*x[11];
	out[12] = 2 * yaux9*x[12];
	out[13] = 2 * yaux9*x[13];
	out[14] = 2 * yaux9*x[14];
	out[15] = 2 * yaux9*x[15];
	out[16] = 2 * yaux9*x[16];
	out[17] = 2 * yaux9*x[17];
	out[18] = 2 * yaux9*x[18];
	out[19] = 2 * yaux9*x[19];
	out[20] = 2 * yaux9*x[20];
}
/** Gradient dl/du **/
void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	ctypeRNum* param = (ctypeRNum*)userparam;
	typeRNum yaux1;

	yaux1 = param[2];

	out[0] = 2 * yaux1*u[0];
	out[1] = 2 * yaux1*u[1];
	out[2] = 2 * yaux1*u[2];
}
/** Gradient dl/dp **/
void dldp(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
}


/** Terminal cost V(T,x(T),p,xdes,userparam) 
    ---------------------------------------- **/
void Vfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
	ctypeRNum* param = (ctypeRNum*)userparam;

	out[0] = param[3] * (POW2(-7.5 + x[9]) + POW2(x[10]) + POW2(x[11]));
}
/** Gradient dV/dx **/
void dVdx(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
	ctypeRNum* param = (ctypeRNum*)userparam;
	typeRNum yaux1;

	yaux1 = param[3];

	out[0] = 0;
	out[1] = 0;
	out[2] = 0;
	out[3] = 0;
	out[4] = 0;
	out[5] = 0;
	out[6] = 0;
	out[7] = 0;
	out[8] = 0;
	out[9] = 2 * yaux1*(-7.5 + x[9]);
	out[10] = 2 * yaux1*x[10];
	out[11] = 2 * yaux1*x[11];
	out[12] = 0;
	out[13] = 0;
	out[14] = 0;
	out[15] = 0;
	out[16] = 0;
	out[17] = 0;
	out[18] = 0;
	out[19] = 0;
	out[20] = 0;
}
/** Gradient dV/dp **/
void dVdp(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
}
/** Gradient dV/dT **/
void dVdT(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
	out[0] = 0;
}

/** Equality constraints g(t,x(t),u(t),p,uperparam) = 0 
    --------------------------------------------------- **/
void gfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian dg/dx multiplied by vector vec, i.e. (dg/dx)^T*vec or vec^T*(dg/dx) **/
void dgdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dg/du multiplied by vector vec, i.e. (dg/du)^T*vec or vec^T*(dg/du) **/
void dgdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dg/dp multiplied by vector vec, i.e. (dg/dp)^T*vec or vec^T*(dg/dp) **/
void dgdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}


/** Inequality constraints h(t,x(t),u(t),p,uperparam) <= 0 
    ------------------------------------------------------ **/
void hfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian dh/dx multiplied by vector vec, i.e. (dh/dx)^T*vec or vec^T*(dg/dx) **/
void dhdx_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dh/du multiplied by vector vec, i.e. (dh/du)^T*vec or vec^T*(dg/du) **/
void dhdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dh/dp multiplied by vector vec, i.e. (dh/dp)^T*vec or vec^T*(dg/dp) **/
void dhdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}


/** Terminal equality constraints gT(T,x(T),p,uperparam) = 0 
    -------------------------------------------------------- **/
void gTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian dgT/dx multiplied by vector vec, i.e. (dgT/dx)^T*vec or vec^T*(dgT/dx) **/
void dgTdx_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dgT/dp multiplied by vector vec, i.e. (dgT/dp)^T*vec or vec^T*(dgT/dp) **/
void dgTdp_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dgT/dT multiplied by vector vec, i.e. (dgT/dT)^T*vec or vec^T*(dgT/dT) **/
void dgTdT_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}


/** Terminal inequality constraints hT(T,x(T),p,uperparam) <= 0 
    ----------------------------------------------------------- **/
void hTfct(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian dhT/dx multiplied by vector vec, i.e. (dhT/dx)^T*vec or vec^T*(dhT/dx) **/
void dhTdx_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dhT/dp multiplied by vector vec, i.e. (dhT/dp)^T*vec or vec^T*(dhT/dp) **/
void dhTdp_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}
/** Jacobian dhT/dT multiplied by vector vec, i.e. (dhT/dT)^T*vec or vec^T*(dhT/dT) **/
void dhTdT_vec(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *vec, typeUSERPARAM *userparam)
{
}


/** Additional functions required for semi-implicit systems 
    M*dx/dt(t) = f(t0+t,x(t),u(t),p) using the solver RODAS 
    ------------------------------------------------------- **/
/** Jacobian df/dx in vector form (column-wise) **/
void dfdx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian df/dx in vector form (column-wise) **/
void dfdxtrans(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian df/dt **/
void dfdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Jacobian d(dH/dx)/dt  **/
void dHdxdt(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *vec, ctypeRNum *p, typeUSERPARAM *userparam)
{
}
/** Mass matrix in vector form (column-wise, either banded or full matrix) **/
void Mfct(typeRNum *out, typeUSERPARAM *userparam)
{
}
/** Transposed mass matrix in vector form (column-wise, either banded or full matrix) **/
void Mtrans(typeRNum *out, typeUSERPARAM *userparam)
{
}
