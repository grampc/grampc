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
	*Nx = 33;
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
	typeRNum yaux100, yaux101, yaux102, yaux103, yaux104, yaux105, yaux106, yaux107, yaux108, yaux109;
	typeRNum yaux110, yaux111, yaux112, yaux113, yaux118, yaux123, yaux129, yaux130, yaux131, yaux132;
	typeRNum yaux133, yaux134, yaux135, yaux136, yaux137, yaux138, yaux139, yaux140, yaux141, yaux142;
	typeRNum yaux143, yaux144, yaux145, yaux150, yaux155, yaux161, yaux163, yaux164, yaux165, yaux166;
	typeRNum yaux167, yaux168, yaux169, yaux170, yaux171, yaux172, yaux173, yaux174, yaux175, yaux176;
	typeRNum yaux177, yaux19, yaux20, yaux21, yaux22, yaux23, yaux24, yaux25, yaux26, yaux27;
	typeRNum yaux28, yaux30, yaux31, yaux32, yaux33, yaux34, yaux35, yaux36, yaux37, yaux38;
	typeRNum yaux39, yaux40, yaux41, yaux42, yaux43, yaux44, yaux45, yaux46, yaux47, yaux52;
	typeRNum yaux53, yaux58, yaux59, yaux65, yaux66, yaux67, yaux68, yaux69, yaux70, yaux71;
	typeRNum yaux72, yaux73, yaux74, yaux75, yaux76, yaux77, yaux78, yaux79, yaux80, yaux81;
	typeRNum yaux86, yaux91, yaux97, yaux98, yaux99;

	yaux19 = x[0];
	yaux31 = x[3];
	yaux21 = x[1];
	yaux23 = x[2];
	yaux20 = yaux19 * yaux19;
	yaux22 = yaux21 * yaux21;
	yaux24 = yaux23 * yaux23;
	yaux25 = yaux20 + yaux22 + yaux24;
	yaux26 = 1 / SQRT(yaux25);
	yaux27 = -0.55*yaux26;
	yaux28 = 6 + yaux27;
	yaux36 = x[4];
	yaux33 = -yaux31;
	yaux34 = yaux19 + yaux33;
	yaux35 = yaux34 * yaux34;
	yaux37 = -yaux36;
	yaux38 = yaux21 + yaux37;
	yaux39 = yaux38 * yaux38;
	yaux40 = x[5];
	yaux41 = -yaux40;
	yaux42 = yaux23 + yaux41;
	yaux43 = yaux42 * yaux42;
	yaux44 = yaux35 + yaux39 + yaux43;
	yaux45 = 1 / SQRT(yaux44);
	yaux46 = -0.55*yaux45;
	yaux47 = 6 + yaux46;
	yaux30 = -yaux19;
	yaux32 = yaux30 + yaux31;
	yaux65 = x[6];
	yaux52 = -yaux21;
	yaux53 = yaux36 + yaux52;
	yaux70 = x[7];
	yaux67 = -yaux65;
	yaux68 = yaux31 + yaux67;
	yaux69 = yaux68 * yaux68;
	yaux71 = -yaux70;
	yaux72 = yaux36 + yaux71;
	yaux73 = yaux72 * yaux72;
	yaux74 = x[8];
	yaux75 = -yaux74;
	yaux76 = yaux40 + yaux75;
	yaux77 = yaux76 * yaux76;
	yaux78 = yaux69 + yaux73 + yaux77;
	yaux79 = 1 / SQRT(yaux78);
	yaux80 = -0.55*yaux79;
	yaux81 = 6 + yaux80;
	yaux58 = -yaux23;
	yaux59 = yaux40 + yaux58;
	yaux66 = yaux33 + yaux65;
	yaux97 = x[9];
	yaux86 = yaux37 + yaux70;
	yaux102 = x[10];
	yaux99 = -yaux97;
	yaux100 = yaux65 + yaux99;
	yaux101 = yaux100 * yaux100;
	yaux103 = -yaux102;
	yaux104 = yaux103 + yaux70;
	yaux105 = yaux104 * yaux104;
	yaux106 = x[11];
	yaux107 = -yaux106;
	yaux108 = yaux107 + yaux74;
	yaux109 = yaux108 * yaux108;
	yaux110 = yaux101 + yaux105 + yaux109;
	yaux111 = 1 / SQRT(yaux110);
	yaux112 = -0.55*yaux111;
	yaux113 = 6 + yaux112;
	yaux91 = yaux41 + yaux74;
	yaux98 = yaux67 + yaux97;
	yaux129 = x[12];
	yaux118 = yaux102 + yaux71;
	yaux134 = x[13];
	yaux131 = -yaux129;
	yaux132 = yaux131 + yaux97;
	yaux133 = yaux132 * yaux132;
	yaux135 = -yaux134;
	yaux136 = yaux102 + yaux135;
	yaux137 = yaux136 * yaux136;
	yaux138 = x[14];
	yaux139 = -yaux138;
	yaux140 = yaux106 + yaux139;
	yaux141 = yaux140 * yaux140;
	yaux142 = yaux133 + yaux137 + yaux141;
	yaux143 = 1 / SQRT(yaux142);
	yaux144 = -0.55*yaux143;
	yaux145 = 6 + yaux144;
	yaux123 = yaux106 + yaux75;
	yaux130 = yaux129 + yaux99;
	yaux161 = x[15];
	yaux150 = yaux103 + yaux134;
	yaux166 = x[16];
	yaux163 = -yaux161;
	yaux164 = yaux129 + yaux163;
	yaux165 = yaux164 * yaux164;
	yaux167 = -yaux166;
	yaux168 = yaux134 + yaux167;
	yaux169 = yaux168 * yaux168;
	yaux170 = x[17];
	yaux171 = -yaux170;
	yaux172 = yaux138 + yaux171;
	yaux173 = yaux172 * yaux172;
	yaux174 = yaux165 + yaux169 + yaux173;
	yaux175 = 1 / SQRT(yaux174);
	yaux176 = -0.55*yaux175;
	yaux177 = 6 + yaux176;
	yaux155 = yaux107 + yaux138;

	out[0] = x[18];
	out[1] = x[19];
	out[2] = x[20];
	out[3] = x[21];
	out[4] = x[22];
	out[5] = x[23];
	out[6] = x[24];
	out[7] = x[25];
	out[8] = x[26];
	out[9] = x[27];
	out[10] = x[28];
	out[11] = x[29];
	out[12] = x[30];
	out[13] = x[31];
	out[14] = x[32];
	out[15] = u[0];
	out[16] = u[1];
	out[17] = u[2];
	out[18] = 13.333333333333334*(-0.1*yaux19*yaux28 + 0.1*yaux32*yaux47);
	out[19] = 13.333333333333334*(-0.1*yaux21*yaux28 + 0.1*yaux47*yaux53);
	out[20] = -9.81 + 13.333333333333334*(-0.1*yaux23*yaux28 + 0.1*yaux47*yaux59);
	out[21] = 13.333333333333334*(-0.1*yaux32*yaux47 + 0.1*yaux66*yaux81);
	out[22] = 13.333333333333334*(-0.1*yaux47*yaux53 + 0.1*yaux81*yaux86);
	out[23] = -9.81 + 13.333333333333334*(-0.1*yaux47*yaux59 + 0.1*yaux81*yaux91);
	out[24] = 13.333333333333334*(-0.1*yaux66*yaux81 + 0.1*yaux113*yaux98);
	out[25] = 13.333333333333334*(0.1*yaux113*yaux118 - 0.1*yaux81*yaux86);
	out[26] = -9.81 + 13.333333333333334*(0.1*yaux113*yaux123 - 0.1*yaux81*yaux91);
	out[27] = 13.333333333333334*(0.1*yaux130*yaux145 - 0.1*yaux113*yaux98);
	out[28] = 13.333333333333334*(-0.1*yaux113*yaux118 + 0.1*yaux145*yaux150);
	out[29] = -9.81 + 13.333333333333334*(-0.1*yaux113*yaux123 + 0.1*yaux145*yaux155);
	out[30] = 13.333333333333334*(-0.1*yaux130*yaux145 + 0.1*(yaux131 + yaux161)*yaux177);
	out[31] = 13.333333333333334*(-0.1*yaux145*yaux150 + 0.1*(yaux135 + yaux166)*yaux177);
	out[32] = -9.81 + 13.333333333333334*(-0.1*yaux145*yaux155 + 0.1*(yaux139 + yaux170)*yaux177);
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
	typeRNum yaux189, yaux19, yaux190, yaux191, yaux192, yaux193, yaux194, yaux195, yaux196, yaux198;
	typeRNum yaux199, yaux2, yaux20, yaux202, yaux203, yaux206, yaux208, yaux209, yaux21, yaux211;
	typeRNum yaux218, yaux22, yaux221, yaux222, yaux224, yaux229, yaux23, yaux236, yaux239, yaux24;
	typeRNum yaux244, yaux245, yaux247, yaux25, yaux254, yaux255, yaux256, yaux257, yaux258, yaux259;
	typeRNum yaux26, yaux260, yaux261, yaux262, yaux263, yaux264, yaux265, yaux266, yaux267, yaux268;
	typeRNum yaux269, yaux270, yaux271, yaux272, yaux273, yaux274, yaux275, yaux277, yaux278, yaux28;
	typeRNum yaux281, yaux282, yaux285, yaux287, yaux288, yaux29, yaux290, yaux297, yaux3, yaux30;
	typeRNum yaux300, yaux301, yaux303, yaux308, yaux31, yaux315, yaux318, yaux32, yaux323, yaux324;
	typeRNum yaux326, yaux33, yaux333, yaux334, yaux335, yaux336, yaux337, yaux338, yaux339, yaux340;
	typeRNum yaux341, yaux342, yaux343, yaux344, yaux345, yaux346, yaux347, yaux348, yaux349, yaux35;
	typeRNum yaux350, yaux351, yaux355, yaux359, yaux36, yaux37, yaux38, yaux392, yaux393, yaux394;
	typeRNum yaux395, yaux4, yaux41, yaux42, yaux43, yaux44, yaux45, yaux48, yaux5, yaux50;
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
	yaux20 = 6 + yaux19;
	yaux22 = -yaux2;
	yaux23 = yaux22 + yaux3;
	yaux24 = POW(yaux17, -1.5);
	yaux33 = POW(yaux32, -1.5);
	yaux43 = -yaux7;
	yaux44 = yaux43 + yaux8;
	yaux51 = -yaux12;
	yaux52 = yaux13 + yaux51;
	yaux28 = vec[18];
	yaux42 = -0.05500000000000001*yaux2*yaux33*yaux7;
	yaux48 = vec[22];
	yaux21 = 0.1*yaux20;
	yaux41 = vec[19];
	yaux35 = 1 / SQRT(yaux32);
	yaux36 = 0.05500000000000001*yaux35;
	yaux37 = -0.1*yaux20;
	yaux1 = vec[21];
	yaux50 = vec[23];
	yaux54 = vec[20];
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
	yaux112 = 6 + yaux111;
	yaux114 = yaux4 + yaux97;
	yaux115 = POW(yaux109, -1.5);
	yaux45 = 0.05500000000000001*yaux24*yaux44*yaux5;
	yaux123 = yaux101 + yaux9;
	yaux56 = 0.05500000000000001*yaux24*yaux5*yaux52;
	yaux130 = yaux105 + yaux14;
	yaux63 = -0.05500000000000001*yaux10*yaux24*yaux44;
	yaux64 = yaux21 + yaux63;
	yaux60 = 0.05500000000000001*yaux10*yaux23*yaux24;
	yaux127 = vec[25];
	yaux113 = 0.1*yaux112;
	yaux67 = 0.05500000000000001*yaux10*yaux24*yaux44;
	yaux119 = -0.1*yaux112;
	yaux96 = vec[24];
	yaux129 = vec[26];
	yaux73 = 0.05500000000000001*yaux10*yaux24*yaux52;
	yaux85 = -0.05500000000000001*yaux15*yaux24*yaux52;
	yaux86 = yaux21 + yaux85;
	yaux77 = 0.05500000000000001*yaux15*yaux23*yaux24;
	yaux80 = 0.05500000000000001*yaux15*yaux24*yaux44;
	yaux89 = 0.05500000000000001*yaux15*yaux24*yaux52;
	yaux116 = -0.05500000000000001*yaux114*yaux115*yaux99;
	yaux117 = yaux113 + yaux116;
	yaux176 = x[9];
	yaux177 = -yaux176;
	yaux178 = yaux177 + yaux97;
	yaux179 = yaux178 * yaux178;
	yaux180 = x[10];
	yaux181 = -yaux180;
	yaux182 = yaux101 + yaux181;
	yaux183 = yaux182 * yaux182;
	yaux184 = x[11];
	yaux185 = -yaux184;
	yaux186 = yaux105 + yaux185;
	yaux187 = yaux186 * yaux186;
	yaux188 = yaux179 + yaux183 + yaux187;
	yaux120 = 0.05500000000000001*yaux114*yaux115*yaux99;
	yaux189 = 1 / SQRT(yaux188);
	yaux190 = -0.55*yaux189;
	yaux191 = 6 + yaux190;
	yaux193 = yaux176 + yaux98;
	yaux194 = POW(yaux188, -1.5);
	yaux124 = 0.05500000000000001*yaux115*yaux123*yaux99;
	yaux202 = yaux102 + yaux180;
	yaux132 = 0.05500000000000001*yaux115*yaux130*yaux99;
	yaux209 = yaux106 + yaux184;
	yaux142 = -0.05500000000000001*yaux103*yaux115*yaux123;
	yaux143 = yaux113 + yaux142;
	yaux139 = 0.05500000000000001*yaux103*yaux114*yaux115;
	yaux206 = vec[28];
	yaux192 = 0.1*yaux191;
	yaux145 = 0.05500000000000001*yaux103*yaux115*yaux123;
	yaux198 = -0.1*yaux191;
	yaux175 = vec[27];
	yaux208 = vec[29];
	yaux150 = 0.05500000000000001*yaux103*yaux115*yaux130;
	yaux165 = -0.05500000000000001*yaux107*yaux115*yaux130;
	yaux166 = yaux113 + yaux165;
	yaux157 = 0.05500000000000001*yaux107*yaux114*yaux115;
	yaux160 = 0.05500000000000001*yaux107*yaux115*yaux123;
	yaux168 = 0.05500000000000001*yaux107*yaux115*yaux130;
	yaux195 = -0.05500000000000001*yaux178*yaux193*yaux194;
	yaux196 = yaux192 + yaux195;
	yaux255 = x[12];
	yaux256 = -yaux255;
	yaux257 = yaux176 + yaux256;
	yaux258 = yaux257 * yaux257;
	yaux259 = x[13];
	yaux260 = -yaux259;
	yaux261 = yaux180 + yaux260;
	yaux262 = yaux261 * yaux261;
	yaux263 = x[14];
	yaux264 = -yaux263;
	yaux265 = yaux184 + yaux264;
	yaux266 = yaux265 * yaux265;
	yaux267 = yaux258 + yaux262 + yaux266;
	yaux199 = 0.05500000000000001*yaux178*yaux193*yaux194;
	yaux268 = 1 / SQRT(yaux267);
	yaux269 = -0.55*yaux268;
	yaux270 = 6 + yaux269;
	yaux272 = yaux177 + yaux255;
	yaux273 = POW(yaux267, -1.5);
	yaux203 = 0.05500000000000001*yaux178*yaux194*yaux202;
	yaux281 = yaux181 + yaux259;
	yaux211 = 0.05500000000000001*yaux178*yaux194*yaux209;
	yaux288 = yaux185 + yaux263;
	yaux221 = -0.05500000000000001*yaux182*yaux194*yaux202;
	yaux222 = yaux192 + yaux221;
	yaux218 = 0.05500000000000001*yaux182*yaux193*yaux194;
	yaux285 = vec[31];
	yaux271 = 0.1*yaux270;
	yaux224 = 0.05500000000000001*yaux182*yaux194*yaux202;
	yaux277 = -0.1*yaux270;
	yaux254 = vec[30];
	yaux287 = vec[32];
	yaux229 = 0.05500000000000001*yaux182*yaux194*yaux209;
	yaux244 = -0.05500000000000001*yaux186*yaux194*yaux209;
	yaux245 = yaux192 + yaux244;
	yaux236 = 0.05500000000000001*yaux186*yaux193*yaux194;
	yaux239 = 0.05500000000000001*yaux186*yaux194*yaux202;
	yaux247 = 0.05500000000000001*yaux186*yaux194*yaux209;
	yaux274 = -0.05500000000000001*yaux257*yaux272*yaux273;
	yaux275 = yaux271 + yaux274;
	yaux278 = 0.05500000000000001*yaux257*yaux272*yaux273;
	yaux333 = x[15];
	yaux334 = -yaux333;
	yaux335 = yaux255 + yaux334;
	yaux336 = yaux335 * yaux335;
	yaux337 = x[16];
	yaux338 = -yaux337;
	yaux339 = yaux259 + yaux338;
	yaux340 = yaux339 * yaux339;
	yaux341 = x[17];
	yaux342 = -yaux341;
	yaux343 = yaux263 + yaux342;
	yaux344 = yaux343 * yaux343;
	yaux345 = yaux336 + yaux340 + yaux344;
	yaux282 = 0.05500000000000001*yaux257*yaux273*yaux281;
	yaux351 = POW(yaux345, -1.5);
	yaux290 = 0.05500000000000001*yaux257*yaux273*yaux288;
	yaux300 = -0.05500000000000001*yaux261*yaux273*yaux281;
	yaux301 = yaux271 + yaux300;
	yaux297 = 0.05500000000000001*yaux261*yaux272*yaux273;
	yaux350 = yaux256 + yaux333;
	yaux303 = 0.05500000000000001*yaux261*yaux273*yaux281;
	yaux346 = 1 / SQRT(yaux345);
	yaux347 = -0.55*yaux346;
	yaux348 = 6 + yaux347;
	yaux349 = -0.1*yaux348;
	yaux355 = yaux260 + yaux337;
	yaux308 = 0.05500000000000001*yaux261*yaux273*yaux288;
	yaux359 = yaux264 + yaux341;
	yaux323 = -0.05500000000000001*yaux265*yaux273*yaux288;
	yaux324 = yaux271 + yaux323;
	yaux315 = 0.05500000000000001*yaux265*yaux272*yaux273;
	yaux318 = 0.05500000000000001*yaux265*yaux273*yaux281;
	yaux326 = 0.05500000000000001*yaux265*yaux273*yaux288;
	yaux392 = SQRT(yaux345);
	yaux393 = -yaux392;
	yaux394 = 0.09166666666666667 + yaux393;
	yaux395 = -10.90909090909091*yaux345*yaux394;

	out[0] = 13.333333333333334*yaux1*yaux26 + 13.333333333333334*yaux28*(-0.6000000000000001 - 0.05500000000000001*yaux29*yaux33 + yaux36 + yaux37 + yaux38) + 13.333333333333334*yaux41*(yaux42 + yaux45) - 0.7333333333333335*yaux24*yaux44*yaux48*yaux5 - 0.7333333333333335*yaux24*yaux5*yaux50*yaux52 + 13.333333333333334*yaux54*(yaux55 + yaux56);
	out[1] = -0.7333333333333335*yaux1*yaux10*yaux23*yaux24 - 0.7333333333333335*yaux10*yaux24*yaux50*yaux52 + 13.333333333333334*yaux28*(yaux42 + yaux60) + 13.333333333333334*yaux48*yaux64 + 13.333333333333334*yaux41*(-0.6000000000000001 - 0.05500000000000001*yaux30*yaux33 + yaux36 + yaux37 + yaux67) + 13.333333333333334*yaux54*(yaux72 + yaux73);
	out[2] = -0.7333333333333335*yaux1*yaux15*yaux23*yaux24 - 0.7333333333333335*yaux15*yaux24*yaux44*yaux48 + 13.333333333333334*yaux28*(yaux55 + yaux77) + 13.333333333333334*yaux41*(yaux72 + yaux80) + 13.333333333333334*yaux50*yaux86 + 13.333333333333334*yaux54*(-0.6000000000000001 - 0.05500000000000001*yaux31*yaux33 + yaux36 + yaux37 + yaux89);
	out[3] = 13.333333333333334*yaux26*yaux28 + 13.333333333333334*yaux1*(yaux119 + yaux120 + yaux37 + yaux38) + 13.333333333333334*(yaux124 + yaux45)*yaux48 - 0.7333333333333335*yaux24*yaux41*yaux44*yaux5 - 0.7333333333333335*yaux24*yaux5*yaux52*yaux54 + 13.333333333333334*yaux50*(yaux132 + yaux56) + 13.333333333333334*yaux117*yaux96 - 0.7333333333333335*yaux115*yaux123*yaux127*yaux99 - 0.7333333333333335*yaux115*yaux129*yaux130*yaux99;
	out[4] = -0.7333333333333335*yaux103*yaux115*yaux129*yaux130 + 13.333333333333334*yaux127*yaux143 - 0.7333333333333335*yaux10*yaux23*yaux24*yaux28 - 0.7333333333333335*yaux10*yaux24*yaux52*yaux54 + 13.333333333333334*yaux1*(yaux139 + yaux60) + 13.333333333333334*yaux41*yaux64 + 13.333333333333334*yaux48*(yaux119 + yaux145 + yaux37 + yaux67) + 13.333333333333334*yaux50*(yaux150 + yaux73) - 0.7333333333333335*yaux103*yaux114*yaux115*yaux96;
	out[5] = -0.7333333333333335*yaux107*yaux115*yaux123*yaux127 + 13.333333333333334*yaux129*yaux166 - 0.7333333333333335*yaux15*yaux23*yaux24*yaux28 - 0.7333333333333335*yaux15*yaux24*yaux41*yaux44 + 13.333333333333334*yaux1*(yaux157 + yaux77) + 13.333333333333334*yaux48*(yaux160 + yaux80) + 13.333333333333334*yaux54*yaux86 + 13.333333333333334*yaux50*(yaux119 + yaux168 + yaux37 + yaux89) - 0.7333333333333335*yaux107*yaux114*yaux115*yaux96;
	out[6] = 13.333333333333334*yaux1*yaux117 + 13.333333333333334*yaux175*yaux196 + 13.333333333333334*yaux127*(yaux124 + yaux203) - 0.7333333333333335*yaux178*yaux194*yaux202*yaux206 - 0.7333333333333335*yaux178*yaux194*yaux208*yaux209 + 13.333333333333334*yaux129*(yaux132 + yaux211) + 13.333333333333334*(yaux119 + yaux120 + yaux198 + yaux199)*yaux96 - 0.7333333333333335*yaux115*yaux123*yaux48*yaux99 - 0.7333333333333335*yaux115*yaux130*yaux50*yaux99;
	out[7] = -0.7333333333333335*yaux1*yaux103*yaux114*yaux115 - 0.7333333333333335*yaux175*yaux182*yaux193*yaux194 - 0.7333333333333335*yaux182*yaux194*yaux208*yaux209 + 13.333333333333334*yaux206*yaux222 + 13.333333333333334*yaux127*(yaux119 + yaux145 + yaux198 + yaux224) + 13.333333333333334*yaux129*(yaux150 + yaux229) + 13.333333333333334*yaux143*yaux48 - 0.7333333333333335*yaux103*yaux115*yaux130*yaux50 + 13.333333333333334*(yaux139 + yaux218)*yaux96;
	out[8] = -0.7333333333333335*yaux1*yaux107*yaux114*yaux115 - 0.7333333333333335*yaux175*yaux186*yaux193*yaux194 - 0.7333333333333335*yaux186*yaux194*yaux202*yaux206 + 13.333333333333334*yaux127*(yaux160 + yaux239) + 13.333333333333334*yaux208*yaux245 + 13.333333333333334*yaux129*(yaux119 + yaux168 + yaux198 + yaux247) - 0.7333333333333335*yaux107*yaux115*yaux123*yaux48 + 13.333333333333334*yaux166*yaux50 + 13.333333333333334*(yaux157 + yaux236)*yaux96;
	out[9] = -0.7333333333333335*yaux127*yaux178*yaux194*yaux202 - 0.7333333333333335*yaux129*yaux178*yaux194*yaux209 + 13.333333333333334*yaux254*yaux275 + 13.333333333333334*yaux175*(yaux198 + yaux199 + yaux277 + yaux278) + 13.333333333333334*yaux206*(yaux203 + yaux282) - 0.7333333333333335*yaux257*yaux273*yaux281*yaux285 - 0.7333333333333335*yaux257*yaux273*yaux287*yaux288 + 13.333333333333334*yaux208*(yaux211 + yaux290) + 13.333333333333334*yaux196*yaux96;
	out[10] = -0.7333333333333335*yaux129*yaux182*yaux194*yaux209 + 13.333333333333334*yaux127*yaux222 - 0.7333333333333335*yaux254*yaux261*yaux272*yaux273 - 0.7333333333333335*yaux261*yaux273*yaux287*yaux288 + 13.333333333333334*yaux175*(yaux218 + yaux297) + 13.333333333333334*yaux285*yaux301 + 13.333333333333334*yaux206*(yaux198 + yaux224 + yaux277 + yaux303) + 13.333333333333334*yaux208*(yaux229 + yaux308) - 0.7333333333333335*yaux182*yaux193*yaux194*yaux96;
	out[11] = -0.7333333333333335*yaux127*yaux186*yaux194*yaux202 + 13.333333333333334*yaux129*yaux245 - 0.7333333333333335*yaux254*yaux265*yaux272*yaux273 - 0.7333333333333335*yaux265*yaux273*yaux281*yaux285 + 13.333333333333334*yaux175*(yaux236 + yaux315) + 13.333333333333334*yaux206*(yaux239 + yaux318) + 13.333333333333334*yaux287*yaux324 + 13.333333333333334*yaux208*(yaux198 + yaux247 + yaux277 + yaux326) - 0.7333333333333335*yaux186*yaux193*yaux194*yaux96;
	out[12] = 13.333333333333334*yaux175*yaux275 - 0.7333333333333335*yaux206*yaux257*yaux273*yaux281 - 0.7333333333333335*yaux208*yaux257*yaux273*yaux288 + 13.333333333333334*yaux254*(yaux277 + yaux278 + yaux349 + 0.05500000000000001*yaux335*yaux350*yaux351) + 13.333333333333334*yaux285*(yaux282 + 0.05500000000000001*yaux335*yaux351*yaux355) + 13.333333333333334*yaux287*(yaux290 + 0.05500000000000001*yaux335*yaux351*yaux359);
	out[13] = -0.7333333333333335*yaux175*yaux261*yaux272*yaux273 - 0.7333333333333335*yaux208*yaux261*yaux273*yaux288 + 13.333333333333334*yaux206*yaux301 + 13.333333333333334*yaux254*(yaux297 + 0.05500000000000001*yaux339*yaux350*yaux351) + 13.333333333333334*yaux285*(yaux277 + yaux303 + yaux349 + 0.05500000000000001*yaux339*yaux351*yaux355) + 13.333333333333334*yaux287*(yaux308 + 0.05500000000000001*yaux339*yaux351*yaux359);
	out[14] = -0.7333333333333335*yaux175*yaux265*yaux272*yaux273 - 0.7333333333333335*yaux206*yaux265*yaux273*yaux281 + 13.333333333333334*yaux208*yaux324 + 13.333333333333334*yaux254*(yaux315 + 0.05500000000000001*yaux343*yaux350*yaux351) + 13.333333333333334*yaux285*(yaux318 + 0.05500000000000001*yaux343*yaux351*yaux355) + 13.333333333333334*yaux287*(yaux277 + yaux326 + yaux349 + 0.05500000000000001*yaux343*yaux351*yaux359);
	out[15] = 0.7333333333333335*yaux351*(1 * yaux285*yaux335*yaux339 + 1 * yaux287*yaux335*yaux343 + 1 * yaux254*(1 * yaux336 + yaux395));
	out[16] = 0.7333333333333335*yaux351*(1 * yaux254*yaux335*yaux339 + 1 * yaux287*yaux339*yaux343 + 1 * yaux285*(1 * yaux340 + yaux395));
	out[17] = -8.000000000000002*yaux351*(-0.09166666666666667*yaux254*yaux335*yaux343 - 0.09166666666666667*yaux285*yaux339*yaux343 + 1 * yaux287*(-0.09166666666666666*yaux344 + 1 * yaux345*yaux394));
	out[18] = vec[0];
	out[19] = vec[1];
	out[20] = vec[2];
	out[21] = vec[3];
	out[22] = vec[4];
	out[23] = vec[5];
	out[24] = vec[6];
	out[25] = vec[7];
	out[26] = vec[8];
	out[27] = vec[9];
	out[28] = vec[10];
	out[29] = vec[11];
	out[30] = vec[12];
	out[31] = vec[13];
	out[32] = vec[14];
}
/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	out[0] = vec[15];
	out[1] = vec[16];
	out[2] = vec[17];
}
/** Jacobian df/dp multiplied by vector vec, i.e. (df/dp)^T*vec or vec^T*(df/dp) **/
void dfdp_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
}


/** Integral cost l(t,x(t),u(t),p,xdes,udes,userparam) 
    -------------------------------------------------- **/
void lfct(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	typeRNum* param = (typeRNum*)userparam;
	typeRNum yaux11;

	yaux11 = x[15];

	/* Number of output variables: 1 */
	out[0] = param[2] * (POW2(u[0]) + POW2(u[1]) + POW2(u[2]))
		+ param[1] * (56.25 - 15 * yaux11 + yaux11 * yaux11 + POW2(x[16]) + POW2(x[17]))
		+ param[0] * (POW2(x[18]) + POW2(x[19]) + POW2(x[20]) + POW2(x[21]) + POW2(x[22]) + POW2(x[23]) + POW2(x[24]) + POW2(x[25]) + POW2(x[26]) + POW2(x[27]) + POW2(x[28]) + POW2(x[29]) + POW2(x[30]) + POW2(x[31]) + POW2(x[32]));
}
/** Gradient dl/dx **/
void dldx(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	typeRNum* param = (typeRNum*)userparam;
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
	out[9] = 0;
	out[10] = 0;
	out[11] = 0;
	out[12] = 0;
	out[13] = 0;
	out[14] = 0;
	out[15] = 2 * yaux1*(-7.5 + x[15]);
	out[16] = 2 * yaux1*x[16];
	out[17] = 2 * yaux1*x[17];
	out[18] = 2 * yaux9*x[18];
	out[19] = 2 * yaux9*x[19];
	out[20] = 2 * yaux9*x[20];
	out[21] = 2 * yaux9*x[21];
	out[22] = 2 * yaux9*x[22];
	out[23] = 2 * yaux9*x[23];
	out[24] = 2 * yaux9*x[24];
	out[25] = 2 * yaux9*x[25];
	out[26] = 2 * yaux9*x[26];
	out[27] = 2 * yaux9*x[27];
	out[28] = 2 * yaux9*x[28];
	out[29] = 2 * yaux9*x[29];
	out[30] = 2 * yaux9*x[30];
	out[31] = 2 * yaux9*x[31];
	out[32] = 2 * yaux9*x[32];
}
/** Gradient dl/du **/
void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	typeRNum* param = (typeRNum*)userparam;
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
	typeRNum* param = (typeRNum*)userparam;

	out[0] = param[3] * (POW2(-7.5 + x[15]) + POW2(x[16]) + POW2(x[17]));
}
/** Gradient dV/dx **/
void dVdx(typeRNum *out, ctypeRNum T, ctypeRNum *x, ctypeRNum *p, ctypeRNum *xdes, typeUSERPARAM *userparam)
{
	typeRNum* param = (typeRNum*)userparam;
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
	out[9] = 0;
	out[10] = 0;
	out[11] = 0;
	out[12] = 0;
	out[13] = 0;
	out[14] = 0;
	out[15] = 2 * yaux1*(-7.5 + x[15]);
	out[16] = 2 * yaux1*x[16];
	out[17] = 2 * yaux1*x[17];
	out[18] = 0;
	out[19] = 0;
	out[20] = 0;
	out[21] = 0;
	out[22] = 0;
	out[23] = 0;
	out[24] = 0;
	out[25] = 0;
	out[26] = 0;
	out[27] = 0;
	out[28] = 0;
	out[29] = 0;
	out[30] = 0;
	out[31] = 0;
	out[32] = 0;
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
