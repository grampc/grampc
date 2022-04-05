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
	*Nx = 57;
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
	typeRNum yaux103, yaux109, yaux110, yaux111, yaux112, yaux113, yaux114, yaux115, yaux116, yaux117;
	typeRNum yaux118, yaux119, yaux120, yaux121, yaux122, yaux123, yaux124, yaux125, yaux130, yaux135;
	typeRNum yaux141, yaux142, yaux143, yaux144, yaux145, yaux146, yaux147, yaux148, yaux149, yaux150;
	typeRNum yaux151, yaux152, yaux153, yaux154, yaux155, yaux156, yaux157, yaux162, yaux167, yaux173;
	typeRNum yaux174, yaux175, yaux176, yaux177, yaux178, yaux179, yaux180, yaux181, yaux182, yaux183;
	typeRNum yaux184, yaux185, yaux186, yaux187, yaux188, yaux189, yaux194, yaux199, yaux205, yaux206;
	typeRNum yaux207, yaux208, yaux209, yaux210, yaux211, yaux212, yaux213, yaux214, yaux215, yaux216;
	typeRNum yaux217, yaux218, yaux219, yaux220, yaux221, yaux226, yaux231, yaux237, yaux238, yaux239;
	typeRNum yaux240, yaux241, yaux242, yaux243, yaux244, yaux245, yaux246, yaux247, yaux248, yaux249;
	typeRNum yaux250, yaux251, yaux252, yaux253, yaux258, yaux263, yaux269, yaux270, yaux271, yaux272;
	typeRNum yaux273, yaux274, yaux275, yaux276, yaux277, yaux278, yaux279, yaux280, yaux281, yaux282;
	typeRNum yaux283, yaux284, yaux285, yaux290, yaux295, yaux301, yaux303, yaux304, yaux305, yaux306;
	typeRNum yaux307, yaux308, yaux309, yaux31, yaux310, yaux311, yaux312, yaux313, yaux314, yaux315;
	typeRNum yaux316, yaux317, yaux32, yaux33, yaux34, yaux35, yaux36, yaux37, yaux38, yaux39;
	typeRNum yaux40, yaux42, yaux43, yaux44, yaux45, yaux46, yaux47, yaux48, yaux49, yaux50;
	typeRNum yaux51, yaux52, yaux53, yaux54, yaux55, yaux56, yaux57, yaux58, yaux59, yaux64;
	typeRNum yaux65, yaux70, yaux71, yaux77, yaux78, yaux79, yaux80, yaux81, yaux82, yaux83;
	typeRNum yaux84, yaux85, yaux86, yaux87, yaux88, yaux89, yaux90, yaux91, yaux92, yaux93;
	typeRNum yaux98;

	yaux31 = x[0];
	yaux43 = x[3];
	yaux33 = x[1];
	yaux35 = x[2];
	yaux32 = yaux31 * yaux31;
	yaux34 = yaux33 * yaux33;
	yaux36 = yaux35 * yaux35;
	yaux37 = yaux32 + yaux34 + yaux36;
	yaux38 = 1 / SQRT(yaux37);
	yaux39 = -0.55*yaux38;
	yaux40 = 10 + yaux39;
	yaux48 = x[4];
	yaux45 = -yaux43;
	yaux46 = yaux31 + yaux45;
	yaux47 = yaux46 * yaux46;
	yaux49 = -yaux48;
	yaux50 = yaux33 + yaux49;
	yaux51 = yaux50 * yaux50;
	yaux52 = x[5];
	yaux53 = -yaux52;
	yaux54 = yaux35 + yaux53;
	yaux55 = yaux54 * yaux54;
	yaux56 = yaux47 + yaux51 + yaux55;
	yaux57 = 1 / SQRT(yaux56);
	yaux58 = -0.55*yaux57;
	yaux59 = 10 + yaux58;
	yaux42 = -yaux31;
	yaux44 = yaux42 + yaux43;
	yaux77 = x[6];
	yaux64 = -yaux33;
	yaux65 = yaux48 + yaux64;
	yaux82 = x[7];
	yaux79 = -yaux77;
	yaux80 = yaux43 + yaux79;
	yaux81 = yaux80 * yaux80;
	yaux83 = -yaux82;
	yaux84 = yaux48 + yaux83;
	yaux85 = yaux84 * yaux84;
	yaux86 = x[8];
	yaux87 = -yaux86;
	yaux88 = yaux52 + yaux87;
	yaux89 = yaux88 * yaux88;
	yaux90 = yaux81 + yaux85 + yaux89;
	yaux91 = 1 / SQRT(yaux90);
	yaux92 = -0.55*yaux91;
	yaux93 = 10 + yaux92;
	yaux70 = -yaux35;
	yaux71 = yaux52 + yaux70;
	yaux78 = yaux45 + yaux77;
	yaux109 = x[9];
	yaux98 = yaux49 + yaux82;
	yaux114 = x[10];
	yaux111 = -yaux109;
	yaux112 = yaux111 + yaux77;
	yaux113 = yaux112 * yaux112;
	yaux115 = -yaux114;
	yaux116 = yaux115 + yaux82;
	yaux117 = yaux116 * yaux116;
	yaux118 = x[11];
	yaux119 = -yaux118;
	yaux120 = yaux119 + yaux86;
	yaux121 = yaux120 * yaux120;
	yaux122 = yaux113 + yaux117 + yaux121;
	yaux123 = 1 / SQRT(yaux122);
	yaux124 = -0.55*yaux123;
	yaux125 = 10 + yaux124;
	yaux103 = yaux53 + yaux86;
	yaux110 = yaux109 + yaux79;
	yaux141 = x[12];
	yaux130 = yaux114 + yaux83;
	yaux146 = x[13];
	yaux143 = -yaux141;
	yaux144 = yaux109 + yaux143;
	yaux145 = yaux144 * yaux144;
	yaux147 = -yaux146;
	yaux148 = yaux114 + yaux147;
	yaux149 = yaux148 * yaux148;
	yaux150 = x[14];
	yaux151 = -yaux150;
	yaux152 = yaux118 + yaux151;
	yaux153 = yaux152 * yaux152;
	yaux154 = yaux145 + yaux149 + yaux153;
	yaux155 = 1 / SQRT(yaux154);
	yaux156 = -0.55*yaux155;
	yaux157 = 10 + yaux156;
	yaux135 = yaux118 + yaux87;
	yaux142 = yaux111 + yaux141;
	yaux173 = x[15];
	yaux162 = yaux115 + yaux146;
	yaux178 = x[16];
	yaux175 = -yaux173;
	yaux176 = yaux141 + yaux175;
	yaux177 = yaux176 * yaux176;
	yaux179 = -yaux178;
	yaux180 = yaux146 + yaux179;
	yaux181 = yaux180 * yaux180;
	yaux182 = x[17];
	yaux183 = -yaux182;
	yaux184 = yaux150 + yaux183;
	yaux185 = yaux184 * yaux184;
	yaux186 = yaux177 + yaux181 + yaux185;
	yaux187 = 1 / SQRT(yaux186);
	yaux188 = -0.55*yaux187;
	yaux189 = 10 + yaux188;
	yaux167 = yaux119 + yaux150;
	yaux174 = yaux143 + yaux173;
	yaux205 = x[18];
	yaux194 = yaux147 + yaux178;
	yaux210 = x[19];
	yaux207 = -yaux205;
	yaux208 = yaux173 + yaux207;
	yaux209 = yaux208 * yaux208;
	yaux211 = -yaux210;
	yaux212 = yaux178 + yaux211;
	yaux213 = yaux212 * yaux212;
	yaux214 = x[20];
	yaux215 = -yaux214;
	yaux216 = yaux182 + yaux215;
	yaux217 = yaux216 * yaux216;
	yaux218 = yaux209 + yaux213 + yaux217;
	yaux219 = 1 / SQRT(yaux218);
	yaux220 = -0.55*yaux219;
	yaux221 = 10 + yaux220;
	yaux199 = yaux151 + yaux182;
	yaux206 = yaux175 + yaux205;
	yaux237 = x[21];
	yaux226 = yaux179 + yaux210;
	yaux242 = x[22];
	yaux239 = -yaux237;
	yaux240 = yaux205 + yaux239;
	yaux241 = yaux240 * yaux240;
	yaux243 = -yaux242;
	yaux244 = yaux210 + yaux243;
	yaux245 = yaux244 * yaux244;
	yaux246 = x[23];
	yaux247 = -yaux246;
	yaux248 = yaux214 + yaux247;
	yaux249 = yaux248 * yaux248;
	yaux250 = yaux241 + yaux245 + yaux249;
	yaux251 = 1 / SQRT(yaux250);
	yaux252 = -0.55*yaux251;
	yaux253 = 10 + yaux252;
	yaux231 = yaux183 + yaux214;
	yaux238 = yaux207 + yaux237;
	yaux269 = x[24];
	yaux258 = yaux211 + yaux242;
	yaux274 = x[25];
	yaux271 = -yaux269;
	yaux272 = yaux237 + yaux271;
	yaux273 = yaux272 * yaux272;
	yaux275 = -yaux274;
	yaux276 = yaux242 + yaux275;
	yaux277 = yaux276 * yaux276;
	yaux278 = x[26];
	yaux279 = -yaux278;
	yaux280 = yaux246 + yaux279;
	yaux281 = yaux280 * yaux280;
	yaux282 = yaux273 + yaux277 + yaux281;
	yaux283 = 1 / SQRT(yaux282);
	yaux284 = -0.55*yaux283;
	yaux285 = 10 + yaux284;
	yaux263 = yaux215 + yaux246;
	yaux270 = yaux239 + yaux269;
	yaux301 = x[27];
	yaux290 = yaux243 + yaux274;
	yaux306 = x[28];
	yaux303 = -yaux301;
	yaux304 = yaux269 + yaux303;
	yaux305 = yaux304 * yaux304;
	yaux307 = -yaux306;
	yaux308 = yaux274 + yaux307;
	yaux309 = yaux308 * yaux308;
	yaux310 = x[29];
	yaux311 = -yaux310;
	yaux312 = yaux278 + yaux311;
	yaux313 = yaux312 * yaux312;
	yaux314 = yaux305 + yaux309 + yaux313;
	yaux315 = 1 / SQRT(yaux314);
	yaux316 = -0.55*yaux315;
	yaux317 = 10 + yaux316;
	yaux295 = yaux247 + yaux278;

	out[0] = x[30];
	out[1] = x[31];
	out[2] = x[32];
	out[3] = x[33];
	out[4] = x[34];
	out[5] = x[35];
	out[6] = x[36];
	out[7] = x[37];
	out[8] = x[38];
	out[9] = x[39];
	out[10] = x[40];
	out[11] = x[41];
	out[12] = x[42];
	out[13] = x[43];
	out[14] = x[44];
	out[15] = x[45];
	out[16] = x[46];
	out[17] = x[47];
	out[18] = x[48];
	out[19] = x[49];
	out[20] = x[50];
	out[21] = x[51];
	out[22] = x[52];
	out[23] = x[53];
	out[24] = x[54];
	out[25] = x[55];
	out[26] = x[56];
	out[27] = u[0];
	out[28] = u[1];
	out[29] = u[2];
	out[30] = 22.22222222222222*(-0.1*yaux31*yaux40 + 0.1*yaux44*yaux59);
	out[31] = 22.22222222222222*(-0.1*yaux33*yaux40 + 0.1*yaux59*yaux65);
	out[32] = -9.81 + 22.22222222222222*(-0.1*yaux35*yaux40 + 0.1*yaux59*yaux71);
	out[33] = 22.22222222222222*(-0.1*yaux44*yaux59 + 0.1*yaux78*yaux93);
	out[34] = 22.22222222222222*(-0.1*yaux59*yaux65 + 0.1*yaux93*yaux98);
	out[35] = -9.81 + 22.22222222222222*(-0.1*yaux59*yaux71 + 0.1*yaux103*yaux93);
	out[36] = 22.22222222222222*(0.1*yaux110*yaux125 - 0.1*yaux78*yaux93);
	out[37] = 22.22222222222222*(0.1*yaux125*yaux130 - 0.1*yaux93*yaux98);
	out[38] = -9.81 + 22.22222222222222*(0.1*yaux125*yaux135 - 0.1*yaux103*yaux93);
	out[39] = 22.22222222222222*(-0.1*yaux110*yaux125 + 0.1*yaux142*yaux157);
	out[40] = 22.22222222222222*(-0.1*yaux125*yaux130 + 0.1*yaux157*yaux162);
	out[41] = -9.81 + 22.22222222222222*(-0.1*yaux125*yaux135 + 0.1*yaux157*yaux167);
	out[42] = 22.22222222222222*(-0.1*yaux142*yaux157 + 0.1*yaux174*yaux189);
	out[43] = 22.22222222222222*(-0.1*yaux157*yaux162 + 0.1*yaux189*yaux194);
	out[44] = -9.81 + 22.22222222222222*(-0.1*yaux157*yaux167 + 0.1*yaux189*yaux199);
	out[45] = 22.22222222222222*(-0.1*yaux174*yaux189 + 0.1*yaux206*yaux221);
	out[46] = 22.22222222222222*(-0.1*yaux189*yaux194 + 0.1*yaux221*yaux226);
	out[47] = -9.81 + 22.22222222222222*(-0.1*yaux189*yaux199 + 0.1*yaux221*yaux231);
	out[48] = 22.22222222222222*(-0.1*yaux206*yaux221 + 0.1*yaux238*yaux253);
	out[49] = 22.22222222222222*(-0.1*yaux221*yaux226 + 0.1*yaux253*yaux258);
	out[50] = -9.81 + 22.22222222222222*(-0.1*yaux221*yaux231 + 0.1*yaux253*yaux263);
	out[51] = 22.22222222222222*(-0.1*yaux238*yaux253 + 0.1*yaux270*yaux285);
	out[52] = 22.22222222222222*(-0.1*yaux253*yaux258 + 0.1*yaux285*yaux290);
	out[53] = -9.81 + 22.22222222222222*(-0.1*yaux253*yaux263 + 0.1*yaux285*yaux295);
	out[54] = 22.22222222222222*(-0.1*yaux270*yaux285 + 0.1*(yaux271 + yaux301)*yaux317);
	out[55] = 22.22222222222222*(-0.1*yaux285*yaux290 + 0.1*(yaux275 + yaux306)*yaux317);
	out[56] = -9.81 + 22.22222222222222*(-0.1*yaux285*yaux295 + 0.1*(yaux279 + yaux310)*yaux317);
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
	typeRNum yaux350, yaux351, yaux352, yaux353, yaux354, yaux356, yaux357, yaux36, yaux360, yaux361;
	typeRNum yaux364, yaux366, yaux367, yaux369, yaux37, yaux376, yaux379, yaux38, yaux380, yaux382;
	typeRNum yaux387, yaux394, yaux397, yaux4, yaux402, yaux403, yaux405, yaux41, yaux412, yaux413;
	typeRNum yaux414, yaux415, yaux416, yaux417, yaux418, yaux419, yaux42, yaux420, yaux421, yaux422;
	typeRNum yaux423, yaux424, yaux425, yaux426, yaux427, yaux428, yaux429, yaux43, yaux430, yaux431;
	typeRNum yaux432, yaux433, yaux435, yaux436, yaux439, yaux44, yaux440, yaux443, yaux445, yaux446;
	typeRNum yaux448, yaux45, yaux455, yaux458, yaux459, yaux461, yaux466, yaux473, yaux476, yaux48;
	typeRNum yaux481, yaux482, yaux484, yaux491, yaux492, yaux493, yaux494, yaux495, yaux496, yaux497;
	typeRNum yaux498, yaux499, yaux5, yaux50, yaux500, yaux501, yaux502, yaux503, yaux504, yaux505;
	typeRNum yaux506, yaux507, yaux508, yaux509, yaux51, yaux510, yaux511, yaux512, yaux514, yaux515;
	typeRNum yaux518, yaux519, yaux52, yaux522, yaux524, yaux525, yaux527, yaux534, yaux537, yaux538;
	typeRNum yaux54, yaux540, yaux545, yaux55, yaux552, yaux555, yaux56, yaux560, yaux561, yaux563;
	typeRNum yaux570, yaux571, yaux572, yaux573, yaux574, yaux575, yaux576, yaux577, yaux578, yaux579;
	typeRNum yaux580, yaux581, yaux582, yaux583, yaux584, yaux585, yaux586, yaux587, yaux588, yaux589;
	typeRNum yaux590, yaux591, yaux593, yaux594, yaux597, yaux598, yaux6, yaux60, yaux601, yaux603;
	typeRNum yaux604, yaux606, yaux613, yaux616, yaux617, yaux619, yaux624, yaux63, yaux631, yaux634;
	typeRNum yaux639, yaux64, yaux640, yaux642, yaux649, yaux650, yaux651, yaux652, yaux653, yaux654;
	typeRNum yaux655, yaux656, yaux657, yaux658, yaux659, yaux660, yaux661, yaux662, yaux663, yaux664;
	typeRNum yaux665, yaux666, yaux667, yaux67, yaux671, yaux675, yaux7, yaux708, yaux709, yaux710;
	typeRNum yaux711, yaux72, yaux73, yaux77, yaux8, yaux80, yaux85, yaux86, yaux89, yaux9;
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
	yaux20 = 10 + yaux19;
	yaux22 = -yaux2;
	yaux23 = yaux22 + yaux3;
	yaux24 = POW(yaux17, -1.5);
	yaux33 = POW(yaux32, -1.5);
	yaux43 = -yaux7;
	yaux44 = yaux43 + yaux8;
	yaux51 = -yaux12;
	yaux52 = yaux13 + yaux51;
	yaux28 = vec[30];
	yaux42 = -0.05500000000000001*yaux2*yaux33*yaux7;
	yaux48 = vec[34];
	yaux21 = 0.1*yaux20;
	yaux41 = vec[31];
	yaux35 = 1 / SQRT(yaux32);
	yaux36 = 0.05500000000000001*yaux35;
	yaux37 = -0.1*yaux20;
	yaux1 = vec[33];
	yaux50 = vec[35];
	yaux54 = vec[32];
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
	yaux112 = 10 + yaux111;
	yaux114 = yaux4 + yaux97;
	yaux115 = POW(yaux109, -1.5);
	yaux45 = 0.05500000000000001*yaux24*yaux44*yaux5;
	yaux123 = yaux101 + yaux9;
	yaux56 = 0.05500000000000001*yaux24*yaux5*yaux52;
	yaux130 = yaux105 + yaux14;
	yaux63 = -0.05500000000000001*yaux10*yaux24*yaux44;
	yaux64 = yaux21 + yaux63;
	yaux60 = 0.05500000000000001*yaux10*yaux23*yaux24;
	yaux127 = vec[37];
	yaux113 = 0.1*yaux112;
	yaux67 = 0.05500000000000001*yaux10*yaux24*yaux44;
	yaux119 = -0.1*yaux112;
	yaux96 = vec[36];
	yaux129 = vec[38];
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
	yaux191 = 10 + yaux190;
	yaux193 = yaux176 + yaux98;
	yaux194 = POW(yaux188, -1.5);
	yaux124 = 0.05500000000000001*yaux115*yaux123*yaux99;
	yaux202 = yaux102 + yaux180;
	yaux132 = 0.05500000000000001*yaux115*yaux130*yaux99;
	yaux209 = yaux106 + yaux184;
	yaux142 = -0.05500000000000001*yaux103*yaux115*yaux123;
	yaux143 = yaux113 + yaux142;
	yaux139 = 0.05500000000000001*yaux103*yaux114*yaux115;
	yaux206 = vec[40];
	yaux192 = 0.1*yaux191;
	yaux145 = 0.05500000000000001*yaux103*yaux115*yaux123;
	yaux198 = -0.1*yaux191;
	yaux175 = vec[39];
	yaux208 = vec[41];
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
	yaux270 = 10 + yaux269;
	yaux272 = yaux177 + yaux255;
	yaux273 = POW(yaux267, -1.5);
	yaux203 = 0.05500000000000001*yaux178*yaux194*yaux202;
	yaux281 = yaux181 + yaux259;
	yaux211 = 0.05500000000000001*yaux178*yaux194*yaux209;
	yaux288 = yaux185 + yaux263;
	yaux221 = -0.05500000000000001*yaux182*yaux194*yaux202;
	yaux222 = yaux192 + yaux221;
	yaux218 = 0.05500000000000001*yaux182*yaux193*yaux194;
	yaux285 = vec[43];
	yaux271 = 0.1*yaux270;
	yaux224 = 0.05500000000000001*yaux182*yaux194*yaux202;
	yaux277 = -0.1*yaux270;
	yaux254 = vec[42];
	yaux287 = vec[44];
	yaux229 = 0.05500000000000001*yaux182*yaux194*yaux209;
	yaux244 = -0.05500000000000001*yaux186*yaux194*yaux209;
	yaux245 = yaux192 + yaux244;
	yaux236 = 0.05500000000000001*yaux186*yaux193*yaux194;
	yaux239 = 0.05500000000000001*yaux186*yaux194*yaux202;
	yaux247 = 0.05500000000000001*yaux186*yaux194*yaux209;
	yaux274 = -0.05500000000000001*yaux257*yaux272*yaux273;
	yaux275 = yaux271 + yaux274;
	yaux334 = x[15];
	yaux335 = -yaux334;
	yaux336 = yaux255 + yaux335;
	yaux337 = yaux336 * yaux336;
	yaux338 = x[16];
	yaux339 = -yaux338;
	yaux340 = yaux259 + yaux339;
	yaux341 = yaux340 * yaux340;
	yaux342 = x[17];
	yaux343 = -yaux342;
	yaux344 = yaux263 + yaux343;
	yaux345 = yaux344 * yaux344;
	yaux346 = yaux337 + yaux341 + yaux345;
	yaux278 = 0.05500000000000001*yaux257*yaux272*yaux273;
	yaux347 = 1 / SQRT(yaux346);
	yaux348 = -0.55*yaux347;
	yaux349 = 10 + yaux348;
	yaux351 = yaux256 + yaux334;
	yaux352 = POW(yaux346, -1.5);
	yaux282 = 0.05500000000000001*yaux257*yaux273*yaux281;
	yaux360 = yaux260 + yaux338;
	yaux290 = 0.05500000000000001*yaux257*yaux273*yaux288;
	yaux367 = yaux264 + yaux342;
	yaux300 = -0.05500000000000001*yaux261*yaux273*yaux281;
	yaux301 = yaux271 + yaux300;
	yaux297 = 0.05500000000000001*yaux261*yaux272*yaux273;
	yaux364 = vec[46];
	yaux350 = 0.1*yaux349;
	yaux303 = 0.05500000000000001*yaux261*yaux273*yaux281;
	yaux356 = -0.1*yaux349;
	yaux333 = vec[45];
	yaux366 = vec[47];
	yaux308 = 0.05500000000000001*yaux261*yaux273*yaux288;
	yaux323 = -0.05500000000000001*yaux265*yaux273*yaux288;
	yaux324 = yaux271 + yaux323;
	yaux315 = 0.05500000000000001*yaux265*yaux272*yaux273;
	yaux318 = 0.05500000000000001*yaux265*yaux273*yaux281;
	yaux326 = 0.05500000000000001*yaux265*yaux273*yaux288;
	yaux353 = -0.05500000000000001*yaux336*yaux351*yaux352;
	yaux354 = yaux350 + yaux353;
	yaux413 = x[18];
	yaux414 = -yaux413;
	yaux415 = yaux334 + yaux414;
	yaux416 = yaux415 * yaux415;
	yaux417 = x[19];
	yaux418 = -yaux417;
	yaux419 = yaux338 + yaux418;
	yaux420 = yaux419 * yaux419;
	yaux421 = x[20];
	yaux422 = -yaux421;
	yaux423 = yaux342 + yaux422;
	yaux424 = yaux423 * yaux423;
	yaux425 = yaux416 + yaux420 + yaux424;
	yaux357 = 0.05500000000000001*yaux336*yaux351*yaux352;
	yaux426 = 1 / SQRT(yaux425);
	yaux427 = -0.55*yaux426;
	yaux428 = 10 + yaux427;
	yaux430 = yaux335 + yaux413;
	yaux431 = POW(yaux425, -1.5);
	yaux361 = 0.05500000000000001*yaux336*yaux352*yaux360;
	yaux439 = yaux339 + yaux417;
	yaux369 = 0.05500000000000001*yaux336*yaux352*yaux367;
	yaux446 = yaux343 + yaux421;
	yaux379 = -0.05500000000000001*yaux340*yaux352*yaux360;
	yaux380 = yaux350 + yaux379;
	yaux376 = 0.05500000000000001*yaux340*yaux351*yaux352;
	yaux443 = vec[49];
	yaux429 = 0.1*yaux428;
	yaux382 = 0.05500000000000001*yaux340*yaux352*yaux360;
	yaux435 = -0.1*yaux428;
	yaux412 = vec[48];
	yaux445 = vec[50];
	yaux387 = 0.05500000000000001*yaux340*yaux352*yaux367;
	yaux402 = -0.05500000000000001*yaux344*yaux352*yaux367;
	yaux403 = yaux350 + yaux402;
	yaux394 = 0.05500000000000001*yaux344*yaux351*yaux352;
	yaux397 = 0.05500000000000001*yaux344*yaux352*yaux360;
	yaux405 = 0.05500000000000001*yaux344*yaux352*yaux367;
	yaux432 = -0.05500000000000001*yaux415*yaux430*yaux431;
	yaux433 = yaux429 + yaux432;
	yaux492 = x[21];
	yaux493 = -yaux492;
	yaux494 = yaux413 + yaux493;
	yaux495 = yaux494 * yaux494;
	yaux496 = x[22];
	yaux497 = -yaux496;
	yaux498 = yaux417 + yaux497;
	yaux499 = yaux498 * yaux498;
	yaux500 = x[23];
	yaux501 = -yaux500;
	yaux502 = yaux421 + yaux501;
	yaux503 = yaux502 * yaux502;
	yaux504 = yaux495 + yaux499 + yaux503;
	yaux436 = 0.05500000000000001*yaux415*yaux430*yaux431;
	yaux505 = 1 / SQRT(yaux504);
	yaux506 = -0.55*yaux505;
	yaux507 = 10 + yaux506;
	yaux509 = yaux414 + yaux492;
	yaux510 = POW(yaux504, -1.5);
	yaux440 = 0.05500000000000001*yaux415*yaux431*yaux439;
	yaux518 = yaux418 + yaux496;
	yaux448 = 0.05500000000000001*yaux415*yaux431*yaux446;
	yaux525 = yaux422 + yaux500;
	yaux458 = -0.05500000000000001*yaux419*yaux431*yaux439;
	yaux459 = yaux429 + yaux458;
	yaux455 = 0.05500000000000001*yaux419*yaux430*yaux431;
	yaux522 = vec[52];
	yaux508 = 0.1*yaux507;
	yaux461 = 0.05500000000000001*yaux419*yaux431*yaux439;
	yaux514 = -0.1*yaux507;
	yaux491 = vec[51];
	yaux524 = vec[53];
	yaux466 = 0.05500000000000001*yaux419*yaux431*yaux446;
	yaux481 = -0.05500000000000001*yaux423*yaux431*yaux446;
	yaux482 = yaux429 + yaux481;
	yaux473 = 0.05500000000000001*yaux423*yaux430*yaux431;
	yaux476 = 0.05500000000000001*yaux423*yaux431*yaux439;
	yaux484 = 0.05500000000000001*yaux423*yaux431*yaux446;
	yaux511 = -0.05500000000000001*yaux494*yaux509*yaux510;
	yaux512 = yaux508 + yaux511;
	yaux571 = x[24];
	yaux572 = -yaux571;
	yaux573 = yaux492 + yaux572;
	yaux574 = yaux573 * yaux573;
	yaux575 = x[25];
	yaux576 = -yaux575;
	yaux577 = yaux496 + yaux576;
	yaux578 = yaux577 * yaux577;
	yaux579 = x[26];
	yaux580 = -yaux579;
	yaux581 = yaux500 + yaux580;
	yaux582 = yaux581 * yaux581;
	yaux583 = yaux574 + yaux578 + yaux582;
	yaux515 = 0.05500000000000001*yaux494*yaux509*yaux510;
	yaux584 = 1 / SQRT(yaux583);
	yaux585 = -0.55*yaux584;
	yaux586 = 10 + yaux585;
	yaux588 = yaux493 + yaux571;
	yaux589 = POW(yaux583, -1.5);
	yaux519 = 0.05500000000000001*yaux494*yaux510*yaux518;
	yaux597 = yaux497 + yaux575;
	yaux527 = 0.05500000000000001*yaux494*yaux510*yaux525;
	yaux604 = yaux501 + yaux579;
	yaux537 = -0.05500000000000001*yaux498*yaux510*yaux518;
	yaux538 = yaux508 + yaux537;
	yaux534 = 0.05500000000000001*yaux498*yaux509*yaux510;
	yaux601 = vec[55];
	yaux587 = 0.1*yaux586;
	yaux540 = 0.05500000000000001*yaux498*yaux510*yaux518;
	yaux593 = -0.1*yaux586;
	yaux570 = vec[54];
	yaux603 = vec[56];
	yaux545 = 0.05500000000000001*yaux498*yaux510*yaux525;
	yaux560 = -0.05500000000000001*yaux502*yaux510*yaux525;
	yaux561 = yaux508 + yaux560;
	yaux552 = 0.05500000000000001*yaux502*yaux509*yaux510;
	yaux555 = 0.05500000000000001*yaux502*yaux510*yaux518;
	yaux563 = 0.05500000000000001*yaux502*yaux510*yaux525;
	yaux590 = -0.05500000000000001*yaux573*yaux588*yaux589;
	yaux591 = yaux587 + yaux590;
	yaux594 = 0.05500000000000001*yaux573*yaux588*yaux589;
	yaux649 = x[27];
	yaux650 = -yaux649;
	yaux651 = yaux571 + yaux650;
	yaux652 = yaux651 * yaux651;
	yaux653 = x[28];
	yaux654 = -yaux653;
	yaux655 = yaux575 + yaux654;
	yaux656 = yaux655 * yaux655;
	yaux657 = x[29];
	yaux658 = -yaux657;
	yaux659 = yaux579 + yaux658;
	yaux660 = yaux659 * yaux659;
	yaux661 = yaux652 + yaux656 + yaux660;
	yaux598 = 0.05500000000000001*yaux573*yaux589*yaux597;
	yaux667 = POW(yaux661, -1.5);
	yaux606 = 0.05500000000000001*yaux573*yaux589*yaux604;
	yaux616 = -0.05500000000000001*yaux577*yaux589*yaux597;
	yaux617 = yaux587 + yaux616;
	yaux613 = 0.05500000000000001*yaux577*yaux588*yaux589;
	yaux666 = yaux572 + yaux649;
	yaux619 = 0.05500000000000001*yaux577*yaux589*yaux597;
	yaux662 = 1 / SQRT(yaux661);
	yaux663 = -0.55*yaux662;
	yaux664 = 10 + yaux663;
	yaux665 = -0.1*yaux664;
	yaux671 = yaux576 + yaux653;
	yaux624 = 0.05500000000000001*yaux577*yaux589*yaux604;
	yaux675 = yaux580 + yaux657;
	yaux639 = -0.05500000000000001*yaux581*yaux589*yaux604;
	yaux640 = yaux587 + yaux639;
	yaux631 = 0.05500000000000001*yaux581*yaux588*yaux589;
	yaux634 = 0.05500000000000001*yaux581*yaux589*yaux597;
	yaux642 = 0.05500000000000001*yaux581*yaux589*yaux604;
	yaux708 = SQRT(yaux661);
	yaux709 = -yaux708;
	yaux710 = 0.05500000000000001 + yaux709;
	yaux711 = -18.18181818181818*yaux661*yaux710;

	out[0] = 22.22222222222222*yaux1*yaux26 + 22.22222222222222*yaux28*(-1. - 0.05500000000000001*yaux29*yaux33 + yaux36 + yaux37 + yaux38) + 22.22222222222222*yaux41*(yaux42 + yaux45) - 1.2222222222222223*yaux24*yaux44*yaux48*yaux5 - 1.2222222222222223*yaux24*yaux5*yaux50*yaux52 + 22.22222222222222*yaux54*(yaux55 + yaux56);
	out[1] = -1.2222222222222223*yaux1*yaux10*yaux23*yaux24 - 1.2222222222222223*yaux10*yaux24*yaux50*yaux52 + 22.22222222222222*yaux28*(yaux42 + yaux60) + 22.22222222222222*yaux48*yaux64 + 22.22222222222222*yaux41*(-1. - 0.05500000000000001*yaux30*yaux33 + yaux36 + yaux37 + yaux67) + 22.22222222222222*yaux54*(yaux72 + yaux73);
	out[2] = -1.2222222222222223*yaux1*yaux15*yaux23*yaux24 - 1.2222222222222223*yaux15*yaux24*yaux44*yaux48 + 22.22222222222222*yaux28*(yaux55 + yaux77) + 22.22222222222222*yaux41*(yaux72 + yaux80) + 22.22222222222222*yaux50*yaux86 + 22.22222222222222*yaux54*(-1. - 0.05500000000000001*yaux31*yaux33 + yaux36 + yaux37 + yaux89);
	out[3] = 22.22222222222222*yaux26*yaux28 + 22.22222222222222*yaux1*(yaux119 + yaux120 + yaux37 + yaux38) + 22.22222222222222*(yaux124 + yaux45)*yaux48 - 1.2222222222222223*yaux24*yaux41*yaux44*yaux5 - 1.2222222222222223*yaux24*yaux5*yaux52*yaux54 + 22.22222222222222*yaux50*(yaux132 + yaux56) + 22.22222222222222*yaux117*yaux96 - 1.2222222222222223*yaux115*yaux123*yaux127*yaux99 - 1.2222222222222223*yaux115*yaux129*yaux130*yaux99;
	out[4] = -1.2222222222222223*yaux103*yaux115*yaux129*yaux130 + 22.22222222222222*yaux127*yaux143 - 1.2222222222222223*yaux10*yaux23*yaux24*yaux28 - 1.2222222222222223*yaux10*yaux24*yaux52*yaux54 + 22.22222222222222*yaux1*(yaux139 + yaux60) + 22.22222222222222*yaux41*yaux64 + 22.22222222222222*yaux48*(yaux119 + yaux145 + yaux37 + yaux67) + 22.22222222222222*yaux50*(yaux150 + yaux73) - 1.2222222222222223*yaux103*yaux114*yaux115*yaux96;
	out[5] = -1.2222222222222223*yaux107*yaux115*yaux123*yaux127 + 22.22222222222222*yaux129*yaux166 - 1.2222222222222223*yaux15*yaux23*yaux24*yaux28 - 1.2222222222222223*yaux15*yaux24*yaux41*yaux44 + 22.22222222222222*yaux1*(yaux157 + yaux77) + 22.22222222222222*yaux48*(yaux160 + yaux80) + 22.22222222222222*yaux54*yaux86 + 22.22222222222222*yaux50*(yaux119 + yaux168 + yaux37 + yaux89) - 1.2222222222222223*yaux107*yaux114*yaux115*yaux96;
	out[6] = 22.22222222222222*yaux1*yaux117 + 22.22222222222222*yaux175*yaux196 + 22.22222222222222*yaux127*(yaux124 + yaux203) - 1.2222222222222223*yaux178*yaux194*yaux202*yaux206 - 1.2222222222222223*yaux178*yaux194*yaux208*yaux209 + 22.22222222222222*yaux129*(yaux132 + yaux211) + 22.22222222222222*(yaux119 + yaux120 + yaux198 + yaux199)*yaux96 - 1.2222222222222223*yaux115*yaux123*yaux48*yaux99 - 1.2222222222222223*yaux115*yaux130*yaux50*yaux99;
	out[7] = -1.2222222222222223*yaux1*yaux103*yaux114*yaux115 - 1.2222222222222223*yaux175*yaux182*yaux193*yaux194 - 1.2222222222222223*yaux182*yaux194*yaux208*yaux209 + 22.22222222222222*yaux206*yaux222 + 22.22222222222222*yaux127*(yaux119 + yaux145 + yaux198 + yaux224) + 22.22222222222222*yaux129*(yaux150 + yaux229) + 22.22222222222222*yaux143*yaux48 - 1.2222222222222223*yaux103*yaux115*yaux130*yaux50 + 22.22222222222222*(yaux139 + yaux218)*yaux96;
	out[8] = -1.2222222222222223*yaux1*yaux107*yaux114*yaux115 - 1.2222222222222223*yaux175*yaux186*yaux193*yaux194 - 1.2222222222222223*yaux186*yaux194*yaux202*yaux206 + 22.22222222222222*yaux127*(yaux160 + yaux239) + 22.22222222222222*yaux208*yaux245 + 22.22222222222222*yaux129*(yaux119 + yaux168 + yaux198 + yaux247) - 1.2222222222222223*yaux107*yaux115*yaux123*yaux48 + 22.22222222222222*yaux166*yaux50 + 22.22222222222222*(yaux157 + yaux236)*yaux96;
	out[9] = -1.2222222222222223*yaux127*yaux178*yaux194*yaux202 - 1.2222222222222223*yaux129*yaux178*yaux194*yaux209 + 22.22222222222222*yaux254*yaux275 + 22.22222222222222*yaux175*(yaux198 + yaux199 + yaux277 + yaux278) + 22.22222222222222*yaux206*(yaux203 + yaux282) - 1.2222222222222223*yaux257*yaux273*yaux281*yaux285 - 1.2222222222222223*yaux257*yaux273*yaux287*yaux288 + 22.22222222222222*yaux208*(yaux211 + yaux290) + 22.22222222222222*yaux196*yaux96;
	out[10] = -1.2222222222222223*yaux129*yaux182*yaux194*yaux209 + 22.22222222222222*yaux127*yaux222 - 1.2222222222222223*yaux254*yaux261*yaux272*yaux273 - 1.2222222222222223*yaux261*yaux273*yaux287*yaux288 + 22.22222222222222*yaux175*(yaux218 + yaux297) + 22.22222222222222*yaux285*yaux301 + 22.22222222222222*yaux206*(yaux198 + yaux224 + yaux277 + yaux303) + 22.22222222222222*yaux208*(yaux229 + yaux308) - 1.2222222222222223*yaux182*yaux193*yaux194*yaux96;
	out[11] = -1.2222222222222223*yaux127*yaux186*yaux194*yaux202 + 22.22222222222222*yaux129*yaux245 - 1.2222222222222223*yaux254*yaux265*yaux272*yaux273 - 1.2222222222222223*yaux265*yaux273*yaux281*yaux285 + 22.22222222222222*yaux175*(yaux236 + yaux315) + 22.22222222222222*yaux206*(yaux239 + yaux318) + 22.22222222222222*yaux287*yaux324 + 22.22222222222222*yaux208*(yaux198 + yaux247 + yaux277 + yaux326) - 1.2222222222222223*yaux186*yaux193*yaux194*yaux96;
	out[12] = 22.22222222222222*yaux175*yaux275 - 1.2222222222222223*yaux206*yaux257*yaux273*yaux281 - 1.2222222222222223*yaux208*yaux257*yaux273*yaux288 + 22.22222222222222*yaux333*yaux354 + 22.22222222222222*yaux254*(yaux277 + yaux278 + yaux356 + yaux357) + 22.22222222222222*yaux285*(yaux282 + yaux361) - 1.2222222222222223*yaux336*yaux352*yaux360*yaux364 - 1.2222222222222223*yaux336*yaux352*yaux366*yaux367 + 22.22222222222222*yaux287*(yaux290 + yaux369);
	out[13] = -1.2222222222222223*yaux175*yaux261*yaux272*yaux273 - 1.2222222222222223*yaux208*yaux261*yaux273*yaux288 + 22.22222222222222*yaux206*yaux301 - 1.2222222222222223*yaux333*yaux340*yaux351*yaux352 - 1.2222222222222223*yaux340*yaux352*yaux366*yaux367 + 22.22222222222222*yaux254*(yaux297 + yaux376) + 22.22222222222222*yaux364*yaux380 + 22.22222222222222*yaux285*(yaux277 + yaux303 + yaux356 + yaux382) + 22.22222222222222*yaux287*(yaux308 + yaux387);
	out[14] = -1.2222222222222223*yaux175*yaux265*yaux272*yaux273 - 1.2222222222222223*yaux206*yaux265*yaux273*yaux281 + 22.22222222222222*yaux208*yaux324 - 1.2222222222222223*yaux333*yaux344*yaux351*yaux352 - 1.2222222222222223*yaux344*yaux352*yaux360*yaux364 + 22.22222222222222*yaux254*(yaux315 + yaux394) + 22.22222222222222*yaux285*(yaux318 + yaux397) + 22.22222222222222*yaux366*yaux403 + 22.22222222222222*yaux287*(yaux277 + yaux326 + yaux356 + yaux405);
	out[15] = 22.22222222222222*yaux254*yaux354 - 1.2222222222222223*yaux285*yaux336*yaux352*yaux360 - 1.2222222222222223*yaux287*yaux336*yaux352*yaux367 + 22.22222222222222*yaux412*yaux433 + 22.22222222222222*yaux333*(yaux356 + yaux357 + yaux435 + yaux436) + 22.22222222222222*yaux364*(yaux361 + yaux440) - 1.2222222222222223*yaux415*yaux431*yaux439*yaux443 - 1.2222222222222223*yaux415*yaux431*yaux445*yaux446 + 22.22222222222222*yaux366*(yaux369 + yaux448);
	out[16] = -1.2222222222222223*yaux254*yaux340*yaux351*yaux352 - 1.2222222222222223*yaux287*yaux340*yaux352*yaux367 + 22.22222222222222*yaux285*yaux380 - 1.2222222222222223*yaux412*yaux419*yaux430*yaux431 - 1.2222222222222223*yaux419*yaux431*yaux445*yaux446 + 22.22222222222222*yaux333*(yaux376 + yaux455) + 22.22222222222222*yaux443*yaux459 + 22.22222222222222*yaux364*(yaux356 + yaux382 + yaux435 + yaux461) + 22.22222222222222*yaux366*(yaux387 + yaux466);
	out[17] = -1.2222222222222223*yaux254*yaux344*yaux351*yaux352 - 1.2222222222222223*yaux285*yaux344*yaux352*yaux360 + 22.22222222222222*yaux287*yaux403 - 1.2222222222222223*yaux412*yaux423*yaux430*yaux431 - 1.2222222222222223*yaux423*yaux431*yaux439*yaux443 + 22.22222222222222*yaux333*(yaux394 + yaux473) + 22.22222222222222*yaux364*(yaux397 + yaux476) + 22.22222222222222*yaux445*yaux482 + 22.22222222222222*yaux366*(yaux356 + yaux405 + yaux435 + yaux484);
	out[18] = 22.22222222222222*yaux333*yaux433 - 1.2222222222222223*yaux364*yaux415*yaux431*yaux439 - 1.2222222222222223*yaux366*yaux415*yaux431*yaux446 + 22.22222222222222*yaux491*yaux512 + 22.22222222222222*yaux412*(yaux435 + yaux436 + yaux514 + yaux515) + 22.22222222222222*yaux443*(yaux440 + yaux519) - 1.2222222222222223*yaux494*yaux510*yaux518*yaux522 - 1.2222222222222223*yaux494*yaux510*yaux524*yaux525 + 22.22222222222222*yaux445*(yaux448 + yaux527);
	out[19] = -1.2222222222222223*yaux333*yaux419*yaux430*yaux431 - 1.2222222222222223*yaux366*yaux419*yaux431*yaux446 + 22.22222222222222*yaux364*yaux459 - 1.2222222222222223*yaux491*yaux498*yaux509*yaux510 - 1.2222222222222223*yaux498*yaux510*yaux524*yaux525 + 22.22222222222222*yaux412*(yaux455 + yaux534) + 22.22222222222222*yaux522*yaux538 + 22.22222222222222*yaux443*(yaux435 + yaux461 + yaux514 + yaux540) + 22.22222222222222*yaux445*(yaux466 + yaux545);
	out[20] = -1.2222222222222223*yaux333*yaux423*yaux430*yaux431 - 1.2222222222222223*yaux364*yaux423*yaux431*yaux439 + 22.22222222222222*yaux366*yaux482 - 1.2222222222222223*yaux491*yaux502*yaux509*yaux510 - 1.2222222222222223*yaux502*yaux510*yaux518*yaux522 + 22.22222222222222*yaux412*(yaux473 + yaux552) + 22.22222222222222*yaux443*(yaux476 + yaux555) + 22.22222222222222*yaux524*yaux561 + 22.22222222222222*yaux445*(yaux435 + yaux484 + yaux514 + yaux563);
	out[21] = 22.22222222222222*yaux412*yaux512 - 1.2222222222222223*yaux443*yaux494*yaux510*yaux518 - 1.2222222222222223*yaux445*yaux494*yaux510*yaux525 + 22.22222222222222*yaux570*yaux591 + 22.22222222222222*yaux491*(yaux514 + yaux515 + yaux593 + yaux594) + 22.22222222222222*yaux522*(yaux519 + yaux598) - 1.2222222222222223*yaux573*yaux589*yaux597*yaux601 - 1.2222222222222223*yaux573*yaux589*yaux603*yaux604 + 22.22222222222222*yaux524*(yaux527 + yaux606);
	out[22] = -1.2222222222222223*yaux412*yaux498*yaux509*yaux510 - 1.2222222222222223*yaux445*yaux498*yaux510*yaux525 + 22.22222222222222*yaux443*yaux538 - 1.2222222222222223*yaux570*yaux577*yaux588*yaux589 - 1.2222222222222223*yaux577*yaux589*yaux603*yaux604 + 22.22222222222222*yaux491*(yaux534 + yaux613) + 22.22222222222222*yaux601*yaux617 + 22.22222222222222*yaux522*(yaux514 + yaux540 + yaux593 + yaux619) + 22.22222222222222*yaux524*(yaux545 + yaux624);
	out[23] = -1.2222222222222223*yaux412*yaux502*yaux509*yaux510 - 1.2222222222222223*yaux443*yaux502*yaux510*yaux518 + 22.22222222222222*yaux445*yaux561 - 1.2222222222222223*yaux570*yaux581*yaux588*yaux589 - 1.2222222222222223*yaux581*yaux589*yaux597*yaux601 + 22.22222222222222*yaux491*(yaux552 + yaux631) + 22.22222222222222*yaux522*(yaux555 + yaux634) + 22.22222222222222*yaux603*yaux640 + 22.22222222222222*yaux524*(yaux514 + yaux563 + yaux593 + yaux642);
	out[24] = 22.22222222222222*yaux491*yaux591 - 1.2222222222222223*yaux522*yaux573*yaux589*yaux597 - 1.2222222222222223*yaux524*yaux573*yaux589*yaux604 + 22.22222222222222*yaux570*(yaux593 + yaux594 + yaux665 + 0.05500000000000001*yaux651*yaux666*yaux667) + 22.22222222222222*yaux601*(yaux598 + 0.05500000000000001*yaux651*yaux667*yaux671) + 22.22222222222222*yaux603*(yaux606 + 0.05500000000000001*yaux651*yaux667*yaux675);
	out[25] = -1.2222222222222223*yaux491*yaux577*yaux588*yaux589 - 1.2222222222222223*yaux524*yaux577*yaux589*yaux604 + 22.22222222222222*yaux522*yaux617 + 22.22222222222222*yaux570*(yaux613 + 0.05500000000000001*yaux655*yaux666*yaux667) + 22.22222222222222*yaux601*(yaux593 + yaux619 + yaux665 + 0.05500000000000001*yaux655*yaux667*yaux671) + 22.22222222222222*yaux603*(yaux624 + 0.05500000000000001*yaux655*yaux667*yaux675);
	out[26] = -1.2222222222222223*yaux491*yaux581*yaux588*yaux589 - 1.2222222222222223*yaux522*yaux581*yaux589*yaux597 + 22.22222222222222*yaux524*yaux640 + 22.22222222222222*yaux570*(yaux631 + 0.05500000000000001*yaux659*yaux666*yaux667) + 22.22222222222222*yaux601*(yaux634 + 0.05500000000000001*yaux659*yaux667*yaux671) + 22.22222222222222*yaux603*(yaux593 + yaux642 + yaux665 + 0.05500000000000001*yaux659*yaux667*yaux675);
	out[27] = 1.2222222222222223*yaux667*(1 * yaux601*yaux651*yaux655 + 1 * yaux603*yaux651*yaux659 + 1 * yaux570*(1 * yaux652 + yaux711));
	out[28] = 1.2222222222222223*yaux667*(1 * yaux570*yaux651*yaux655 + 1 * yaux603*yaux655*yaux659 + 1 * yaux601*(1 * yaux656 + yaux711));
	out[29] = -22.22222222222222*yaux667*(-0.05500000000000001*yaux570*yaux651*yaux659 - 0.05500000000000001*yaux601*yaux655*yaux659 + 1 * yaux603*(-0.05500000000000001*yaux660 + 1 * yaux661*yaux710));
	out[30] = vec[0];
	out[31] = vec[1];
	out[32] = vec[2];
	out[33] = vec[3];
	out[34] = vec[4];
	out[35] = vec[5];
	out[36] = vec[6];
	out[37] = vec[7];
	out[38] = vec[8];
	out[39] = vec[9];
	out[40] = vec[10];
	out[41] = vec[11];
	out[42] = vec[12];
	out[43] = vec[13];
	out[44] = vec[14];
	out[45] = vec[15];
	out[46] = vec[16];
	out[47] = vec[17];
	out[48] = vec[18];
	out[49] = vec[19];
	out[50] = vec[20];
	out[51] = vec[21];
	out[52] = vec[22];
	out[53] = vec[23];
	out[54] = vec[24];
	out[55] = vec[25];
	out[56] = vec[26];
}
/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	out[0] = vec[27];
	out[1] = vec[28];
	out[2] = vec[29];
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

	yaux11 = x[27];


	out[0] = param[2] * (POW2(u[0]) + POW2(u[1]) + POW2(u[2]))
		+ param[1] * (56.25 - 15 * yaux11 + yaux11 * yaux11 + POW2(x[28]) + POW2(x[29]))
		+ param[0] * (POW2(x[30]) + POW2(x[31]) + POW2(x[32]) + POW2(x[33]) + POW2(x[34]) + POW2(x[35])
			+ POW2(x[36]) + POW2(x[37]) + POW2(x[38]) + POW2(x[39]) + POW2(x[40]) + POW2(x[41]) + POW2(x[42])
			+ POW2(x[43]) + POW2(x[44]) + POW2(x[45]) + POW2(x[46]) + POW2(x[47]) + POW2(x[48]) + POW2(x[49])
			+ POW2(x[50]) + POW2(x[51]) + POW2(x[52]) + POW2(x[53]) + POW2(x[54]) + POW2(x[55]) + POW2(x[56]));
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
	out[15] = 0;
	out[16] = 0;
	out[17] = 0;
	out[18] = 0;
	out[19] = 0;
	out[20] = 0;
	out[21] = 0;
	out[22] = 0;
	out[23] = 0;
	out[24] = 0;
	out[25] = 0;
	out[26] = 0;
	out[27] = 2 * yaux1*(-7.5 + x[27]);
	out[28] = 2 * yaux1*x[28];
	out[29] = 2 * yaux1*x[29];
	out[30] = 2 * yaux9*x[30];
	out[31] = 2 * yaux9*x[31];
	out[32] = 2 * yaux9*x[32];
	out[33] = 2 * yaux9*x[33];
	out[34] = 2 * yaux9*x[34];
	out[35] = 2 * yaux9*x[35];
	out[36] = 2 * yaux9*x[36];
	out[37] = 2 * yaux9*x[37];
	out[38] = 2 * yaux9*x[38];
	out[39] = 2 * yaux9*x[39];
	out[40] = 2 * yaux9*x[40];
	out[41] = 2 * yaux9*x[41];
	out[42] = 2 * yaux9*x[42];
	out[43] = 2 * yaux9*x[43];
	out[44] = 2 * yaux9*x[44];
	out[45] = 2 * yaux9*x[45];
	out[46] = 2 * yaux9*x[46];
	out[47] = 2 * yaux9*x[47];
	out[48] = 2 * yaux9*x[48];
	out[49] = 2 * yaux9*x[49];
	out[50] = 2 * yaux9*x[50];
	out[51] = 2 * yaux9*x[51];
	out[52] = 2 * yaux9*x[52];
	out[53] = 2 * yaux9*x[53];
	out[54] = 2 * yaux9*x[54];
	out[55] = 2 * yaux9*x[55];
	out[56] = 2 * yaux9*x[56];
}
/** Gradient dl/du **/
void dldu(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *u, ctypeRNum *p, ctypeRNum *xdes, ctypeRNum *udes, typeUSERPARAM *userparam)
{
	typeRNum* param = (typeRNum*)userparam;

	typeRNum yaux1 = param[2];

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

	out[0] = param[3] * (POW2(-7.5 + x[27]) + POW2(x[28]) + POW2(x[29]));
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
	out[15] = 0;
	out[16] = 0;
	out[17] = 0;
	out[18] = 0;
	out[19] = 0;
	out[20] = 0;
	out[21] = 0;
	out[22] = 0;
	out[23] = 0;
	out[24] = 0;
	out[25] = 0;
	out[26] = 0;
	out[27] = 2 * yaux1*(-7.5 + x[27]);
	out[28] = 2 * yaux1*x[28];
	out[29] = 2 * yaux1*x[29];
	out[30] = 0;
	out[31] = 0;
	out[32] = 0;
	out[33] = 0;
	out[34] = 0;
	out[35] = 0;
	out[36] = 0;
	out[37] = 0;
	out[38] = 0;
	out[39] = 0;
	out[40] = 0;
	out[41] = 0;
	out[42] = 0;
	out[43] = 0;
	out[44] = 0;
	out[45] = 0;
	out[46] = 0;
	out[47] = 0;
	out[48] = 0;
	out[49] = 0;
	out[50] = 0;
	out[51] = 0;
	out[52] = 0;
	out[53] = 0;
	out[54] = 0;
	out[55] = 0;
	out[56] = 0;
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
