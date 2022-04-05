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
	*Nx = 69;
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
	typeRNum yaux104, yaux109, yaux115, yaux116, yaux117, yaux118, yaux119, yaux120, yaux121, yaux122;
	typeRNum yaux123, yaux124, yaux125, yaux126, yaux127, yaux128, yaux129, yaux130, yaux131, yaux136;
	typeRNum yaux141, yaux147, yaux148, yaux149, yaux150, yaux151, yaux152, yaux153, yaux154, yaux155;
	typeRNum yaux156, yaux157, yaux158, yaux159, yaux160, yaux161, yaux162, yaux163, yaux168, yaux173;
	typeRNum yaux179, yaux180, yaux181, yaux182, yaux183, yaux184, yaux185, yaux186, yaux187, yaux188;
	typeRNum yaux189, yaux190, yaux191, yaux192, yaux193, yaux194, yaux195, yaux200, yaux205, yaux211;
	typeRNum yaux212, yaux213, yaux214, yaux215, yaux216, yaux217, yaux218, yaux219, yaux220, yaux221;
	typeRNum yaux222, yaux223, yaux224, yaux225, yaux226, yaux227, yaux232, yaux237, yaux243, yaux244;
	typeRNum yaux245, yaux246, yaux247, yaux248, yaux249, yaux250, yaux251, yaux252, yaux253, yaux254;
	typeRNum yaux255, yaux256, yaux257, yaux258, yaux259, yaux264, yaux269, yaux275, yaux276, yaux277;
	typeRNum yaux278, yaux279, yaux280, yaux281, yaux282, yaux283, yaux284, yaux285, yaux286, yaux287;
	typeRNum yaux288, yaux289, yaux290, yaux291, yaux296, yaux301, yaux307, yaux308, yaux309, yaux310;
	typeRNum yaux311, yaux312, yaux313, yaux314, yaux315, yaux316, yaux317, yaux318, yaux319, yaux320;
	typeRNum yaux321, yaux322, yaux323, yaux328, yaux333, yaux339, yaux340, yaux341, yaux342, yaux343;
	typeRNum yaux344, yaux345, yaux346, yaux347, yaux348, yaux349, yaux350, yaux351, yaux352, yaux353;
	typeRNum yaux354, yaux355, yaux360, yaux365, yaux37, yaux371, yaux373, yaux374, yaux375, yaux376;
	typeRNum yaux377, yaux378, yaux379, yaux38, yaux380, yaux381, yaux382, yaux383, yaux384, yaux385;
	typeRNum yaux386, yaux387, yaux39, yaux40, yaux41, yaux42, yaux43, yaux44, yaux45, yaux46;
	typeRNum yaux48, yaux49, yaux50, yaux51, yaux52, yaux53, yaux54, yaux55, yaux56, yaux57;
	typeRNum yaux58, yaux59, yaux60, yaux61, yaux62, yaux63, yaux64, yaux65, yaux70, yaux71;
	typeRNum yaux76, yaux77, yaux83, yaux84, yaux85, yaux86, yaux87, yaux88, yaux89, yaux90;
	typeRNum yaux91, yaux92, yaux93, yaux94, yaux95, yaux96, yaux97, yaux98, yaux99;

	yaux37 = x[0];
	yaux49 = x[3];
	yaux39 = x[1];
	yaux41 = x[2];
	yaux38 = yaux37 * yaux37;
	yaux40 = yaux39 * yaux39;
	yaux42 = yaux41 * yaux41;
	yaux43 = yaux38 + yaux40 + yaux42;
	yaux44 = 1 / SQRT(yaux43);
	yaux45 = -0.55*yaux44;
	yaux46 = 12 + yaux45;
	yaux54 = x[4];
	yaux51 = -yaux49;
	yaux52 = yaux37 + yaux51;
	yaux53 = yaux52 * yaux52;
	yaux55 = -yaux54;
	yaux56 = yaux39 + yaux55;
	yaux57 = yaux56 * yaux56;
	yaux58 = x[5];
	yaux59 = -yaux58;
	yaux60 = yaux41 + yaux59;
	yaux61 = yaux60 * yaux60;
	yaux62 = yaux53 + yaux57 + yaux61;
	yaux63 = 1 / SQRT(yaux62);
	yaux64 = -0.55*yaux63;
	yaux65 = 12 + yaux64;
	yaux48 = -yaux37;
	yaux50 = yaux48 + yaux49;
	yaux83 = x[6];
	yaux70 = -yaux39;
	yaux71 = yaux54 + yaux70;
	yaux88 = x[7];
	yaux85 = -yaux83;
	yaux86 = yaux49 + yaux85;
	yaux87 = yaux86 * yaux86;
	yaux89 = -yaux88;
	yaux90 = yaux54 + yaux89;
	yaux91 = yaux90 * yaux90;
	yaux92 = x[8];
	yaux93 = -yaux92;
	yaux94 = yaux58 + yaux93;
	yaux95 = yaux94 * yaux94;
	yaux96 = yaux87 + yaux91 + yaux95;
	yaux97 = 1 / SQRT(yaux96);
	yaux98 = -0.55*yaux97;
	yaux99 = 12 + yaux98;
	yaux76 = -yaux41;
	yaux77 = yaux58 + yaux76;
	yaux84 = yaux51 + yaux83;
	yaux115 = x[9];
	yaux104 = yaux55 + yaux88;
	yaux120 = x[10];
	yaux117 = -yaux115;
	yaux118 = yaux117 + yaux83;
	yaux119 = yaux118 * yaux118;
	yaux121 = -yaux120;
	yaux122 = yaux121 + yaux88;
	yaux123 = yaux122 * yaux122;
	yaux124 = x[11];
	yaux125 = -yaux124;
	yaux126 = yaux125 + yaux92;
	yaux127 = yaux126 * yaux126;
	yaux128 = yaux119 + yaux123 + yaux127;
	yaux129 = 1 / SQRT(yaux128);
	yaux130 = -0.55*yaux129;
	yaux131 = 12 + yaux130;
	yaux109 = yaux59 + yaux92;
	yaux116 = yaux115 + yaux85;
	yaux147 = x[12];
	yaux136 = yaux120 + yaux89;
	yaux152 = x[13];
	yaux149 = -yaux147;
	yaux150 = yaux115 + yaux149;
	yaux151 = yaux150 * yaux150;
	yaux153 = -yaux152;
	yaux154 = yaux120 + yaux153;
	yaux155 = yaux154 * yaux154;
	yaux156 = x[14];
	yaux157 = -yaux156;
	yaux158 = yaux124 + yaux157;
	yaux159 = yaux158 * yaux158;
	yaux160 = yaux151 + yaux155 + yaux159;
	yaux161 = 1 / SQRT(yaux160);
	yaux162 = -0.55*yaux161;
	yaux163 = 12 + yaux162;
	yaux141 = yaux124 + yaux93;
	yaux148 = yaux117 + yaux147;
	yaux179 = x[15];
	yaux168 = yaux121 + yaux152;
	yaux184 = x[16];
	yaux181 = -yaux179;
	yaux182 = yaux147 + yaux181;
	yaux183 = yaux182 * yaux182;
	yaux185 = -yaux184;
	yaux186 = yaux152 + yaux185;
	yaux187 = yaux186 * yaux186;
	yaux188 = x[17];
	yaux189 = -yaux188;
	yaux190 = yaux156 + yaux189;
	yaux191 = yaux190 * yaux190;
	yaux192 = yaux183 + yaux187 + yaux191;
	yaux193 = 1 / SQRT(yaux192);
	yaux194 = -0.55*yaux193;
	yaux195 = 12 + yaux194;
	yaux173 = yaux125 + yaux156;
	yaux180 = yaux149 + yaux179;
	yaux211 = x[18];
	yaux200 = yaux153 + yaux184;
	yaux216 = x[19];
	yaux213 = -yaux211;
	yaux214 = yaux179 + yaux213;
	yaux215 = yaux214 * yaux214;
	yaux217 = -yaux216;
	yaux218 = yaux184 + yaux217;
	yaux219 = yaux218 * yaux218;
	yaux220 = x[20];
	yaux221 = -yaux220;
	yaux222 = yaux188 + yaux221;
	yaux223 = yaux222 * yaux222;
	yaux224 = yaux215 + yaux219 + yaux223;
	yaux225 = 1 / SQRT(yaux224);
	yaux226 = -0.55*yaux225;
	yaux227 = 12 + yaux226;
	yaux205 = yaux157 + yaux188;
	yaux212 = yaux181 + yaux211;
	yaux243 = x[21];
	yaux232 = yaux185 + yaux216;
	yaux248 = x[22];
	yaux245 = -yaux243;
	yaux246 = yaux211 + yaux245;
	yaux247 = yaux246 * yaux246;
	yaux249 = -yaux248;
	yaux250 = yaux216 + yaux249;
	yaux251 = yaux250 * yaux250;
	yaux252 = x[23];
	yaux253 = -yaux252;
	yaux254 = yaux220 + yaux253;
	yaux255 = yaux254 * yaux254;
	yaux256 = yaux247 + yaux251 + yaux255;
	yaux257 = 1 / SQRT(yaux256);
	yaux258 = -0.55*yaux257;
	yaux259 = 12 + yaux258;
	yaux237 = yaux189 + yaux220;
	yaux244 = yaux213 + yaux243;
	yaux275 = x[24];
	yaux264 = yaux217 + yaux248;
	yaux280 = x[25];
	yaux277 = -yaux275;
	yaux278 = yaux243 + yaux277;
	yaux279 = yaux278 * yaux278;
	yaux281 = -yaux280;
	yaux282 = yaux248 + yaux281;
	yaux283 = yaux282 * yaux282;
	yaux284 = x[26];
	yaux285 = -yaux284;
	yaux286 = yaux252 + yaux285;
	yaux287 = yaux286 * yaux286;
	yaux288 = yaux279 + yaux283 + yaux287;
	yaux289 = 1 / SQRT(yaux288);
	yaux290 = -0.55*yaux289;
	yaux291 = 12 + yaux290;
	yaux269 = yaux221 + yaux252;
	yaux276 = yaux245 + yaux275;
	yaux307 = x[27];
	yaux296 = yaux249 + yaux280;
	yaux312 = x[28];
	yaux309 = -yaux307;
	yaux310 = yaux275 + yaux309;
	yaux311 = yaux310 * yaux310;
	yaux313 = -yaux312;
	yaux314 = yaux280 + yaux313;
	yaux315 = yaux314 * yaux314;
	yaux316 = x[29];
	yaux317 = -yaux316;
	yaux318 = yaux284 + yaux317;
	yaux319 = yaux318 * yaux318;
	yaux320 = yaux311 + yaux315 + yaux319;
	yaux321 = 1 / SQRT(yaux320);
	yaux322 = -0.55*yaux321;
	yaux323 = 12 + yaux322;
	yaux301 = yaux253 + yaux284;
	yaux308 = yaux277 + yaux307;
	yaux339 = x[30];
	yaux328 = yaux281 + yaux312;
	yaux344 = x[31];
	yaux341 = -yaux339;
	yaux342 = yaux307 + yaux341;
	yaux343 = yaux342 * yaux342;
	yaux345 = -yaux344;
	yaux346 = yaux312 + yaux345;
	yaux347 = yaux346 * yaux346;
	yaux348 = x[32];
	yaux349 = -yaux348;
	yaux350 = yaux316 + yaux349;
	yaux351 = yaux350 * yaux350;
	yaux352 = yaux343 + yaux347 + yaux351;
	yaux353 = 1 / SQRT(yaux352);
	yaux354 = -0.55*yaux353;
	yaux355 = 12 + yaux354;
	yaux333 = yaux285 + yaux316;
	yaux340 = yaux309 + yaux339;
	yaux371 = x[33];
	yaux360 = yaux313 + yaux344;
	yaux376 = x[34];
	yaux373 = -yaux371;
	yaux374 = yaux339 + yaux373;
	yaux375 = yaux374 * yaux374;
	yaux377 = -yaux376;
	yaux378 = yaux344 + yaux377;
	yaux379 = yaux378 * yaux378;
	yaux380 = x[35];
	yaux381 = -yaux380;
	yaux382 = yaux348 + yaux381;
	yaux383 = yaux382 * yaux382;
	yaux384 = yaux375 + yaux379 + yaux383;
	yaux385 = 1 / SQRT(yaux384);
	yaux386 = -0.55*yaux385;
	yaux387 = 12 + yaux386;
	yaux365 = yaux317 + yaux348;

	out[0] = x[36];
	out[1] = x[37];
	out[2] = x[38];
	out[3] = x[39];
	out[4] = x[40];
	out[5] = x[41];
	out[6] = x[42];
	out[7] = x[43];
	out[8] = x[44];
	out[9] = x[45];
	out[10] = x[46];
	out[11] = x[47];
	out[12] = x[48];
	out[13] = x[49];
	out[14] = x[50];
	out[15] = x[51];
	out[16] = x[52];
	out[17] = x[53];
	out[18] = x[54];
	out[19] = x[55];
	out[20] = x[56];
	out[21] = x[57];
	out[22] = x[58];
	out[23] = x[59];
	out[24] = x[60];
	out[25] = x[61];
	out[26] = x[62];
	out[27] = x[63];
	out[28] = x[64];
	out[29] = x[65];
	out[30] = x[66];
	out[31] = x[67];
	out[32] = x[68];
	out[33] = u[0];
	out[34] = u[1];
	out[35] = u[2];
	out[36] = 26.666666666666668*(-0.1*yaux37*yaux46 + 0.1*yaux50*yaux65);
	out[37] = 26.666666666666668*(-0.1*yaux39*yaux46 + 0.1*yaux65*yaux71);
	out[38] = -9.81 + 26.666666666666668*(-0.1*yaux41*yaux46 + 0.1*yaux65*yaux77);
	out[39] = 26.666666666666668*(-0.1*yaux50*yaux65 + 0.1*yaux84*yaux99);
	out[40] = 26.666666666666668*(-0.1*yaux65*yaux71 + 0.1*yaux104*yaux99);
	out[41] = -9.81 + 26.666666666666668*(-0.1*yaux65*yaux77 + 0.1*yaux109*yaux99);
	out[42] = 26.666666666666668*(0.1*yaux116*yaux131 - 0.1*yaux84*yaux99);
	out[43] = 26.666666666666668*(0.1*yaux131*yaux136 - 0.1*yaux104*yaux99);
	out[44] = -9.81 + 26.666666666666668*(0.1*yaux131*yaux141 - 0.1*yaux109*yaux99);
	out[45] = 26.666666666666668*(-0.1*yaux116*yaux131 + 0.1*yaux148*yaux163);
	out[46] = 26.666666666666668*(-0.1*yaux131*yaux136 + 0.1*yaux163*yaux168);
	out[47] = -9.81 + 26.666666666666668*(-0.1*yaux131*yaux141 + 0.1*yaux163*yaux173);
	out[48] = 26.666666666666668*(-0.1*yaux148*yaux163 + 0.1*yaux180*yaux195);
	out[49] = 26.666666666666668*(-0.1*yaux163*yaux168 + 0.1*yaux195*yaux200);
	out[50] = -9.81 + 26.666666666666668*(-0.1*yaux163*yaux173 + 0.1*yaux195*yaux205);
	out[51] = 26.666666666666668*(-0.1*yaux180*yaux195 + 0.1*yaux212*yaux227);
	out[52] = 26.666666666666668*(-0.1*yaux195*yaux200 + 0.1*yaux227*yaux232);
	out[53] = -9.81 + 26.666666666666668*(-0.1*yaux195*yaux205 + 0.1*yaux227*yaux237);
	out[54] = 26.666666666666668*(-0.1*yaux212*yaux227 + 0.1*yaux244*yaux259);
	out[55] = 26.666666666666668*(-0.1*yaux227*yaux232 + 0.1*yaux259*yaux264);
	out[56] = -9.81 + 26.666666666666668*(-0.1*yaux227*yaux237 + 0.1*yaux259*yaux269);
	out[57] = 26.666666666666668*(-0.1*yaux244*yaux259 + 0.1*yaux276*yaux291);
	out[58] = 26.666666666666668*(-0.1*yaux259*yaux264 + 0.1*yaux291*yaux296);
	out[59] = -9.81 + 26.666666666666668*(-0.1*yaux259*yaux269 + 0.1*yaux291*yaux301);
	out[60] = 26.666666666666668*(-0.1*yaux276*yaux291 + 0.1*yaux308*yaux323);
	out[61] = 26.666666666666668*(-0.1*yaux291*yaux296 + 0.1*yaux323*yaux328);
	out[62] = -9.81 + 26.666666666666668*(-0.1*yaux291*yaux301 + 0.1*yaux323*yaux333);
	out[63] = 26.666666666666668*(-0.1*yaux308*yaux323 + 0.1*yaux340*yaux355);
	out[64] = 26.666666666666668*(-0.1*yaux323*yaux328 + 0.1*yaux355*yaux360);
	out[65] = -9.81 + 26.666666666666668*(-0.1*yaux323*yaux333 + 0.1*yaux355*yaux365);
	out[66] = 26.666666666666668*(-0.1*yaux340*yaux355 + 0.1*(yaux341 + yaux371)*yaux387);
	out[67] = 26.666666666666668*(-0.1*yaux355*yaux360 + 0.1*(yaux345 + yaux376)*yaux387);
	out[68] = -9.81 + 26.666666666666668*(-0.1*yaux355*yaux365 + 0.1*(yaux349 + yaux380)*yaux387);
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
	typeRNum yaux665, yaux666, yaux667, yaux668, yaux669, yaux67, yaux670, yaux672, yaux673, yaux676;
	typeRNum yaux677, yaux680, yaux682, yaux683, yaux685, yaux692, yaux695, yaux696, yaux698, yaux7;
	typeRNum yaux703, yaux710, yaux713, yaux718, yaux719, yaux72, yaux721, yaux728, yaux729, yaux73;
	typeRNum yaux730, yaux731, yaux732, yaux733, yaux734, yaux735, yaux736, yaux737, yaux738, yaux739;
	typeRNum yaux740, yaux741, yaux742, yaux743, yaux744, yaux745, yaux746, yaux747, yaux748, yaux749;
	typeRNum yaux751, yaux752, yaux755, yaux756, yaux759, yaux761, yaux762, yaux764, yaux77, yaux771;
	typeRNum yaux774, yaux775, yaux777, yaux782, yaux789, yaux792, yaux797, yaux798, yaux8, yaux80;
	typeRNum yaux800, yaux807, yaux808, yaux809, yaux810, yaux811, yaux812, yaux813, yaux814, yaux815;
	typeRNum yaux816, yaux817, yaux818, yaux819, yaux820, yaux821, yaux822, yaux823, yaux824, yaux825;
	typeRNum yaux829, yaux833, yaux85, yaux86, yaux866, yaux867, yaux868, yaux869, yaux89, yaux9;
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
	yaux20 = 12 + yaux19;
	yaux22 = -yaux2;
	yaux23 = yaux22 + yaux3;
	yaux24 = POW(yaux17, -1.5);
	yaux33 = POW(yaux32, -1.5);
	yaux43 = -yaux7;
	yaux44 = yaux43 + yaux8;
	yaux51 = -yaux12;
	yaux52 = yaux13 + yaux51;
	yaux28 = vec[36];
	yaux42 = -0.05500000000000001*yaux2*yaux33*yaux7;
	yaux48 = vec[40];
	yaux21 = 0.1*yaux20;
	yaux41 = vec[37];
	yaux35 = 1 / SQRT(yaux32);
	yaux36 = 0.05500000000000001*yaux35;
	yaux37 = -0.1*yaux20;
	yaux1 = vec[39];
	yaux50 = vec[41];
	yaux54 = vec[38];
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
	yaux112 = 12 + yaux111;
	yaux114 = yaux4 + yaux97;
	yaux115 = POW(yaux109, -1.5);
	yaux45 = 0.05500000000000001*yaux24*yaux44*yaux5;
	yaux123 = yaux101 + yaux9;
	yaux56 = 0.05500000000000001*yaux24*yaux5*yaux52;
	yaux130 = yaux105 + yaux14;
	yaux63 = -0.05500000000000001*yaux10*yaux24*yaux44;
	yaux64 = yaux21 + yaux63;
	yaux60 = 0.05500000000000001*yaux10*yaux23*yaux24;
	yaux127 = vec[43];
	yaux113 = 0.1*yaux112;
	yaux67 = 0.05500000000000001*yaux10*yaux24*yaux44;
	yaux119 = -0.1*yaux112;
	yaux96 = vec[42];
	yaux129 = vec[44];
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
	yaux191 = 12 + yaux190;
	yaux193 = yaux176 + yaux98;
	yaux194 = POW(yaux188, -1.5);
	yaux124 = 0.05500000000000001*yaux115*yaux123*yaux99;
	yaux202 = yaux102 + yaux180;
	yaux132 = 0.05500000000000001*yaux115*yaux130*yaux99;
	yaux209 = yaux106 + yaux184;
	yaux142 = -0.05500000000000001*yaux103*yaux115*yaux123;
	yaux143 = yaux113 + yaux142;
	yaux139 = 0.05500000000000001*yaux103*yaux114*yaux115;
	yaux206 = vec[46];
	yaux192 = 0.1*yaux191;
	yaux145 = 0.05500000000000001*yaux103*yaux115*yaux123;
	yaux198 = -0.1*yaux191;
	yaux175 = vec[45];
	yaux208 = vec[47];
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
	yaux270 = 12 + yaux269;
	yaux272 = yaux177 + yaux255;
	yaux273 = POW(yaux267, -1.5);
	yaux203 = 0.05500000000000001*yaux178*yaux194*yaux202;
	yaux281 = yaux181 + yaux259;
	yaux211 = 0.05500000000000001*yaux178*yaux194*yaux209;
	yaux288 = yaux185 + yaux263;
	yaux221 = -0.05500000000000001*yaux182*yaux194*yaux202;
	yaux222 = yaux192 + yaux221;
	yaux218 = 0.05500000000000001*yaux182*yaux193*yaux194;
	yaux285 = vec[49];
	yaux271 = 0.1*yaux270;
	yaux224 = 0.05500000000000001*yaux182*yaux194*yaux202;
	yaux277 = -0.1*yaux270;
	yaux254 = vec[48];
	yaux287 = vec[50];
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
	yaux349 = 12 + yaux348;
	yaux351 = yaux256 + yaux334;
	yaux352 = POW(yaux346, -1.5);
	yaux282 = 0.05500000000000001*yaux257*yaux273*yaux281;
	yaux360 = yaux260 + yaux338;
	yaux290 = 0.05500000000000001*yaux257*yaux273*yaux288;
	yaux367 = yaux264 + yaux342;
	yaux300 = -0.05500000000000001*yaux261*yaux273*yaux281;
	yaux301 = yaux271 + yaux300;
	yaux297 = 0.05500000000000001*yaux261*yaux272*yaux273;
	yaux364 = vec[52];
	yaux350 = 0.1*yaux349;
	yaux303 = 0.05500000000000001*yaux261*yaux273*yaux281;
	yaux356 = -0.1*yaux349;
	yaux333 = vec[51];
	yaux366 = vec[53];
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
	yaux428 = 12 + yaux427;
	yaux430 = yaux335 + yaux413;
	yaux431 = POW(yaux425, -1.5);
	yaux361 = 0.05500000000000001*yaux336*yaux352*yaux360;
	yaux439 = yaux339 + yaux417;
	yaux369 = 0.05500000000000001*yaux336*yaux352*yaux367;
	yaux446 = yaux343 + yaux421;
	yaux379 = -0.05500000000000001*yaux340*yaux352*yaux360;
	yaux380 = yaux350 + yaux379;
	yaux376 = 0.05500000000000001*yaux340*yaux351*yaux352;
	yaux443 = vec[55];
	yaux429 = 0.1*yaux428;
	yaux382 = 0.05500000000000001*yaux340*yaux352*yaux360;
	yaux435 = -0.1*yaux428;
	yaux412 = vec[54];
	yaux445 = vec[56];
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
	yaux507 = 12 + yaux506;
	yaux509 = yaux414 + yaux492;
	yaux510 = POW(yaux504, -1.5);
	yaux440 = 0.05500000000000001*yaux415*yaux431*yaux439;
	yaux518 = yaux418 + yaux496;
	yaux448 = 0.05500000000000001*yaux415*yaux431*yaux446;
	yaux525 = yaux422 + yaux500;
	yaux458 = -0.05500000000000001*yaux419*yaux431*yaux439;
	yaux459 = yaux429 + yaux458;
	yaux455 = 0.05500000000000001*yaux419*yaux430*yaux431;
	yaux522 = vec[58];
	yaux508 = 0.1*yaux507;
	yaux461 = 0.05500000000000001*yaux419*yaux431*yaux439;
	yaux514 = -0.1*yaux507;
	yaux491 = vec[57];
	yaux524 = vec[59];
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
	yaux586 = 12 + yaux585;
	yaux588 = yaux493 + yaux571;
	yaux589 = POW(yaux583, -1.5);
	yaux519 = 0.05500000000000001*yaux494*yaux510*yaux518;
	yaux597 = yaux497 + yaux575;
	yaux527 = 0.05500000000000001*yaux494*yaux510*yaux525;
	yaux604 = yaux501 + yaux579;
	yaux537 = -0.05500000000000001*yaux498*yaux510*yaux518;
	yaux538 = yaux508 + yaux537;
	yaux534 = 0.05500000000000001*yaux498*yaux509*yaux510;
	yaux601 = vec[61];
	yaux587 = 0.1*yaux586;
	yaux540 = 0.05500000000000001*yaux498*yaux510*yaux518;
	yaux593 = -0.1*yaux586;
	yaux570 = vec[60];
	yaux603 = vec[62];
	yaux545 = 0.05500000000000001*yaux498*yaux510*yaux525;
	yaux560 = -0.05500000000000001*yaux502*yaux510*yaux525;
	yaux561 = yaux508 + yaux560;
	yaux552 = 0.05500000000000001*yaux502*yaux509*yaux510;
	yaux555 = 0.05500000000000001*yaux502*yaux510*yaux518;
	yaux563 = 0.05500000000000001*yaux502*yaux510*yaux525;
	yaux590 = -0.05500000000000001*yaux573*yaux588*yaux589;
	yaux591 = yaux587 + yaux590;
	yaux650 = x[27];
	yaux651 = -yaux650;
	yaux652 = yaux571 + yaux651;
	yaux653 = yaux652 * yaux652;
	yaux654 = x[28];
	yaux655 = -yaux654;
	yaux656 = yaux575 + yaux655;
	yaux657 = yaux656 * yaux656;
	yaux658 = x[29];
	yaux659 = -yaux658;
	yaux660 = yaux579 + yaux659;
	yaux661 = yaux660 * yaux660;
	yaux662 = yaux653 + yaux657 + yaux661;
	yaux594 = 0.05500000000000001*yaux573*yaux588*yaux589;
	yaux663 = 1 / SQRT(yaux662);
	yaux664 = -0.55*yaux663;
	yaux665 = 12 + yaux664;
	yaux667 = yaux572 + yaux650;
	yaux668 = POW(yaux662, -1.5);
	yaux598 = 0.05500000000000001*yaux573*yaux589*yaux597;
	yaux676 = yaux576 + yaux654;
	yaux606 = 0.05500000000000001*yaux573*yaux589*yaux604;
	yaux683 = yaux580 + yaux658;
	yaux616 = -0.05500000000000001*yaux577*yaux589*yaux597;
	yaux617 = yaux587 + yaux616;
	yaux613 = 0.05500000000000001*yaux577*yaux588*yaux589;
	yaux680 = vec[64];
	yaux666 = 0.1*yaux665;
	yaux619 = 0.05500000000000001*yaux577*yaux589*yaux597;
	yaux672 = -0.1*yaux665;
	yaux649 = vec[63];
	yaux682 = vec[65];
	yaux624 = 0.05500000000000001*yaux577*yaux589*yaux604;
	yaux639 = -0.05500000000000001*yaux581*yaux589*yaux604;
	yaux640 = yaux587 + yaux639;
	yaux631 = 0.05500000000000001*yaux581*yaux588*yaux589;
	yaux634 = 0.05500000000000001*yaux581*yaux589*yaux597;
	yaux642 = 0.05500000000000001*yaux581*yaux589*yaux604;
	yaux669 = -0.05500000000000001*yaux652*yaux667*yaux668;
	yaux670 = yaux666 + yaux669;
	yaux729 = x[30];
	yaux730 = -yaux729;
	yaux731 = yaux650 + yaux730;
	yaux732 = yaux731 * yaux731;
	yaux733 = x[31];
	yaux734 = -yaux733;
	yaux735 = yaux654 + yaux734;
	yaux736 = yaux735 * yaux735;
	yaux737 = x[32];
	yaux738 = -yaux737;
	yaux739 = yaux658 + yaux738;
	yaux740 = yaux739 * yaux739;
	yaux741 = yaux732 + yaux736 + yaux740;
	yaux673 = 0.05500000000000001*yaux652*yaux667*yaux668;
	yaux742 = 1 / SQRT(yaux741);
	yaux743 = -0.55*yaux742;
	yaux744 = 12 + yaux743;
	yaux746 = yaux651 + yaux729;
	yaux747 = POW(yaux741, -1.5);
	yaux677 = 0.05500000000000001*yaux652*yaux668*yaux676;
	yaux755 = yaux655 + yaux733;
	yaux685 = 0.05500000000000001*yaux652*yaux668*yaux683;
	yaux762 = yaux659 + yaux737;
	yaux695 = -0.05500000000000001*yaux656*yaux668*yaux676;
	yaux696 = yaux666 + yaux695;
	yaux692 = 0.05500000000000001*yaux656*yaux667*yaux668;
	yaux759 = vec[67];
	yaux745 = 0.1*yaux744;
	yaux698 = 0.05500000000000001*yaux656*yaux668*yaux676;
	yaux751 = -0.1*yaux744;
	yaux728 = vec[66];
	yaux761 = vec[68];
	yaux703 = 0.05500000000000001*yaux656*yaux668*yaux683;
	yaux718 = -0.05500000000000001*yaux660*yaux668*yaux683;
	yaux719 = yaux666 + yaux718;
	yaux710 = 0.05500000000000001*yaux660*yaux667*yaux668;
	yaux713 = 0.05500000000000001*yaux660*yaux668*yaux676;
	yaux721 = 0.05500000000000001*yaux660*yaux668*yaux683;
	yaux748 = -0.05500000000000001*yaux731*yaux746*yaux747;
	yaux749 = yaux745 + yaux748;
	yaux752 = 0.05500000000000001*yaux731*yaux746*yaux747;
	yaux807 = x[33];
	yaux808 = -yaux807;
	yaux809 = yaux729 + yaux808;
	yaux810 = yaux809 * yaux809;
	yaux811 = x[34];
	yaux812 = -yaux811;
	yaux813 = yaux733 + yaux812;
	yaux814 = yaux813 * yaux813;
	yaux815 = x[35];
	yaux816 = -yaux815;
	yaux817 = yaux737 + yaux816;
	yaux818 = yaux817 * yaux817;
	yaux819 = yaux810 + yaux814 + yaux818;
	yaux756 = 0.05500000000000001*yaux731*yaux747*yaux755;
	yaux825 = POW(yaux819, -1.5);
	yaux764 = 0.05500000000000001*yaux731*yaux747*yaux762;
	yaux774 = -0.05500000000000001*yaux735*yaux747*yaux755;
	yaux775 = yaux745 + yaux774;
	yaux771 = 0.05500000000000001*yaux735*yaux746*yaux747;
	yaux824 = yaux730 + yaux807;
	yaux777 = 0.05500000000000001*yaux735*yaux747*yaux755;
	yaux820 = 1 / SQRT(yaux819);
	yaux821 = -0.55*yaux820;
	yaux822 = 12 + yaux821;
	yaux823 = -0.1*yaux822;
	yaux829 = yaux734 + yaux811;
	yaux782 = 0.05500000000000001*yaux735*yaux747*yaux762;
	yaux833 = yaux738 + yaux815;
	yaux797 = -0.05500000000000001*yaux739*yaux747*yaux762;
	yaux798 = yaux745 + yaux797;
	yaux789 = 0.05500000000000001*yaux739*yaux746*yaux747;
	yaux792 = 0.05500000000000001*yaux739*yaux747*yaux755;
	yaux800 = 0.05500000000000001*yaux739*yaux747*yaux762;
	yaux866 = SQRT(yaux819);
	yaux867 = -yaux866;
	yaux868 = 0.04583333333333334 + yaux867;
	yaux869 = -21.81818181818182*yaux819*yaux868;

	out[0] = 26.666666666666668*yaux1*yaux26 + 26.666666666666668*yaux28*(-1.2000000000000002 - 0.05500000000000001*yaux29*yaux33 + yaux36 + yaux37 + yaux38) + 26.666666666666668*yaux41*(yaux42 + yaux45) - 1.466666666666667*yaux24*yaux44*yaux48*yaux5 - 1.466666666666667*yaux24*yaux5*yaux50*yaux52 + 26.666666666666668*yaux54*(yaux55 + yaux56);
	out[1] = -1.466666666666667*yaux1*yaux10*yaux23*yaux24 - 1.466666666666667*yaux10*yaux24*yaux50*yaux52 + 26.666666666666668*yaux28*(yaux42 + yaux60) + 26.666666666666668*yaux48*yaux64 + 26.666666666666668*yaux41*(-1.2000000000000002 - 0.05500000000000001*yaux30*yaux33 + yaux36 + yaux37 + yaux67) + 26.666666666666668*yaux54*(yaux72 + yaux73);
	out[2] = -1.466666666666667*yaux1*yaux15*yaux23*yaux24 - 1.466666666666667*yaux15*yaux24*yaux44*yaux48 + 26.666666666666668*yaux28*(yaux55 + yaux77) + 26.666666666666668*yaux41*(yaux72 + yaux80) + 26.666666666666668*yaux50*yaux86 + 26.666666666666668*yaux54*(-1.2000000000000002 - 0.05500000000000001*yaux31*yaux33 + yaux36 + yaux37 + yaux89);
	out[3] = 26.666666666666668*yaux26*yaux28 + 26.666666666666668*yaux1*(yaux119 + yaux120 + yaux37 + yaux38) + 26.666666666666668*(yaux124 + yaux45)*yaux48 - 1.466666666666667*yaux24*yaux41*yaux44*yaux5 - 1.466666666666667*yaux24*yaux5*yaux52*yaux54 + 26.666666666666668*yaux50*(yaux132 + yaux56) + 26.666666666666668*yaux117*yaux96 - 1.466666666666667*yaux115*yaux123*yaux127*yaux99 - 1.466666666666667*yaux115*yaux129*yaux130*yaux99;
	out[4] = -1.466666666666667*yaux103*yaux115*yaux129*yaux130 + 26.666666666666668*yaux127*yaux143 - 1.466666666666667*yaux10*yaux23*yaux24*yaux28 - 1.466666666666667*yaux10*yaux24*yaux52*yaux54 + 26.666666666666668*yaux1*(yaux139 + yaux60) + 26.666666666666668*yaux41*yaux64 + 26.666666666666668*yaux48*(yaux119 + yaux145 + yaux37 + yaux67) + 26.666666666666668*yaux50*(yaux150 + yaux73) - 1.466666666666667*yaux103*yaux114*yaux115*yaux96;
	out[5] = -1.466666666666667*yaux107*yaux115*yaux123*yaux127 + 26.666666666666668*yaux129*yaux166 - 1.466666666666667*yaux15*yaux23*yaux24*yaux28 - 1.466666666666667*yaux15*yaux24*yaux41*yaux44 + 26.666666666666668*yaux1*(yaux157 + yaux77) + 26.666666666666668*yaux48*(yaux160 + yaux80) + 26.666666666666668*yaux54*yaux86 + 26.666666666666668*yaux50*(yaux119 + yaux168 + yaux37 + yaux89) - 1.466666666666667*yaux107*yaux114*yaux115*yaux96;
	out[6] = 26.666666666666668*yaux1*yaux117 + 26.666666666666668*yaux175*yaux196 + 26.666666666666668*yaux127*(yaux124 + yaux203) - 1.466666666666667*yaux178*yaux194*yaux202*yaux206 - 1.466666666666667*yaux178*yaux194*yaux208*yaux209 + 26.666666666666668*yaux129*(yaux132 + yaux211) + 26.666666666666668*(yaux119 + yaux120 + yaux198 + yaux199)*yaux96 - 1.466666666666667*yaux115*yaux123*yaux48*yaux99 - 1.466666666666667*yaux115*yaux130*yaux50*yaux99;
	out[7] = -1.466666666666667*yaux1*yaux103*yaux114*yaux115 - 1.466666666666667*yaux175*yaux182*yaux193*yaux194 - 1.466666666666667*yaux182*yaux194*yaux208*yaux209 + 26.666666666666668*yaux206*yaux222 + 26.666666666666668*yaux127*(yaux119 + yaux145 + yaux198 + yaux224) + 26.666666666666668*yaux129*(yaux150 + yaux229) + 26.666666666666668*yaux143*yaux48 - 1.466666666666667*yaux103*yaux115*yaux130*yaux50 + 26.666666666666668*(yaux139 + yaux218)*yaux96;
	out[8] = -1.466666666666667*yaux1*yaux107*yaux114*yaux115 - 1.466666666666667*yaux175*yaux186*yaux193*yaux194 - 1.466666666666667*yaux186*yaux194*yaux202*yaux206 + 26.666666666666668*yaux127*(yaux160 + yaux239) + 26.666666666666668*yaux208*yaux245 + 26.666666666666668*yaux129*(yaux119 + yaux168 + yaux198 + yaux247) - 1.466666666666667*yaux107*yaux115*yaux123*yaux48 + 26.666666666666668*yaux166*yaux50 + 26.666666666666668*(yaux157 + yaux236)*yaux96;
	out[9] = -1.466666666666667*yaux127*yaux178*yaux194*yaux202 - 1.466666666666667*yaux129*yaux178*yaux194*yaux209 + 26.666666666666668*yaux254*yaux275 + 26.666666666666668*yaux175*(yaux198 + yaux199 + yaux277 + yaux278) + 26.666666666666668*yaux206*(yaux203 + yaux282) - 1.466666666666667*yaux257*yaux273*yaux281*yaux285 - 1.466666666666667*yaux257*yaux273*yaux287*yaux288 + 26.666666666666668*yaux208*(yaux211 + yaux290) + 26.666666666666668*yaux196*yaux96;
	out[10] = -1.466666666666667*yaux129*yaux182*yaux194*yaux209 + 26.666666666666668*yaux127*yaux222 - 1.466666666666667*yaux254*yaux261*yaux272*yaux273 - 1.466666666666667*yaux261*yaux273*yaux287*yaux288 + 26.666666666666668*yaux175*(yaux218 + yaux297) + 26.666666666666668*yaux285*yaux301 + 26.666666666666668*yaux206*(yaux198 + yaux224 + yaux277 + yaux303) + 26.666666666666668*yaux208*(yaux229 + yaux308) - 1.466666666666667*yaux182*yaux193*yaux194*yaux96;
	out[11] = -1.466666666666667*yaux127*yaux186*yaux194*yaux202 + 26.666666666666668*yaux129*yaux245 - 1.466666666666667*yaux254*yaux265*yaux272*yaux273 - 1.466666666666667*yaux265*yaux273*yaux281*yaux285 + 26.666666666666668*yaux175*(yaux236 + yaux315) + 26.666666666666668*yaux206*(yaux239 + yaux318) + 26.666666666666668*yaux287*yaux324 + 26.666666666666668*yaux208*(yaux198 + yaux247 + yaux277 + yaux326) - 1.466666666666667*yaux186*yaux193*yaux194*yaux96;
	out[12] = 26.666666666666668*yaux175*yaux275 - 1.466666666666667*yaux206*yaux257*yaux273*yaux281 - 1.466666666666667*yaux208*yaux257*yaux273*yaux288 + 26.666666666666668*yaux333*yaux354 + 26.666666666666668*yaux254*(yaux277 + yaux278 + yaux356 + yaux357) + 26.666666666666668*yaux285*(yaux282 + yaux361) - 1.466666666666667*yaux336*yaux352*yaux360*yaux364 - 1.466666666666667*yaux336*yaux352*yaux366*yaux367 + 26.666666666666668*yaux287*(yaux290 + yaux369);
	out[13] = -1.466666666666667*yaux175*yaux261*yaux272*yaux273 - 1.466666666666667*yaux208*yaux261*yaux273*yaux288 + 26.666666666666668*yaux206*yaux301 - 1.466666666666667*yaux333*yaux340*yaux351*yaux352 - 1.466666666666667*yaux340*yaux352*yaux366*yaux367 + 26.666666666666668*yaux254*(yaux297 + yaux376) + 26.666666666666668*yaux364*yaux380 + 26.666666666666668*yaux285*(yaux277 + yaux303 + yaux356 + yaux382) + 26.666666666666668*yaux287*(yaux308 + yaux387);
	out[14] = -1.466666666666667*yaux175*yaux265*yaux272*yaux273 - 1.466666666666667*yaux206*yaux265*yaux273*yaux281 + 26.666666666666668*yaux208*yaux324 - 1.466666666666667*yaux333*yaux344*yaux351*yaux352 - 1.466666666666667*yaux344*yaux352*yaux360*yaux364 + 26.666666666666668*yaux254*(yaux315 + yaux394) + 26.666666666666668*yaux285*(yaux318 + yaux397) + 26.666666666666668*yaux366*yaux403 + 26.666666666666668*yaux287*(yaux277 + yaux326 + yaux356 + yaux405);
	out[15] = 26.666666666666668*yaux254*yaux354 - 1.466666666666667*yaux285*yaux336*yaux352*yaux360 - 1.466666666666667*yaux287*yaux336*yaux352*yaux367 + 26.666666666666668*yaux412*yaux433 + 26.666666666666668*yaux333*(yaux356 + yaux357 + yaux435 + yaux436) + 26.666666666666668*yaux364*(yaux361 + yaux440) - 1.466666666666667*yaux415*yaux431*yaux439*yaux443 - 1.466666666666667*yaux415*yaux431*yaux445*yaux446 + 26.666666666666668*yaux366*(yaux369 + yaux448);
	out[16] = -1.466666666666667*yaux254*yaux340*yaux351*yaux352 - 1.466666666666667*yaux287*yaux340*yaux352*yaux367 + 26.666666666666668*yaux285*yaux380 - 1.466666666666667*yaux412*yaux419*yaux430*yaux431 - 1.466666666666667*yaux419*yaux431*yaux445*yaux446 + 26.666666666666668*yaux333*(yaux376 + yaux455) + 26.666666666666668*yaux443*yaux459 + 26.666666666666668*yaux364*(yaux356 + yaux382 + yaux435 + yaux461) + 26.666666666666668*yaux366*(yaux387 + yaux466);
	out[17] = -1.466666666666667*yaux254*yaux344*yaux351*yaux352 - 1.466666666666667*yaux285*yaux344*yaux352*yaux360 + 26.666666666666668*yaux287*yaux403 - 1.466666666666667*yaux412*yaux423*yaux430*yaux431 - 1.466666666666667*yaux423*yaux431*yaux439*yaux443 + 26.666666666666668*yaux333*(yaux394 + yaux473) + 26.666666666666668*yaux364*(yaux397 + yaux476) + 26.666666666666668*yaux445*yaux482 + 26.666666666666668*yaux366*(yaux356 + yaux405 + yaux435 + yaux484);
	out[18] = 26.666666666666668*yaux333*yaux433 - 1.466666666666667*yaux364*yaux415*yaux431*yaux439 - 1.466666666666667*yaux366*yaux415*yaux431*yaux446 + 26.666666666666668*yaux491*yaux512 + 26.666666666666668*yaux412*(yaux435 + yaux436 + yaux514 + yaux515) + 26.666666666666668*yaux443*(yaux440 + yaux519) - 1.466666666666667*yaux494*yaux510*yaux518*yaux522 - 1.466666666666667*yaux494*yaux510*yaux524*yaux525 + 26.666666666666668*yaux445*(yaux448 + yaux527);
	out[19] = -1.466666666666667*yaux333*yaux419*yaux430*yaux431 - 1.466666666666667*yaux366*yaux419*yaux431*yaux446 + 26.666666666666668*yaux364*yaux459 - 1.466666666666667*yaux491*yaux498*yaux509*yaux510 - 1.466666666666667*yaux498*yaux510*yaux524*yaux525 + 26.666666666666668*yaux412*(yaux455 + yaux534) + 26.666666666666668*yaux522*yaux538 + 26.666666666666668*yaux443*(yaux435 + yaux461 + yaux514 + yaux540) + 26.666666666666668*yaux445*(yaux466 + yaux545);
	out[20] = -1.466666666666667*yaux333*yaux423*yaux430*yaux431 - 1.466666666666667*yaux364*yaux423*yaux431*yaux439 + 26.666666666666668*yaux366*yaux482 - 1.466666666666667*yaux491*yaux502*yaux509*yaux510 - 1.466666666666667*yaux502*yaux510*yaux518*yaux522 + 26.666666666666668*yaux412*(yaux473 + yaux552) + 26.666666666666668*yaux443*(yaux476 + yaux555) + 26.666666666666668*yaux524*yaux561 + 26.666666666666668*yaux445*(yaux435 + yaux484 + yaux514 + yaux563);
	out[21] = 26.666666666666668*yaux412*yaux512 - 1.466666666666667*yaux443*yaux494*yaux510*yaux518 - 1.466666666666667*yaux445*yaux494*yaux510*yaux525 + 26.666666666666668*yaux570*yaux591 + 26.666666666666668*yaux491*(yaux514 + yaux515 + yaux593 + yaux594) + 26.666666666666668*yaux522*(yaux519 + yaux598) - 1.466666666666667*yaux573*yaux589*yaux597*yaux601 - 1.466666666666667*yaux573*yaux589*yaux603*yaux604 + 26.666666666666668*yaux524*(yaux527 + yaux606);
	out[22] = -1.466666666666667*yaux412*yaux498*yaux509*yaux510 - 1.466666666666667*yaux445*yaux498*yaux510*yaux525 + 26.666666666666668*yaux443*yaux538 - 1.466666666666667*yaux570*yaux577*yaux588*yaux589 - 1.466666666666667*yaux577*yaux589*yaux603*yaux604 + 26.666666666666668*yaux491*(yaux534 + yaux613) + 26.666666666666668*yaux601*yaux617 + 26.666666666666668*yaux522*(yaux514 + yaux540 + yaux593 + yaux619) + 26.666666666666668*yaux524*(yaux545 + yaux624);
	out[23] = -1.466666666666667*yaux412*yaux502*yaux509*yaux510 - 1.466666666666667*yaux443*yaux502*yaux510*yaux518 + 26.666666666666668*yaux445*yaux561 - 1.466666666666667*yaux570*yaux581*yaux588*yaux589 - 1.466666666666667*yaux581*yaux589*yaux597*yaux601 + 26.666666666666668*yaux491*(yaux552 + yaux631) + 26.666666666666668*yaux522*(yaux555 + yaux634) + 26.666666666666668*yaux603*yaux640 + 26.666666666666668*yaux524*(yaux514 + yaux563 + yaux593 + yaux642);
	out[24] = 26.666666666666668*yaux491*yaux591 - 1.466666666666667*yaux522*yaux573*yaux589*yaux597 - 1.466666666666667*yaux524*yaux573*yaux589*yaux604 + 26.666666666666668*yaux649*yaux670 + 26.666666666666668*yaux570*(yaux593 + yaux594 + yaux672 + yaux673) + 26.666666666666668*yaux601*(yaux598 + yaux677) - 1.466666666666667*yaux652*yaux668*yaux676*yaux680 - 1.466666666666667*yaux652*yaux668*yaux682*yaux683 + 26.666666666666668*yaux603*(yaux606 + yaux685);
	out[25] = -1.466666666666667*yaux491*yaux577*yaux588*yaux589 - 1.466666666666667*yaux524*yaux577*yaux589*yaux604 + 26.666666666666668*yaux522*yaux617 - 1.466666666666667*yaux649*yaux656*yaux667*yaux668 - 1.466666666666667*yaux656*yaux668*yaux682*yaux683 + 26.666666666666668*yaux570*(yaux613 + yaux692) + 26.666666666666668*yaux680*yaux696 + 26.666666666666668*yaux601*(yaux593 + yaux619 + yaux672 + yaux698) + 26.666666666666668*yaux603*(yaux624 + yaux703);
	out[26] = -1.466666666666667*yaux491*yaux581*yaux588*yaux589 - 1.466666666666667*yaux522*yaux581*yaux589*yaux597 + 26.666666666666668*yaux524*yaux640 - 1.466666666666667*yaux649*yaux660*yaux667*yaux668 - 1.466666666666667*yaux660*yaux668*yaux676*yaux680 + 26.666666666666668*yaux570*(yaux631 + yaux710) + 26.666666666666668*yaux601*(yaux634 + yaux713) + 26.666666666666668*yaux682*yaux719 + 26.666666666666668*yaux603*(yaux593 + yaux642 + yaux672 + yaux721);
	out[27] = 26.666666666666668*yaux570*yaux670 - 1.466666666666667*yaux601*yaux652*yaux668*yaux676 - 1.466666666666667*yaux603*yaux652*yaux668*yaux683 + 26.666666666666668*yaux728*yaux749 + 26.666666666666668*yaux649*(yaux672 + yaux673 + yaux751 + yaux752) + 26.666666666666668*yaux680*(yaux677 + yaux756) - 1.466666666666667*yaux731*yaux747*yaux755*yaux759 - 1.466666666666667*yaux731*yaux747*yaux761*yaux762 + 26.666666666666668*yaux682*(yaux685 + yaux764);
	out[28] = -1.466666666666667*yaux570*yaux656*yaux667*yaux668 - 1.466666666666667*yaux603*yaux656*yaux668*yaux683 + 26.666666666666668*yaux601*yaux696 - 1.466666666666667*yaux728*yaux735*yaux746*yaux747 - 1.466666666666667*yaux735*yaux747*yaux761*yaux762 + 26.666666666666668*yaux649*(yaux692 + yaux771) + 26.666666666666668*yaux759*yaux775 + 26.666666666666668*yaux680*(yaux672 + yaux698 + yaux751 + yaux777) + 26.666666666666668*yaux682*(yaux703 + yaux782);
	out[29] = -1.466666666666667*yaux570*yaux660*yaux667*yaux668 - 1.466666666666667*yaux601*yaux660*yaux668*yaux676 + 26.666666666666668*yaux603*yaux719 - 1.466666666666667*yaux728*yaux739*yaux746*yaux747 - 1.466666666666667*yaux739*yaux747*yaux755*yaux759 + 26.666666666666668*yaux649*(yaux710 + yaux789) + 26.666666666666668*yaux680*(yaux713 + yaux792) + 26.666666666666668*yaux761*yaux798 + 26.666666666666668*yaux682*(yaux672 + yaux721 + yaux751 + yaux800);
	out[30] = 26.666666666666668*yaux649*yaux749 - 1.466666666666667*yaux680*yaux731*yaux747*yaux755 - 1.466666666666667*yaux682*yaux731*yaux747*yaux762 + 26.666666666666668*yaux728*(yaux751 + yaux752 + yaux823 + 0.05500000000000001*yaux809*yaux824*yaux825) + 26.666666666666668*yaux759*(yaux756 + 0.05500000000000001*yaux809*yaux825*yaux829) + 26.666666666666668*yaux761*(yaux764 + 0.05500000000000001*yaux809*yaux825*yaux833);
	out[31] = -1.466666666666667*yaux649*yaux735*yaux746*yaux747 - 1.466666666666667*yaux682*yaux735*yaux747*yaux762 + 26.666666666666668*yaux680*yaux775 + 26.666666666666668*yaux728*(yaux771 + 0.05500000000000001*yaux813*yaux824*yaux825) + 26.666666666666668*yaux759*(yaux751 + yaux777 + yaux823 + 0.05500000000000001*yaux813*yaux825*yaux829) + 26.666666666666668*yaux761*(yaux782 + 0.05500000000000001*yaux813*yaux825*yaux833);
	out[32] = -1.466666666666667*yaux649*yaux739*yaux746*yaux747 - 1.466666666666667*yaux680*yaux739*yaux747*yaux755 + 26.666666666666668*yaux682*yaux798 + 26.666666666666668*yaux728*(yaux789 + 0.05500000000000001*yaux817*yaux824*yaux825) + 26.666666666666668*yaux759*(yaux792 + 0.05500000000000001*yaux817*yaux825*yaux829) + 26.666666666666668*yaux761*(yaux751 + yaux800 + yaux823 + 0.05500000000000001*yaux817*yaux825*yaux833);
	out[33] = 1.466666666666667*yaux825*(1 * yaux759*yaux809*yaux813 + 1 * yaux761*yaux809*yaux817 + 1 * yaux728*(1 * yaux810 + yaux869));
	out[34] = 1.466666666666667*yaux825*(1 * yaux728*yaux809*yaux813 + 1 * yaux761*yaux813*yaux817 + 1 * yaux759*(1 * yaux814 + yaux869));
	out[35] = -32.00000000000001*yaux825*(-0.04583333333333334*yaux728*yaux809*yaux817 - 0.04583333333333334*yaux759*yaux813*yaux817 + 1 * yaux761*(-0.04583333333333333*yaux818 + 1 * yaux819*yaux868));
	out[36] = vec[0];
	out[37] = vec[1];
	out[38] = vec[2];
	out[39] = vec[3];
	out[40] = vec[4];
	out[41] = vec[5];
	out[42] = vec[6];
	out[43] = vec[7];
	out[44] = vec[8];
	out[45] = vec[9];
	out[46] = vec[10];
	out[47] = vec[11];
	out[48] = vec[12];
	out[49] = vec[13];
	out[50] = vec[14];
	out[51] = vec[15];
	out[52] = vec[16];
	out[53] = vec[17];
	out[54] = vec[18];
	out[55] = vec[19];
	out[56] = vec[20];
	out[57] = vec[21];
	out[58] = vec[22];
	out[59] = vec[23];
	out[60] = vec[24];
	out[61] = vec[25];
	out[62] = vec[26];
	out[63] = vec[27];
	out[64] = vec[28];
	out[65] = vec[29];
	out[66] = vec[30];
	out[67] = vec[31];
	out[68] = vec[32];
}
/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	out[0] = vec[33];
	out[1] = vec[34];
	out[2] = vec[35];
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

	yaux11 = x[33];

	out[0] = param[2] * (POW2(u[0]) + POW2(u[1]) + POW2(u[2]))
		+ param[1] * (56.25 - 15 * yaux11 + yaux11 * yaux11 + POW2(x[34]) + POW2(x[35]))
		+ param[0] * (POW2(x[36]) + POW2(x[37]) + POW2(x[38]) + POW2(x[39]) + POW2(x[40]) + POW2(x[41]) + POW2(x[42])
			+ POW2(x[43]) + POW2(x[44]) + POW2(x[45]) + POW2(x[46]) + POW2(x[47]) + POW2(x[48]) + POW2(x[49]) + POW2(x[50])
			+ POW2(x[51]) + POW2(x[52]) + POW2(x[53]) + POW2(x[54]) + POW2(x[55]) + POW2(x[56]) + POW2(x[57]) + POW2(x[58])
			+ POW2(x[59]) + POW2(x[60]) + POW2(x[61]) + POW2(x[62]) + POW2(x[63]) + POW2(x[64]) + POW2(x[65]) + POW2(x[66])
			+ POW2(x[67]) + POW2(x[68]));
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
	out[27] = 0;
	out[28] = 0;
	out[29] = 0;
	out[30] = 0;
	out[31] = 0;
	out[32] = 0;
	out[33] = 2 * yaux1*(-7.5 + x[33]);
	out[34] = 2 * yaux1*x[34];
	out[35] = 2 * yaux1*x[35];
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
	out[57] = 2 * yaux9*x[57];
	out[58] = 2 * yaux9*x[58];
	out[59] = 2 * yaux9*x[59];
	out[60] = 2 * yaux9*x[60];
	out[61] = 2 * yaux9*x[61];
	out[62] = 2 * yaux9*x[62];
	out[63] = 2 * yaux9*x[63];
	out[64] = 2 * yaux9*x[64];
	out[65] = 2 * yaux9*x[65];
	out[66] = 2 * yaux9*x[66];
	out[67] = 2 * yaux9*x[67];
	out[68] = 2 * yaux9*x[68];
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

	out[0] = param[3] * (POW2(-7.5 + x[33]) + POW2(x[34]) + POW2(x[35]));
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
	out[27] = 0;
	out[28] = 0;
	out[29] = 0;
	out[30] = 0;
	out[31] = 0;
	out[32] = 0;
	out[33] = 2 * yaux1*(-7.5 + x[33]);
	out[34] = 2 * yaux1*x[34];
	out[35] = 2 * yaux1*x[35];
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
	out[57] = 0;
	out[58] = 0;
	out[59] = 0;
	out[60] = 0;
	out[61] = 0;
	out[62] = 0;
	out[63] = 0;
	out[64] = 0;
	out[65] = 0;
	out[66] = 0;
	out[67] = 0;
	out[68] = 0;
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
