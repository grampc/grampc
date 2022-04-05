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
	*Nx = 45;
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
	typeRNum yaux103, yaux104, yaux105, yaux106, yaux107, yaux108, yaux109, yaux110, yaux111, yaux112;
	typeRNum yaux113, yaux114, yaux115, yaux116, yaux117, yaux118, yaux119, yaux124, yaux129, yaux135;
	typeRNum yaux136, yaux137, yaux138, yaux139, yaux140, yaux141, yaux142, yaux143, yaux144, yaux145;
	typeRNum yaux146, yaux147, yaux148, yaux149, yaux150, yaux151, yaux156, yaux161, yaux167, yaux168;
	typeRNum yaux169, yaux170, yaux171, yaux172, yaux173, yaux174, yaux175, yaux176, yaux177, yaux178;
	typeRNum yaux179, yaux180, yaux181, yaux182, yaux183, yaux188, yaux193, yaux199, yaux200, yaux201;
	typeRNum yaux202, yaux203, yaux204, yaux205, yaux206, yaux207, yaux208, yaux209, yaux210, yaux211;
	typeRNum yaux212, yaux213, yaux214, yaux215, yaux220, yaux225, yaux231, yaux233, yaux234, yaux235;
	typeRNum yaux236, yaux237, yaux238, yaux239, yaux240, yaux241, yaux242, yaux243, yaux244, yaux245;
	typeRNum yaux246, yaux247, yaux25, yaux26, yaux27, yaux28, yaux29, yaux30, yaux31, yaux32;
	typeRNum yaux33, yaux34, yaux36, yaux37, yaux38, yaux39, yaux40, yaux41, yaux42, yaux43;
	typeRNum yaux44, yaux45, yaux46, yaux47, yaux48, yaux49, yaux50, yaux51, yaux52, yaux53;
	typeRNum yaux58, yaux59, yaux64, yaux65, yaux71, yaux72, yaux73, yaux74, yaux75, yaux76;
	typeRNum yaux77, yaux78, yaux79, yaux80, yaux81, yaux82, yaux83, yaux84, yaux85, yaux86;
	typeRNum yaux87, yaux92, yaux97;

	yaux25 = x[0];
	yaux37 = x[3];
	yaux27 = x[1];
	yaux29 = x[2];
	yaux26 = yaux25 * yaux25;
	yaux28 = yaux27 * yaux27;
	yaux30 = yaux29 * yaux29;
	yaux31 = yaux26 + yaux28 + yaux30;
	yaux32 = 1 / SQRT(yaux31);
	yaux33 = -0.55*yaux32;
	yaux34 = 8 + yaux33;
	yaux42 = x[4];
	yaux39 = -yaux37;
	yaux40 = yaux25 + yaux39;
	yaux41 = yaux40 * yaux40;
	yaux43 = -yaux42;
	yaux44 = yaux27 + yaux43;
	yaux45 = yaux44 * yaux44;
	yaux46 = x[5];
	yaux47 = -yaux46;
	yaux48 = yaux29 + yaux47;
	yaux49 = yaux48 * yaux48;
	yaux50 = yaux41 + yaux45 + yaux49;
	yaux51 = 1 / SQRT(yaux50);
	yaux52 = -0.55*yaux51;
	yaux53 = 8 + yaux52;
	yaux36 = -yaux25;
	yaux38 = yaux36 + yaux37;
	yaux71 = x[6];
	yaux58 = -yaux27;
	yaux59 = yaux42 + yaux58;
	yaux76 = x[7];
	yaux73 = -yaux71;
	yaux74 = yaux37 + yaux73;
	yaux75 = yaux74 * yaux74;
	yaux77 = -yaux76;
	yaux78 = yaux42 + yaux77;
	yaux79 = yaux78 * yaux78;
	yaux80 = x[8];
	yaux81 = -yaux80;
	yaux82 = yaux46 + yaux81;
	yaux83 = yaux82 * yaux82;
	yaux84 = yaux75 + yaux79 + yaux83;
	yaux85 = 1 / SQRT(yaux84);
	yaux86 = -0.55*yaux85;
	yaux87 = 8 + yaux86;
	yaux64 = -yaux29;
	yaux65 = yaux46 + yaux64;
	yaux72 = yaux39 + yaux71;
	yaux103 = x[9];
	yaux92 = yaux43 + yaux76;
	yaux108 = x[10];
	yaux105 = -yaux103;
	yaux106 = yaux105 + yaux71;
	yaux107 = yaux106 * yaux106;
	yaux109 = -yaux108;
	yaux110 = yaux109 + yaux76;
	yaux111 = yaux110 * yaux110;
	yaux112 = x[11];
	yaux113 = -yaux112;
	yaux114 = yaux113 + yaux80;
	yaux115 = yaux114 * yaux114;
	yaux116 = yaux107 + yaux111 + yaux115;
	yaux117 = 1 / SQRT(yaux116);
	yaux118 = -0.55*yaux117;
	yaux119 = 8 + yaux118;
	yaux97 = yaux47 + yaux80;
	yaux104 = yaux103 + yaux73;
	yaux135 = x[12];
	yaux124 = yaux108 + yaux77;
	yaux140 = x[13];
	yaux137 = -yaux135;
	yaux138 = yaux103 + yaux137;
	yaux139 = yaux138 * yaux138;
	yaux141 = -yaux140;
	yaux142 = yaux108 + yaux141;
	yaux143 = yaux142 * yaux142;
	yaux144 = x[14];
	yaux145 = -yaux144;
	yaux146 = yaux112 + yaux145;
	yaux147 = yaux146 * yaux146;
	yaux148 = yaux139 + yaux143 + yaux147;
	yaux149 = 1 / SQRT(yaux148);
	yaux150 = -0.55*yaux149;
	yaux151 = 8 + yaux150;
	yaux129 = yaux112 + yaux81;
	yaux136 = yaux105 + yaux135;
	yaux167 = x[15];
	yaux156 = yaux109 + yaux140;
	yaux172 = x[16];
	yaux169 = -yaux167;
	yaux170 = yaux135 + yaux169;
	yaux171 = yaux170 * yaux170;
	yaux173 = -yaux172;
	yaux174 = yaux140 + yaux173;
	yaux175 = yaux174 * yaux174;
	yaux176 = x[17];
	yaux177 = -yaux176;
	yaux178 = yaux144 + yaux177;
	yaux179 = yaux178 * yaux178;
	yaux180 = yaux171 + yaux175 + yaux179;
	yaux181 = 1 / SQRT(yaux180);
	yaux182 = -0.55*yaux181;
	yaux183 = 8 + yaux182;
	yaux161 = yaux113 + yaux144;
	yaux168 = yaux137 + yaux167;
	yaux199 = x[18];
	yaux188 = yaux141 + yaux172;
	yaux204 = x[19];
	yaux201 = -yaux199;
	yaux202 = yaux167 + yaux201;
	yaux203 = yaux202 * yaux202;
	yaux205 = -yaux204;
	yaux206 = yaux172 + yaux205;
	yaux207 = yaux206 * yaux206;
	yaux208 = x[20];
	yaux209 = -yaux208;
	yaux210 = yaux176 + yaux209;
	yaux211 = yaux210 * yaux210;
	yaux212 = yaux203 + yaux207 + yaux211;
	yaux213 = 1 / SQRT(yaux212);
	yaux214 = -0.55*yaux213;
	yaux215 = 8 + yaux214;
	yaux193 = yaux145 + yaux176;
	yaux200 = yaux169 + yaux199;
	yaux231 = x[21];
	yaux220 = yaux173 + yaux204;
	yaux236 = x[22];
	yaux233 = -yaux231;
	yaux234 = yaux199 + yaux233;
	yaux235 = yaux234 * yaux234;
	yaux237 = -yaux236;
	yaux238 = yaux204 + yaux237;
	yaux239 = yaux238 * yaux238;
	yaux240 = x[23];
	yaux241 = -yaux240;
	yaux242 = yaux208 + yaux241;
	yaux243 = yaux242 * yaux242;
	yaux244 = yaux235 + yaux239 + yaux243;
	yaux245 = 1 / SQRT(yaux244);
	yaux246 = -0.55*yaux245;
	yaux247 = 8 + yaux246;
	yaux225 = yaux177 + yaux208;

	out[0] = x[24];
	out[1] = x[25];
	out[2] = x[26];
	out[3] = x[27];
	out[4] = x[28];
	out[5] = x[29];
	out[6] = x[30];
	out[7] = x[31];
	out[8] = x[32];
	out[9] = x[33];
	out[10] = x[34];
	out[11] = x[35];
	out[12] = x[36];
	out[13] = x[37];
	out[14] = x[38];
	out[15] = x[39];
	out[16] = x[40];
	out[17] = x[41];
	out[18] = x[42];
	out[19] = x[43];
	out[20] = x[44];
	out[21] = u[0];
	out[22] = u[1];
	out[23] = u[2];
	out[24] = 17.77777777777778*(-0.1*yaux25*yaux34 + 0.1*yaux38*yaux53);
	out[25] = 17.77777777777778*(-0.1*yaux27*yaux34 + 0.1*yaux53*yaux59);
	out[26] = -9.81 + 17.77777777777778*(-0.1*yaux29*yaux34 + 0.1*yaux53*yaux65);
	out[27] = 17.77777777777778*(-0.1*yaux38*yaux53 + 0.1*yaux72*yaux87);
	out[28] = 17.77777777777778*(-0.1*yaux53*yaux59 + 0.1*yaux87*yaux92);
	out[29] = -9.81 + 17.77777777777778*(-0.1*yaux53*yaux65 + 0.1*yaux87*yaux97);
	out[30] = 17.77777777777778*(0.1*yaux104*yaux119 - 0.1*yaux72*yaux87);
	out[31] = 17.77777777777778*(0.1*yaux119*yaux124 - 0.1*yaux87*yaux92);
	out[32] = -9.81 + 17.77777777777778*(0.1*yaux119*yaux129 - 0.1*yaux87*yaux97);
	out[33] = 17.77777777777778*(-0.1*yaux104*yaux119 + 0.1*yaux136*yaux151);
	out[34] = 17.77777777777778*(-0.1*yaux119*yaux124 + 0.1*yaux151*yaux156);
	out[35] = -9.81 + 17.77777777777778*(-0.1*yaux119*yaux129 + 0.1*yaux151*yaux161);
	out[36] = 17.77777777777778*(-0.1*yaux136*yaux151 + 0.1*yaux168*yaux183);
	out[37] = 17.77777777777778*(-0.1*yaux151*yaux156 + 0.1*yaux183*yaux188);
	out[38] = -9.81 + 17.77777777777778*(-0.1*yaux151*yaux161 + 0.1*yaux183*yaux193);
	out[39] = 17.77777777777778*(-0.1*yaux168*yaux183 + 0.1*yaux200*yaux215);
	out[40] = 17.77777777777778*(-0.1*yaux183*yaux188 + 0.1*yaux215*yaux220);
	out[41] = -9.81 + 17.77777777777778*(-0.1*yaux183*yaux193 + 0.1*yaux215*yaux225);
	out[42] = 17.77777777777778*(-0.1*yaux200*yaux215 + 0.1*(yaux201 + yaux231)*yaux247);
	out[43] = 17.77777777777778*(-0.1*yaux215*yaux220 + 0.1*(yaux205 + yaux236)*yaux247);
	out[44] = -9.81 + 17.77777777777778*(-0.1*yaux215*yaux225 + 0.1*(yaux209 + yaux240)*yaux247);
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
	typeRNum yaux506, yaux507, yaux508, yaux509, yaux51, yaux513, yaux517, yaux52, yaux54, yaux55;
	typeRNum yaux550, yaux551, yaux552, yaux553, yaux56, yaux6, yaux60, yaux63, yaux64, yaux67;
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
	yaux20 = 8 + yaux19;
	yaux22 = -yaux2;
	yaux23 = yaux22 + yaux3;
	yaux24 = POW(yaux17, -1.5);
	yaux33 = POW(yaux32, -1.5);
	yaux43 = -yaux7;
	yaux44 = yaux43 + yaux8;
	yaux51 = -yaux12;
	yaux52 = yaux13 + yaux51;
	yaux28 = vec[24];
	yaux42 = -0.05500000000000001*yaux2*yaux33*yaux7;
	yaux48 = vec[28];
	yaux21 = 0.1*yaux20;
	yaux41 = vec[25];
	yaux35 = 1 / SQRT(yaux32);
	yaux36 = 0.05500000000000001*yaux35;
	yaux37 = -0.1*yaux20;
	yaux1 = vec[27];
	yaux50 = vec[29];
	yaux54 = vec[26];
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
	yaux112 = 8 + yaux111;
	yaux114 = yaux4 + yaux97;
	yaux115 = POW(yaux109, -1.5);
	yaux45 = 0.05500000000000001*yaux24*yaux44*yaux5;
	yaux123 = yaux101 + yaux9;
	yaux56 = 0.05500000000000001*yaux24*yaux5*yaux52;
	yaux130 = yaux105 + yaux14;
	yaux63 = -0.05500000000000001*yaux10*yaux24*yaux44;
	yaux64 = yaux21 + yaux63;
	yaux60 = 0.05500000000000001*yaux10*yaux23*yaux24;
	yaux127 = vec[31];
	yaux113 = 0.1*yaux112;
	yaux67 = 0.05500000000000001*yaux10*yaux24*yaux44;
	yaux119 = -0.1*yaux112;
	yaux96 = vec[30];
	yaux129 = vec[32];
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
	yaux191 = 8 + yaux190;
	yaux193 = yaux176 + yaux98;
	yaux194 = POW(yaux188, -1.5);
	yaux124 = 0.05500000000000001*yaux115*yaux123*yaux99;
	yaux202 = yaux102 + yaux180;
	yaux132 = 0.05500000000000001*yaux115*yaux130*yaux99;
	yaux209 = yaux106 + yaux184;
	yaux142 = -0.05500000000000001*yaux103*yaux115*yaux123;
	yaux143 = yaux113 + yaux142;
	yaux139 = 0.05500000000000001*yaux103*yaux114*yaux115;
	yaux206 = vec[34];
	yaux192 = 0.1*yaux191;
	yaux145 = 0.05500000000000001*yaux103*yaux115*yaux123;
	yaux198 = -0.1*yaux191;
	yaux175 = vec[33];
	yaux208 = vec[35];
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
	yaux270 = 8 + yaux269;
	yaux272 = yaux177 + yaux255;
	yaux273 = POW(yaux267, -1.5);
	yaux203 = 0.05500000000000001*yaux178*yaux194*yaux202;
	yaux281 = yaux181 + yaux259;
	yaux211 = 0.05500000000000001*yaux178*yaux194*yaux209;
	yaux288 = yaux185 + yaux263;
	yaux221 = -0.05500000000000001*yaux182*yaux194*yaux202;
	yaux222 = yaux192 + yaux221;
	yaux218 = 0.05500000000000001*yaux182*yaux193*yaux194;
	yaux285 = vec[37];
	yaux271 = 0.1*yaux270;
	yaux224 = 0.05500000000000001*yaux182*yaux194*yaux202;
	yaux277 = -0.1*yaux270;
	yaux254 = vec[36];
	yaux287 = vec[38];
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
	yaux349 = 8 + yaux348;
	yaux351 = yaux256 + yaux334;
	yaux352 = POW(yaux346, -1.5);
	yaux282 = 0.05500000000000001*yaux257*yaux273*yaux281;
	yaux360 = yaux260 + yaux338;
	yaux290 = 0.05500000000000001*yaux257*yaux273*yaux288;
	yaux367 = yaux264 + yaux342;
	yaux300 = -0.05500000000000001*yaux261*yaux273*yaux281;
	yaux301 = yaux271 + yaux300;
	yaux297 = 0.05500000000000001*yaux261*yaux272*yaux273;
	yaux364 = vec[40];
	yaux350 = 0.1*yaux349;
	yaux303 = 0.05500000000000001*yaux261*yaux273*yaux281;
	yaux356 = -0.1*yaux349;
	yaux333 = vec[39];
	yaux366 = vec[41];
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
	yaux428 = 8 + yaux427;
	yaux430 = yaux335 + yaux413;
	yaux431 = POW(yaux425, -1.5);
	yaux361 = 0.05500000000000001*yaux336*yaux352*yaux360;
	yaux439 = yaux339 + yaux417;
	yaux369 = 0.05500000000000001*yaux336*yaux352*yaux367;
	yaux446 = yaux343 + yaux421;
	yaux379 = -0.05500000000000001*yaux340*yaux352*yaux360;
	yaux380 = yaux350 + yaux379;
	yaux376 = 0.05500000000000001*yaux340*yaux351*yaux352;
	yaux443 = vec[43];
	yaux429 = 0.1*yaux428;
	yaux382 = 0.05500000000000001*yaux340*yaux352*yaux360;
	yaux435 = -0.1*yaux428;
	yaux412 = vec[42];
	yaux445 = vec[44];
	yaux387 = 0.05500000000000001*yaux340*yaux352*yaux367;
	yaux402 = -0.05500000000000001*yaux344*yaux352*yaux367;
	yaux403 = yaux350 + yaux402;
	yaux394 = 0.05500000000000001*yaux344*yaux351*yaux352;
	yaux397 = 0.05500000000000001*yaux344*yaux352*yaux360;
	yaux405 = 0.05500000000000001*yaux344*yaux352*yaux367;
	yaux432 = -0.05500000000000001*yaux415*yaux430*yaux431;
	yaux433 = yaux429 + yaux432;
	yaux436 = 0.05500000000000001*yaux415*yaux430*yaux431;
	yaux491 = x[21];
	yaux492 = -yaux491;
	yaux493 = yaux413 + yaux492;
	yaux494 = yaux493 * yaux493;
	yaux495 = x[22];
	yaux496 = -yaux495;
	yaux497 = yaux417 + yaux496;
	yaux498 = yaux497 * yaux497;
	yaux499 = x[23];
	yaux500 = -yaux499;
	yaux501 = yaux421 + yaux500;
	yaux502 = yaux501 * yaux501;
	yaux503 = yaux494 + yaux498 + yaux502;
	yaux440 = 0.05500000000000001*yaux415*yaux431*yaux439;
	yaux509 = POW(yaux503, -1.5);
	yaux448 = 0.05500000000000001*yaux415*yaux431*yaux446;
	yaux458 = -0.05500000000000001*yaux419*yaux431*yaux439;
	yaux459 = yaux429 + yaux458;
	yaux455 = 0.05500000000000001*yaux419*yaux430*yaux431;
	yaux508 = yaux414 + yaux491;
	yaux461 = 0.05500000000000001*yaux419*yaux431*yaux439;
	yaux504 = 1 / SQRT(yaux503);
	yaux505 = -0.55*yaux504;
	yaux506 = 8 + yaux505;
	yaux507 = -0.1*yaux506;
	yaux513 = yaux418 + yaux495;
	yaux466 = 0.05500000000000001*yaux419*yaux431*yaux446;
	yaux517 = yaux422 + yaux499;
	yaux481 = -0.05500000000000001*yaux423*yaux431*yaux446;
	yaux482 = yaux429 + yaux481;
	yaux473 = 0.05500000000000001*yaux423*yaux430*yaux431;
	yaux476 = 0.05500000000000001*yaux423*yaux431*yaux439;
	yaux484 = 0.05500000000000001*yaux423*yaux431*yaux446;
	yaux550 = SQRT(yaux503);
	yaux551 = -yaux550;
	yaux552 = 0.06875 + yaux551;
	yaux553 = -14.545454545454545*yaux503*yaux552;

	out[0] = 17.77777777777778*yaux1*yaux26 + 17.77777777777778*yaux28*(-0.8 - 0.05500000000000001*yaux29*yaux33 + yaux36 + yaux37 + yaux38) + 17.77777777777778*yaux41*(yaux42 + yaux45) - 0.977777777777778*yaux24*yaux44*yaux48*yaux5 - 0.977777777777778*yaux24*yaux5*yaux50*yaux52 + 17.77777777777778*yaux54*(yaux55 + yaux56);
	out[1] = -0.977777777777778*yaux1*yaux10*yaux23*yaux24 - 0.977777777777778*yaux10*yaux24*yaux50*yaux52 + 17.77777777777778*yaux28*(yaux42 + yaux60) + 17.77777777777778*yaux48*yaux64 + 17.77777777777778*yaux41*(-0.8 - 0.05500000000000001*yaux30*yaux33 + yaux36 + yaux37 + yaux67) + 17.77777777777778*yaux54*(yaux72 + yaux73);
	out[2] = -0.977777777777778*yaux1*yaux15*yaux23*yaux24 - 0.977777777777778*yaux15*yaux24*yaux44*yaux48 + 17.77777777777778*yaux28*(yaux55 + yaux77) + 17.77777777777778*yaux41*(yaux72 + yaux80) + 17.77777777777778*yaux50*yaux86 + 17.77777777777778*yaux54*(-0.8 - 0.05500000000000001*yaux31*yaux33 + yaux36 + yaux37 + yaux89);
	out[3] = 17.77777777777778*yaux26*yaux28 + 17.77777777777778*yaux1*(yaux119 + yaux120 + yaux37 + yaux38) + 17.77777777777778*(yaux124 + yaux45)*yaux48 - 0.977777777777778*yaux24*yaux41*yaux44*yaux5 - 0.977777777777778*yaux24*yaux5*yaux52*yaux54 + 17.77777777777778*yaux50*(yaux132 + yaux56) + 17.77777777777778*yaux117*yaux96 - 0.977777777777778*yaux115*yaux123*yaux127*yaux99 - 0.977777777777778*yaux115*yaux129*yaux130*yaux99;
	out[4] = -0.977777777777778*yaux103*yaux115*yaux129*yaux130 + 17.77777777777778*yaux127*yaux143 - 0.977777777777778*yaux10*yaux23*yaux24*yaux28 - 0.977777777777778*yaux10*yaux24*yaux52*yaux54 + 17.77777777777778*yaux1*(yaux139 + yaux60) + 17.77777777777778*yaux41*yaux64 + 17.77777777777778*yaux48*(yaux119 + yaux145 + yaux37 + yaux67) + 17.77777777777778*yaux50*(yaux150 + yaux73) - 0.977777777777778*yaux103*yaux114*yaux115*yaux96;
	out[5] = -0.977777777777778*yaux107*yaux115*yaux123*yaux127 + 17.77777777777778*yaux129*yaux166 - 0.977777777777778*yaux15*yaux23*yaux24*yaux28 - 0.977777777777778*yaux15*yaux24*yaux41*yaux44 + 17.77777777777778*yaux1*(yaux157 + yaux77) + 17.77777777777778*yaux48*(yaux160 + yaux80) + 17.77777777777778*yaux54*yaux86 + 17.77777777777778*yaux50*(yaux119 + yaux168 + yaux37 + yaux89) - 0.977777777777778*yaux107*yaux114*yaux115*yaux96;
	out[6] = 17.77777777777778*yaux1*yaux117 + 17.77777777777778*yaux175*yaux196 + 17.77777777777778*yaux127*(yaux124 + yaux203) - 0.977777777777778*yaux178*yaux194*yaux202*yaux206 - 0.977777777777778*yaux178*yaux194*yaux208*yaux209 + 17.77777777777778*yaux129*(yaux132 + yaux211) + 17.77777777777778*(yaux119 + yaux120 + yaux198 + yaux199)*yaux96 - 0.977777777777778*yaux115*yaux123*yaux48*yaux99 - 0.977777777777778*yaux115*yaux130*yaux50*yaux99;
	out[7] = -0.977777777777778*yaux1*yaux103*yaux114*yaux115 - 0.977777777777778*yaux175*yaux182*yaux193*yaux194 - 0.977777777777778*yaux182*yaux194*yaux208*yaux209 + 17.77777777777778*yaux206*yaux222 + 17.77777777777778*yaux127*(yaux119 + yaux145 + yaux198 + yaux224) + 17.77777777777778*yaux129*(yaux150 + yaux229) + 17.77777777777778*yaux143*yaux48 - 0.977777777777778*yaux103*yaux115*yaux130*yaux50 + 17.77777777777778*(yaux139 + yaux218)*yaux96;
	out[8] = -0.977777777777778*yaux1*yaux107*yaux114*yaux115 - 0.977777777777778*yaux175*yaux186*yaux193*yaux194 - 0.977777777777778*yaux186*yaux194*yaux202*yaux206 + 17.77777777777778*yaux127*(yaux160 + yaux239) + 17.77777777777778*yaux208*yaux245 + 17.77777777777778*yaux129*(yaux119 + yaux168 + yaux198 + yaux247) - 0.977777777777778*yaux107*yaux115*yaux123*yaux48 + 17.77777777777778*yaux166*yaux50 + 17.77777777777778*(yaux157 + yaux236)*yaux96;
	out[9] = -0.977777777777778*yaux127*yaux178*yaux194*yaux202 - 0.977777777777778*yaux129*yaux178*yaux194*yaux209 + 17.77777777777778*yaux254*yaux275 + 17.77777777777778*yaux175*(yaux198 + yaux199 + yaux277 + yaux278) + 17.77777777777778*yaux206*(yaux203 + yaux282) - 0.977777777777778*yaux257*yaux273*yaux281*yaux285 - 0.977777777777778*yaux257*yaux273*yaux287*yaux288 + 17.77777777777778*yaux208*(yaux211 + yaux290) + 17.77777777777778*yaux196*yaux96;
	out[10] = -0.977777777777778*yaux129*yaux182*yaux194*yaux209 + 17.77777777777778*yaux127*yaux222 - 0.977777777777778*yaux254*yaux261*yaux272*yaux273 - 0.977777777777778*yaux261*yaux273*yaux287*yaux288 + 17.77777777777778*yaux175*(yaux218 + yaux297) + 17.77777777777778*yaux285*yaux301 + 17.77777777777778*yaux206*(yaux198 + yaux224 + yaux277 + yaux303) + 17.77777777777778*yaux208*(yaux229 + yaux308) - 0.977777777777778*yaux182*yaux193*yaux194*yaux96;
	out[11] = -0.977777777777778*yaux127*yaux186*yaux194*yaux202 + 17.77777777777778*yaux129*yaux245 - 0.977777777777778*yaux254*yaux265*yaux272*yaux273 - 0.977777777777778*yaux265*yaux273*yaux281*yaux285 + 17.77777777777778*yaux175*(yaux236 + yaux315) + 17.77777777777778*yaux206*(yaux239 + yaux318) + 17.77777777777778*yaux287*yaux324 + 17.77777777777778*yaux208*(yaux198 + yaux247 + yaux277 + yaux326) - 0.977777777777778*yaux186*yaux193*yaux194*yaux96;
	out[12] = 17.77777777777778*yaux175*yaux275 - 0.977777777777778*yaux206*yaux257*yaux273*yaux281 - 0.977777777777778*yaux208*yaux257*yaux273*yaux288 + 17.77777777777778*yaux333*yaux354 + 17.77777777777778*yaux254*(yaux277 + yaux278 + yaux356 + yaux357) + 17.77777777777778*yaux285*(yaux282 + yaux361) - 0.977777777777778*yaux336*yaux352*yaux360*yaux364 - 0.977777777777778*yaux336*yaux352*yaux366*yaux367 + 17.77777777777778*yaux287*(yaux290 + yaux369);
	out[13] = -0.977777777777778*yaux175*yaux261*yaux272*yaux273 - 0.977777777777778*yaux208*yaux261*yaux273*yaux288 + 17.77777777777778*yaux206*yaux301 - 0.977777777777778*yaux333*yaux340*yaux351*yaux352 - 0.977777777777778*yaux340*yaux352*yaux366*yaux367 + 17.77777777777778*yaux254*(yaux297 + yaux376) + 17.77777777777778*yaux364*yaux380 + 17.77777777777778*yaux285*(yaux277 + yaux303 + yaux356 + yaux382) + 17.77777777777778*yaux287*(yaux308 + yaux387);
	out[14] = -0.977777777777778*yaux175*yaux265*yaux272*yaux273 - 0.977777777777778*yaux206*yaux265*yaux273*yaux281 + 17.77777777777778*yaux208*yaux324 - 0.977777777777778*yaux333*yaux344*yaux351*yaux352 - 0.977777777777778*yaux344*yaux352*yaux360*yaux364 + 17.77777777777778*yaux254*(yaux315 + yaux394) + 17.77777777777778*yaux285*(yaux318 + yaux397) + 17.77777777777778*yaux366*yaux403 + 17.77777777777778*yaux287*(yaux277 + yaux326 + yaux356 + yaux405);
	out[15] = 17.77777777777778*yaux254*yaux354 - 0.977777777777778*yaux285*yaux336*yaux352*yaux360 - 0.977777777777778*yaux287*yaux336*yaux352*yaux367 + 17.77777777777778*yaux412*yaux433 + 17.77777777777778*yaux333*(yaux356 + yaux357 + yaux435 + yaux436) + 17.77777777777778*yaux364*(yaux361 + yaux440) - 0.977777777777778*yaux415*yaux431*yaux439*yaux443 - 0.977777777777778*yaux415*yaux431*yaux445*yaux446 + 17.77777777777778*yaux366*(yaux369 + yaux448);
	out[16] = -0.977777777777778*yaux254*yaux340*yaux351*yaux352 - 0.977777777777778*yaux287*yaux340*yaux352*yaux367 + 17.77777777777778*yaux285*yaux380 - 0.977777777777778*yaux412*yaux419*yaux430*yaux431 - 0.977777777777778*yaux419*yaux431*yaux445*yaux446 + 17.77777777777778*yaux333*(yaux376 + yaux455) + 17.77777777777778*yaux443*yaux459 + 17.77777777777778*yaux364*(yaux356 + yaux382 + yaux435 + yaux461) + 17.77777777777778*yaux366*(yaux387 + yaux466);
	out[17] = -0.977777777777778*yaux254*yaux344*yaux351*yaux352 - 0.977777777777778*yaux285*yaux344*yaux352*yaux360 + 17.77777777777778*yaux287*yaux403 - 0.977777777777778*yaux412*yaux423*yaux430*yaux431 - 0.977777777777778*yaux423*yaux431*yaux439*yaux443 + 17.77777777777778*yaux333*(yaux394 + yaux473) + 17.77777777777778*yaux364*(yaux397 + yaux476) + 17.77777777777778*yaux445*yaux482 + 17.77777777777778*yaux366*(yaux356 + yaux405 + yaux435 + yaux484);
	out[18] = 17.77777777777778*yaux333*yaux433 - 0.977777777777778*yaux364*yaux415*yaux431*yaux439 - 0.977777777777778*yaux366*yaux415*yaux431*yaux446 + 17.77777777777778*yaux412*(yaux435 + yaux436 + yaux507 + 0.05500000000000001*yaux493*yaux508*yaux509) + 17.77777777777778*yaux443*(yaux440 + 0.05500000000000001*yaux493*yaux509*yaux513) + 17.77777777777778*yaux445*(yaux448 + 0.05500000000000001*yaux493*yaux509*yaux517);
	out[19] = -0.977777777777778*yaux333*yaux419*yaux430*yaux431 - 0.977777777777778*yaux366*yaux419*yaux431*yaux446 + 17.77777777777778*yaux364*yaux459 + 17.77777777777778*yaux412*(yaux455 + 0.05500000000000001*yaux497*yaux508*yaux509) + 17.77777777777778*yaux443*(yaux435 + yaux461 + yaux507 + 0.05500000000000001*yaux497*yaux509*yaux513) + 17.77777777777778*yaux445*(yaux466 + 0.05500000000000001*yaux497*yaux509*yaux517);
	out[20] = -0.977777777777778*yaux333*yaux423*yaux430*yaux431 - 0.977777777777778*yaux364*yaux423*yaux431*yaux439 + 17.77777777777778*yaux366*yaux482 + 17.77777777777778*yaux412*(yaux473 + 0.05500000000000001*yaux501*yaux508*yaux509) + 17.77777777777778*yaux443*(yaux476 + 0.05500000000000001*yaux501*yaux509*yaux513) + 17.77777777777778*yaux445*(yaux435 + yaux484 + yaux507 + 0.05500000000000001*yaux501*yaux509*yaux517);
	out[21] = 0.977777777777778*yaux509*(1 * yaux443*yaux493*yaux497 + 1 * yaux445*yaux493*yaux501 + 1 * yaux412*(1 * yaux494 + yaux553));
	out[22] = 0.977777777777778*yaux509*(1 * yaux412*yaux493*yaux497 + 1 * yaux445*yaux497*yaux501 + 1 * yaux443*(1 * yaux498 + yaux553));
	out[23] = -14.222222222222223*yaux509*(-0.06875*yaux412*yaux493*yaux501 - 0.06875*yaux443*yaux497*yaux501 + 1 * yaux445*(-0.06875*yaux502 + 1 * yaux503*yaux552));
	out[24] = vec[0];
	out[25] = vec[1];
	out[26] = vec[2];
	out[27] = vec[3];
	out[28] = vec[4];
	out[29] = vec[5];
	out[30] = vec[6];
	out[31] = vec[7];
	out[32] = vec[8];
	out[33] = vec[9];
	out[34] = vec[10];
	out[35] = vec[11];
	out[36] = vec[12];
	out[37] = vec[13];
	out[38] = vec[14];
	out[39] = vec[15];
	out[40] = vec[16];
	out[41] = vec[17];
	out[42] = vec[18];
	out[43] = vec[19];
	out[44] = vec[20];
}
/** Jacobian df/du multiplied by vector vec, i.e. (df/du)^T*vec or vec^T*(df/du) **/
void dfdu_vec(typeRNum *out, ctypeRNum t, ctypeRNum *x, ctypeRNum *vec, ctypeRNum *u, ctypeRNum *p, typeUSERPARAM *userparam)
{
	out[0] = vec[21];
	out[1] = vec[22];
	out[2] = vec[23];
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

	yaux11 = x[21];

	out[0] = param[2] * (POW2(u[0]) + POW2(u[1]) + POW2(u[2]))
		+ param[1] * (56.25 - 15 * yaux11 + yaux11 * yaux11 + POW2(x[22]) + POW2(x[23]))
		+ param[0] * (POW2(x[24]) + POW2(x[25]) + POW2(x[26]) + POW2(x[27]) + POW2(x[28]) + POW2(x[29])
			+ POW2(x[30]) + POW2(x[31]) + POW2(x[32]) + POW2(x[33]) + POW2(x[34]) + POW2(x[35]) + POW2(x[36])
			+ POW2(x[37]) + POW2(x[38]) + POW2(x[39]) + POW2(x[40]) + POW2(x[41]) + POW2(x[42]) + POW2(x[43]) + POW2(x[44]));
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
	out[21] = 2 * yaux1*(-7.5 + x[21]);
	out[22] = 2 * yaux1*x[22];
	out[23] = 2 * yaux1*x[23];
	out[24] = 2 * yaux9*x[24];
	out[25] = 2 * yaux9*x[25];
	out[26] = 2 * yaux9*x[26];
	out[27] = 2 * yaux9*x[27];
	out[28] = 2 * yaux9*x[28];
	out[29] = 2 * yaux9*x[29];
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

	out[0] = param[3] * (POW2(-7.5 + x[21]) + POW2(x[22]) + POW2(x[23]));
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
	out[21] = 2 * yaux1*(-7.5 + x[21]);
	out[22] = 2 * yaux1*x[22];
	out[23] = 2 * yaux1*x[23];
	out[24] = 0;
	out[25] = 0;
	out[26] = 0;
	out[27] = 0;
	out[28] = 0;
	out[29] = 0;
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
