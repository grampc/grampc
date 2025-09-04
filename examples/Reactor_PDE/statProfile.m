% N=11;
% p=[2,0,1,1,0.2]
% [x0,u0]     = statProfile(p,1,N);
% [xdes,udes] = statProfile(p,2,N);

function [xStat,uStat] = statProfile(p,yStat,N)
%
% This file is part of GRAMPC - (https://github.com/grampc/grampc)
%
% GRAMPC -- A software framework for embedded nonlinear model predictive
% control using a gradient-based augmented Lagrangian approach
%
% Copyright 2014-2025 by Knut Graichen, Andreas Voelz, Thore Wietzke,
% Tobias Englert (<v2.3), Felix Mesmer (<v2.3), Soenke Rhein (<v2.3),
% Bartosz Kaepernick (<v2.0), Tilman Utz (<v2.0).
% All rights reserved.
%
% GRAMPC is distributed under the BSD-3-Clause license, see LICENSE.txt
%
% File: statProfile.m
% Authors: Bartosz Kaepernick, Knut Graichen, Tilman Utz
% Date: February 2014
% Version: v1.0
%
% MATLAB/Simulink file to determine initial conditions for the reactor 
% example of GRAMPC.
%

z       = linspace(0,1,N);
solinit = bvpinit(z,[yStat,0.0]);

opt    = bvpset('FJacobian',@fjac);
solbvp = bvp4c(@ode,@bc,solinit,opt,p,yStat);

sol = deval(solbvp,z);

xStat = sol(1,:);
uStat = 0*sol(1,end) + sol(2,end);

end


function out = ode(z,y,p,yStat)

q0  = p(1);
q1  = p(2);
nu  = p(3);
r0  = p(4);
r1  = p(5);

out = [ y(2);
	      (nu*y(2) - q1*y(2)*y(2) - (r0 + r1*y(1))*y(1))/(q0 + q1*y(1)); ];

end


function out = bc(ya,yb,p,yStat)

out = [ ya(1)-yStat;
	      ya(2); ];

end


function out = fjac(z,y,p,yStat)

q0  = p(1);
q1  = p(2);
nu  = p(3);
r0  = p(4);
r1  = p(5);

out(1,1) = 0.0; 
out(1,2) = 1.0;
out(2,1) = -(r0 + 2 * r1 * y(1)) / (q0 + q1 * y(1)) + (q1 * y(2) ^ 2 - nu * y(2) + r0 * y(1) + r1 * y(1) ^ 2) / (q0 + q1 * y(1)) ^ 2 * q1;
out(2,2) = -(2 * q1 * y(2) - nu) / (q0 + q1 * y(1));

end
