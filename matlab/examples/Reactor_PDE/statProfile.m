function [xStat,uStat] = statProfile(p,yStat,N)
  %
  % This file is part of GRAMPC.
  %
  % GRAMPC - a gradient-based MPC software for real-time applications
  %
  % Copyright (C) 2014 by Bartosz Kaepernick, Knut Graichen, Tilman Utz
  % Developed at the Institute of Measurement, Control, and
  % Microtechnology, University of Ulm. All rights reserved.
  %
  % GRAMPC is free software: you can redistribute it and/or modify
  % it under the terms of the GNU Lesser General Public License as 
  % published by the Free Software Foundation, either version 3 of 
  % the License, or (at your option) any later version.
  %
  % GRAMPC is distributed in the hope that it will be useful,
  % but WITHOUT ANY WARRANTY; without even the implied warranty of
  % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  % GNU Lesser General Public License for more details.
  %
  % You should have received a copy of the GNU Lesser General Public 
  % License along with GRAMPC. If not, see <http://www.gnu.org/licenses/>.
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
uStat = sol(1,end) + sol(2,end);

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
