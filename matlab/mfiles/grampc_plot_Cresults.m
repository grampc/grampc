function vec = grampc_plot_Cresults(filepath,figNr,OCP,decimation)
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

if nargin <1 || isempty(filepath)
    filepath= [cd '\res\'];
end
if nargin <2 || isempty(figNr)
    figNr = 200;
end
if nargin < 3 || isempty(OCP)
    OCP = 0;
end
if nargin < 4 || isempty(decimation)
    decimation = 1;
end

vec.x = load([filepath 'xvec.txt'])';
vec.u = load([filepath 'uvec.txt'])';
vec.t = load([filepath 'tvec.txt'])';
vec.status = load([filepath 'status.txt'])';
vec.iter = load([filepath 'itervec.txt'])';

if OCP ==0
    vec.J = load([filepath 'Jvec.txt'])';
    vec.p = load([filepath 'pvec.txt'])';
    vec.Ncfct = load([filepath 'Ncfctvec.txt'])';
    vec.Npen = load([filepath 'Npenvec.txt'])';
    vec.Thor = load([filepath 'Thorvec.txt'])';
else
    vec.J = load([filepath 'Jiter.txt']);
    vec.p = load([filepath 'piter.txt'])';
    vec.Ncfct = load([filepath 'Ncfctiter.txt'])';
    vec.Npen = load([filepath 'Npeniter.txt'])';
    vec.Thor = load([filepath 'Thoriter.txt'])';
    vec.cfct = load([filepath 'cfctvec.txt'])';
    vec.lag = load([filepath 'lagvec.txt'])';
    vec.multiter = load([filepath 'multitervec.txt'])';
    vec.pen = load([filepath 'penvec.txt'])';
    auglagiter = 1:length(vec.multiter);
end


figure(figNr)
sh(1) = subplot(2,3,1);
plot(vec.t(1:decimation:end),vec.x(:,1:decimation:end)),xlabel('time'),title('States')

sh(3) = subplot(2,3,3);
plot(vec.t(1:decimation:end),vec.u(:,1:decimation:end)),xlabel('time'),title('Controls')


if OCP ==0
sh(2) = subplot(2,3,2);
plot(vec.t(1:decimation:end),vec.J(:,1:decimation:end)),xlabel('time'),title('Costs')

sh(4) = subplot(2,3,4);
plot(vec.t(1:decimation:end),vec.Ncfct(1:decimation:end)),xlabel('time'),title('||cfct||')

sh(5) = subplot(2,3,5);
plot(vec.t(1:decimation:end),vec.Thor(1:decimation:end)),xlabel('time'),title('Thor')

sh(6) = subplot(2,3,6);
plot(vec.t(1:decimation:end),vec.Npen(1:decimation:end)),xlabel('time'),title('||pen||')
    
    
else
sh2(1) = subplot(2,3,2);
plot(auglagiter,vec.J(:,1:decimation:end)),xlabel('aug lag iterations'),title('Costs')

sh(2) = subplot(2,4,5);
plot(vec.t(1:decimation:end),vec.cfct(:,1:decimation:end)),xlabel('time'),title('Constraints')

sh(4) = subplot(2,4,6);
plot(vec.t(1:decimation:end),vec.lag(:,1:decimation:end)),xlabel('time'),title('Lagrangian Multiplier')

sh(5) = subplot(2,4,7);
plot(vec.t(1:decimation:end),vec.pen(:,1:decimation:end)),xlabel('time'),title('Penalty Parameter')

sh2(2) = subplot(2,4,8);
stem(auglagiter,vec.iter(:,1:decimation:end)),xlabel('aug lag iterations'),title('gradient iterations')

linkaxes(sh2,'x')
xlim(sh2,[auglagiter(1)-0.25 auglagiter(end)+0.25])
end

linkaxes(sh,'x')
xlim(sh,[vec.t(1) vec.t(end)])
end