%--- University of Washington, Department of Aeronautics & Astronautics ---
%---------- Advanced Dynamics, Validation & Control Research Lab ----------
%
% Numerical integrator for the RABT algorithm (used with rabt.m)
%
% Author: Peter Uth
% Created: July 2017
%--------------------------------------------------------------------------
function [xdot,g,y] = rabt_iter(x,u,p)
global fun step type

if strcmp(type,'stepsize') || ~exist('type','var')
    h = p(end);     % for varying step size (default)
elseif strcmp(type,'tolerance')
    h = 0.1;        % for varying zero tolerance
end

xdot = fun(x,u,p);

% let state be initial conditions
xt = p(1:end-1);

% preallocate sizes
k1 = zeros(numel(xt),step);
k2 = zeros(numel(xt),step);
k3 = zeros(numel(xt),step);
k4 = zeros(numel(xt),step);
xtiter = zeros(numel(xt),step);

% 4th order Runge Kutta
for j = 1:step
  k1(:,j) = fun(xt,u,p);
  k2(:,j) = fun(xt+0.5*h*k1(:,j),u,p); 
  k3(:,j) = fun(xt+0.5*h*k2(:,j),u,p);
  k4(:,j) = fun(xt+h*k2(:,j),u,p);
  xtiter(:,j) = xt+(h/6)*(k1(:,j)+2*k2(:,j)+2*k3(:,j)+k4(:,j));
  xt = xtiter(:,j);
end
xtend = xtiter(:,j);

if strcmp(type,'stepsize') || ~exist('type','var')
    g = xtend - x;                      % for varying step size (default)
elseif strcmp(type,'tolerance')
    g = xtend - x - [p(end);p(end)];    % for varying zero tolerance
end

y = [];
   
end
