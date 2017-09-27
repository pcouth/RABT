%--- University of Washington, Department of Aeronautics & Astronautics ---
%---------- Advanced Dynamics, Validation & Control Research Lab ----------
%
% A script to test the RABT algorithm on 3 examples
%--------------------------------------------------------------------------
close all; clear all; clc

%% Example 1: Van der Pol Oscillator
F1 = @(x,u,p) [-x(2);x(1)+(x(1)^2-1)*x(2)];
plimit1 = [-3 3;-3 3;0.1 1];     
steps1 = 100;                   
equil1 = [0;0];                  
ztol1 = 0.001; 

tic
basin1 = rabt(F1,equil1,plimit1,steps1,ztol1);
toc

open('ex1_plane.fig')
hold on
scatter(basin1(1,:),basin1(2,:),15,'b','filled')
box on, xlim([-3 3]), ylim([-3 3])

%% Example 2
F2 = @(x,u,p) [-0.42*x(1)-1.05*x(2)-2.3*x(1)^2-0.5*x(1)*x(2)-x(1)^3;1.98*x(1)+x(1)*x(2)];
plimit2 = [-1.5 1.5;-2.5 1.5;0.1 1]; 
steps2 = 100;                   
equil2 = [0;0];                  
ztol2 = 0.001;

tic
basin2 = rabt(F2,equil2,plimit2,steps2,ztol2);
toc

open('ex2_plane.fig')
hold on
scatter(basin2(1,:),basin2(2,:),15,'b','filled')
box on, xlim([-1.5 1.5]), ylim([-2.5 1.5])

%% Example 3
F3 = @(x,u,p) [-0.5*x(1)-x(2)-x(1)^3;-0.5*x(1)-x(1)^2];
plimit3 = [-2 2;-2 2;0.1 1];
steps3 = 200;
equil3 = [-0.5;0.375];
ztol3 = 0.001;

tic
basin3 = rabt(F3,equil3,plimit3,steps3,ztol3);
toc

open('ex3_plane.fig')
hold on
scatter(basin3(1,:),basin3(2,:),15,'b','filled')
box on, xlim([-2 2]), ylim([-2 2])
