%    begin                : November 2020
%    authors              : Rachele Nebbia Colomba, Chiara Sammarco, Giorgio Simonini
%    copyright            : Dipartimento di Ingegneria dell`Informazione (DII) Universita´ di pisa    
%    email                : rachelenebbia <at> gmail <dot> com

%%Description: Init file for Non-linear analysis and control of a rocket Falcone 9 system
%              here we define the parameters and initial conditions to run the simulations

%Clear workspace
clc; close all; clear all;

display('Falcon9 V1.1 loading parameters...')
% RocketSystem Parameters
g_par = 9.81;                       % m/s^2
M_par = 5e5;                        % Kg
l_par = 34;                         % m semi-axis
J_par = (M_par*4*l_par^2)/12;       % Nm            (assumption: uniform mass)
d_par = 3.66;                       % m diameter

% Fdrag = 0.5*rho*cd*A*v^2; 

cdy_par = 0.20;                             % drag coefficient y axis
cdx_par = 1;                                % drag coefficient y axis
Ax_par= d_par*2*l_par;                      % impact area
Ay_par = pi*(d_par/2)^2;                    % impact area
rho_par = 1.225;                            % Kg/m^3 air density   
mu_x_par = 0.5*rho_par*cdx_par*Ax_par;      % N*s/m viscous friction coefficient
mu_y_par= 0.5*rho_par*cdy_par*Ay_par;       % N*s/m viscous friction 
T_par =  M_par*g_par;                   	% Nm thrusther constant propulsion 
% T_par = 5885; 
n_sys = 7;                                  % system degree


%% Initial conditions dyn feedback control
T0 = M_par*g_par;                          % Nm thrusther propulsion 
x0 = 2;                                    % m initial offset x-axis
y0 = -5;                                   % m initial offset y-axis
theta0 = deg2rad(15);                      % ° initial offset heading angle
muy0 = 0.1;                                % N*s/m initial viscous friction 
M0 = 0.5*M_par;                            % kg uncertain mass 
J0 = 0.8*J_par;                            % kg*m^3 uncertain moment of inertia 

display('Loading done, Enjoy!')