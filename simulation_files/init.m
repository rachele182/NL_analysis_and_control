%% Initialization file; 
% 
% clc; close all; clear all;

display('Falcon9 V1.1 loading parameters...')
% System Parameters
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
rho_par = 1.225;                            % Kg/m^3    
mu_x_par = 0.5*rho_par*cdx_par*Ax_par;      % N*s/m viscous friction coefficient
mu_y_par= 0.5*rho_par*cdy_par*Ay_par;       % N*s/m viscous friction 
T_par =  M_par*g_par;                   	% thrusther constant propulsion 
% T_par = 5885; 
n_sys = 7;


%% Initial conditions dyn feedback 
T0 = M_par*g_par;
x0 = 2; %[m]
y0 = -5; %[m]
theta0 = deg2rad(15); %[gradi]

muy0 = 0.1;
M0 = 0.5*M_par;
J0 = 0.8*J_par;
display('Loading done, Enjoy!')