%% Non linear System 
 run init.m %loading parameters

%% BASIC ANALYSIS
% Case 0: State: q = [x dx theta dtheta y dy] 
         % Inputs = [T phi]
         % Outputs = [x y theta]
syms x dx theta dtheta y dy T phi real
% Parameters
syms mu_x mu_y l M g J real
% Equations
x_dot = dx;
dx_dot = M^-1*(-T*sin(theta+phi) - mu_x*dx);
theta_dot = dtheta;
dtheta_dot = J^-1*(-T*l*sin(phi));
y_dot = dy;
dy_dot = M^-1*(T*cos(theta+phi) - mu_y*dy - M*g);

q = [x dx theta dtheta y dy]';
f = [x_dot dx_dot theta_dot dtheta_dot y_dot dy_dot]';
u = [T phi]';
h = [x; theta; y];
%Linearization around the equilibrium point qeq = zeros(6,1); ueq = [T0 0]
syms T0 real
A = jacobian(f,q');
B = jacobian(f,u');
C = jacobian(h,q'); 
A = subs(A, {x dx theta dtheta y dy phi T},{0 0 0 0 0 0 0 T_par});
B = subs(B, {x dx theta dtheta y dy phi T},{0 0 0 0 0 0 0 T_par});
C = subs(C, {x dx theta dtheta y dy phi T},{0 0 0 0 0 0 0 T_par});

% check controllability of approximate linear system with inputs (T,phi)
Ac = subs(A,{x dx theta dtheta y dy phi M mu_x mu_y},{0 0 0 0 0 0 0 M_par 0.5 0.5});
Bc = subs(B, {x dx theta dtheta y dy phi M J l},{0 0 0 0 0 0 0 M_par J_par l_par});
R = ctrb(Ac,Bc);
r = rank(R); 
% check observability linear system with outputs (x,theta,y)
O = obsv(Ac,C);
o = rank(O);

%In this case the approximated linear system is completely controllable and
%observable (sufficient condition for local observability and
%controllability non linear system)

%% CONTROLLABILITY ANALYSIS
% Case 1: No control of propulsion, T constant.  
% State: q = [x x_dot theta theta_dot y y_dot] 
         % Inputs = [phi]
         
u1 = [phi];
A1 = jacobian(f,q');
B1 = jacobian(f,u1');
A1 = subs(A1,{x dx theta dtheta y dy phi},{zeros(1,7)});
B1 = subs(B1, {x dx theta dtheta y dy phi},{zeros(1,7)});

% check controllability of approximate linear system with inputs (phi)
Ac1 = subs(A1,{x dx theta dtheta y dy phi M T mu_x mu_y},{0 0 0 0 0 0 0 M_par T_par 0.5 0.5});
Bc1 = subs(B1, {x dx theta dtheta y dy phi M T J l mu_x mu_y},{0 0 0 0 0 0 0 M_par T_par J_par l_par 0.5 0.5});
R1 = ctrb(Ac1,Bc1);
r1 = rank(R1);
% --> in this case no complete controllability of approximated linear
% system (in particular y,ydot not reachable)



%% OBSERVABILITY ANALYSIS
% Case 2: Measurements: y theta theta_dot
h2 = [theta;theta_dot;y]; 
C2 = jacobian(h2,q');
%Linearization around qeq = 0 , ueq = [0,T0]
C2 = subs(C2, {x dx theta dtheta y dy phi},{zeros(1,7)});
O2 = obsv(Ac,C2);
o2 = rank(O2);
%Approx. linear system not completely observable, o2 = 4 --> x,xdot not observable 

%% Case 3: Measurements:
h3 = [0.5*(x^2+y^2); theta];
C3 = jacobian(h3,q');
%Linearization around qeq = 0 , ueq = [0,T0]
C3 = subs(C3, {x dx theta dtheta y dy phi},{zeros(1,7)});
O3 = obsv(Ac,C3);
o3 = rank(O3);
%Approx. linear system not completely observable, o3 =2 --> only
%theta,theta_dot observable



