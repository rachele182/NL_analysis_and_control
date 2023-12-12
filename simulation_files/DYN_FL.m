%% DYNAMIC FEEDBACK LINEARIZATION
%% input T is delayed 

run init.m
n_sys = 8;

% Variables
syms x dx theta dtheta y dy T phi real 
% Parameters
syms mu_x mu_y l M grav J dphi dT real 

q = [x dx theta dtheta y dy phi T]';
% Inputs
u = [dT,dphi];

f = [   
	dx;
    (M^-1)*(-mu_x*dx) + T*M^-1*(-sin(theta+phi));
    dtheta;
    T*J^-1*(-l*sin(phi));
    dy;
    T*M^-1*(cos(theta+phi)) + M^-1*( - mu_y*dy -M*grav);
    0;
    0]; 
   
g1 = [0 0 0 0 0 0 0 1]';
g2 = [	0;
  		0;
        0;
        0;
        0;
        0;
        1;
		0];
g = [g1 , g2]; 

% outputs
% h1 = (x^2+y^2)/2;
h1 = y;
h2 = theta;
h = [h1; h2];
%% Check conditions to have complete feedback linearization for MIMO systems

% CONDITION : singularity + involutive

k=1;
tol = n_sys;
n_in = size(g,2);
gamma = {g1, g2};
gamma_mat = [g1, g2];
cond = true;
while k < n_sys
	[gamma, gamma_mat] = filtr_feed(gamma, gamma_mat, f, q, n_in);

	gamma_mat = simplify(gamma_mat);
	r = rank(gamma_mat);

	% singularity check
	gamma_mat_num = subs(gamma_mat,{x dx theta dtheta y dy phi},{zeros(1,7)});
	r_num = rank(gamma_mat_num);
	if r ~= r_num
		cond = false;
		fprintf("gamma singular on step: %i\n", k)
	end

	% involutivity check
    if k>1
        invol_mat=[];
        for index1 = 1:(size(gamma,2)-1)
            for index2 = (index1+1):size(gamma,2)
                lie_tmp = lie_b(gamma{index1},gamma{index2},q);
                invol_mat = [invol_mat, lie_tmp];
            end
        end
        invol_mat_num = subs(invol_mat,{x dx theta dtheta y dy phi},{zeros(1,7)});
        r_invol = rank([gamma_mat_num, invol_mat_num]);
        if r_invol ~= r_num
            cond = false;
            fprintf("gamma not involutive on step: %i\n\n", k)
        end
    end

	k = k+1;
end
% rank condition
if cond == true
	r = rank(gamma_mat);
	if r == n_sys
		fprintf("complete linearization is possible")
	else
		fprintf("go on with partial linearization")
	end
end




%% Feedback Linearization Input-Output

% derivative series
lg1lfh1 = dederni(h1, g1, q);
lg2lfh1 = dederni(h1, g2, q);
lg1lfh2 = dederni(h2, g1, q);
lg2lfh2 = dederni(h2, g2, q);

index = 1;
y1{1} = dederni(h1, f, q);
while lg1lfh1 == 0 && lg2lfh1 == 0 
    y1{index+1} = dederni(y1{index}, f, q);
    lg1lfh1 = dederni(y1{index}, g1, q);
    lg2lfh1 = dederni(y1{index}, g2, q);
    index = index +1;
end

index = 1;
y2{1} = dederni(h2, f, q);
while lg1lfh2 == 0 && lg2lfh2 == 0
    y2{index+1} = dederni(y2{index}, f, q);
    lg1lfh2 = dederni(y2{index}, g1, q);
    lg2lfh2 = dederni(y2{index}, g2, q);
    index = index +1;
end

r1 = size(y1,2);
r2 = size(y2,2);

Gamma = [y1{r1}; y2{r2}];
E = [lg1lfh1, lg2lfh1;
     lg1lfh2, lg2lfh2]

E_num = subs(E, {x dx theta dtheta y dy phi},{zeros(1,7)})

r_tot = r1 + r2;
if det(E_num) ~= 0
    fprintf('the relative degree is: %i\n', r_tot)
    if r_tot == n_sys
        fprintf('the system is completely linearizable')
    else
        fprintf('the system is partially linearizable')
    end
else
    fprintf('It is not possible to linearize the system')
end

%% Variable change
csi1 = [h1];
for index = 1:r1-1
    csi1 = [csi1; y1{index}];
end
csi2 = [h2];
for index = 1:r2-1
    csi2 = [csi2; y2{index}];
end
csi = [csi1; csi2]
zita = [x; dx]
Z = [csi;zita]
jac = jacobian(Z,q);
jac_num = subs(jac, {x dx theta dtheta y dy phi},{zeros(1,7)});
d = det(jac_num) %check variable change

%% Normal form
A1 = [0 1 0; 0 0 1; 0 0 0];
A2 = [0 1 0; 0 0 1; 0 0 0];
A = blkdiag(A1,A2);
B1 = [0 0 1]';
B2 = [0 0 1]';
B = blkdiag(B1,B2);
C1 = [1 0 0];
C2 = [1 0 0];
C = blkdiag(C1,C2);

zita_dot = jacobian(zita,q)*f; 
%da sostiture a q le nuove variabili csi,zita
zita_dot = [zita(2); -M^-1*mu_x*zita(2) + T*M^-1*(-sin(csi(4) + asin(-csi(6)*J/(T*l))))];

% zero dinamica
syms s1 s2 v1 v2 real
s = [s1;s2];
s_dot = [   s2;
            -M^-1*mu_x*s2];
V = 0.5*s1^2 + 0.5*s2^2;
V_dot = jacobian(V, s)*s_dot;

% control
v = [v1; v2];
u = -(E^-1)*Gamma + (E^-1)*v;
% v= K*csi; 2x6

%% Functions

function [gamma,gamma_mat] = filtr_feed(gamma_in, gamma_mat_in, f, q, n_in)
    gamma_mat = gamma_mat_in;
    gamma = gamma_in;
    r = rank(gamma_mat);
	index = size(gamma_in,2) - n_in + 1;
	while index <= size(gamma_in,2)
		lie_tmp = lie_b(f,gamma_in{index},q);
		gamma_mat = [gamma_mat, lie_tmp];
		gamma{end+1} = lie_tmp;
		index = index + 1;
	end
end


function Lgf = lie_b(g,f,q)
    Lgf = jacobian(f,q)*g - jacobian(g,q)*f;
end

function der = dederni(h, f, q)
    der = jacobian(h,q)*f;
end