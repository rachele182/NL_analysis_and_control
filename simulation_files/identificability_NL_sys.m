run init.m


% creation of system in state affine control form
% state vector
syms x dx theta dtheta y dy  phi real 
% Parameters
syms mu_x mu_y l M grav J T dphi real

q = [x dx theta dtheta y dy phi mu_x mu_y]';
u = [T, dphi]; %inputs
n_sys = size(q,1);

f = [   
	dx;
    (M^-1)*(-mu_x*dx);
    dtheta;
    0;
    dy;
    M^-1*( -mu_y*dy -M*grav);
    0;
    0;
    0]; 

g1 = [
	0;
    -(M^-1)*sin(theta+phi);
    0;
	-(J^-1)*(l*sin(phi));
    0;
	(M^-1)*cos(theta+phi);
    0;
    0;
    0];

g2 = [
	0;
	0;
	0;
	0;
	0;
	0;
	1;
    0;
    0];


%% 1° case:
% definition of outputs  
h1 = 0.5*(x^2+y^2);
h2 = theta;
%Initialization of steps
dh1 = jacobian(h1,q); 
dh2 = jacobian(h2,q);
delta = {f, g1, g2};
omega = {dh1, dh2};
omega_mat = [dh1; dh2];
q0 = [0; 0; 0; 0; 0; 0; 0; 0; 0]; 
r{1} = rank(omega_mat);


%% iterative filtration method
k=1;
tol = n_sys;
omega_mat_prec = zeros(size(omega_mat));

while k<tol
    [omega,omega_mat,r{k+1}] = filtration(omega, omega_mat, delta, q);
    omega_mat_num = subs(omega_mat, {q(1) q(2) q(3) q(4) q(5) q(6) q(7)}, {zeros(1,7)});
    r_num = rank(omega_mat_num);
    omega_mat_prec_num = subs(omega_mat_prec, {q(1) q(2) q(3) q(4) q(5) q(6) q(7)}, {zeros(1,7)});
    r_num_prec = rank(omega_mat_prec_num);
    if r_num == r{k+1} && r_num_prec == r{k} %check if deltak deltak+1 not singular
        if r{k+1} == r{k}    %check stop condition
           k = tol;
        end 
    end
    k = k+1;
    omega_mat_prec = omega_mat;
end
omega_mat = simplify(omega_mat);



%% Utils
%Filtration
function [omega,omega_mat,r] = filtration(omega_in,omega_mat_in,delta_in,q)
    omega_mat = omega_mat_in;
    omega = omega_in;
    r = rank(omega_mat);
    for index1 = 1:size(omega_in,2)
        for index2 = 1:size(delta_in,2)
            lie_tmp = lie_b(omega_in{index1},delta_in{index2},q);
            omega_mat = [omega_mat; lie_tmp];
			r = rank(omega_mat);
            %r_new = rank(omega_mat);
            if is_zero(lie_tmp) == 0                     % check if new vector is null
                omega{end+1} = lie_tmp;
            else
                omega_mat = omega_mat(1:end-1,:);
            end
        end
    end
end



function Lvw = lie_b(w,v,q)
    Lvw = v'*jacobian(w',q)' + w*jacobian(v,q);
end


function ret = is_zero(vect)
	ret = 1;
	for bla = 1:length(vect)
		if vect(bla) ~= 0
			ret = 0;
		end
	end
end
