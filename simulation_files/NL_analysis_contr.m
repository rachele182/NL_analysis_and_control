%% CONTROLLABILITY ANALYSIS FOR NON LIN SYS

run init.m

% creation of system in state affine control form
% state vector
syms x dx theta dtheta y dy T phi real 
% Parameters
syms mu_x mu_y l M grav J dphi real

q = [x dx theta dtheta y dy phi]';

f = [   dx;
        M^-1*(-T*sin(theta+phi) - mu_x*dx);
        dtheta;
        J^-1*(-T*l*sin(phi));
        dy;
        M^-1*(T*cos(theta+phi) - mu_y*dy -M*grav);
        0]; 
    
% u = dphi;
u = [dphi];

g = [0 0 0 0 0 0 1]';


%% Substitute parameters values in delta
% f = subs(f,{mu_x mu_y l M grav J T},{mu_x_par mu_y_par l_par M_par g_par J_par T_par{);


%% filtration method for non-lin sys
% this method demonstrates the small-time local accessibility of NON-LIN sys; 
% delta_0 = span{g}; delta = span{f,g}; delta_1 =  delta_0 + [delta_0,delta];
% [delta_0,delta] = Lgf ....

% initialization
q0 = zeros(7,1); 
delta_lb{1} = g;
delta = {f,g};
delta_mat = [delta_lb{1}];
r{1} = rank(delta_mat);

% iterative steps filtration
k=1;
tol = n_sys;
delta_mat_prec = zeros(size(delta_mat));
while k<tol
    [delta_lb,delta_mat,r{k+1}] = filtration(delta_lb, delta_mat, delta,q');
    delta_mat_num = subs(delta_mat, {q(1) q(2) q(3) q(4) q(5) q(6) q(7)}, {zeros(1,7)});
    r_num = rank(delta_mat_num);
    delta_mat_prec_num = subs(delta_mat_prec, {q(1) q(2) q(3) q(4) q(5) q(6) q(7)}, {zeros(1,7)});
    r_num_prec = rank(delta_mat_prec_num);
    if r_num == r{k+1} && r_num_prec == r{k} %check if deltak deltak+1 not singular
        if r{k+1} == r{k}    %check stop condition
           k = tol;
        end 
    end
    k = k+1;
    delta_mat_prec = delta_mat;
end
delta_mat = simplify(delta_mat);


function [delta_lb,delta_mat,r] = filtration(delta_lb_in,delta_mat_in,delta_in,q)
    delta_mat = delta_mat_in;
    delta_lb = delta_lb_in;
    r = rank(delta_mat);
    for index1 = 1:size(delta_lb_in,2)
        for index2 = 1:size(delta_in,2)
            lie_tmp = lie_b(delta_lb_in{index1},delta_in{index2},q);
            delta_mat = [delta_mat, lie_tmp];
            r = rank(delta_mat);
            if is_zero(lie_tmp) == 0                     % check if new vector is null
                delta_lb{end+1} = lie_tmp;
            else
                delta_mat = delta_mat(:,1:end-1);
            end
        end
    end
end


function Lgf = lie_b(g,f,q)
    Lgf = jacobian(f,q)*g - jacobian(g,q)*f;
end

function ret = is_zero(vect)
	ret = 1;
	for bla = 1:length(vect)
		if vect(bla) ~= 0
			ret = 0;
		end
	end
end
    