function [Delta_k, int_k] = projectionSubproblem(Z_k, P_k, params)
%UNTITLED4 Summary of this function goes here
% Detailed explanation goes here
N = params.N;
Delta_k = [];
int_k = [];
ops = sdpsettings('solver','mosek','cachesolvers',1,'verbose',0);

%delta = [1cx, 2cy, 3cdotx, 4cdoty, 5rc1, 6rc2, 7fx1, 8fy1, 9fx2, 10fy2, 11rdotc1, 12rdotc2]
%delta = [delta_c∈R2;delta_cdot∈R2;delta_rci∈R2;delta_f∈R4;delta_rdotci∈R2]

G_k = params.projG_k;
%For cartpole
%G_k(end, end) = 0;


parfor (i = 1:(N - 1), 8)
    delta = sdpvar(params.dim, 1);
    intVar = binvar(params.orthDim, 1);
 
    %consensus cost
    obj = (Z_k(:, i) - delta + P_k(:, i))' * G_k * (Z_k(:, i) - delta + P_k(:, i));

    constr = [];
    % orthogonality constraints:  f>= 0, rdotci orthogonal to f

    %TODO: constraints on integer variables
    %constr = [constr, intVar(1) + intVar(2) >= 1];
    
    % zero force during swing
    constr = [constr, params.Adelta_delta * delta + params.Adelta_int * intVar >= params.bdelta];

    optimize(constr,obj,ops);
    Delta_k = [Delta_k, value(delta)];
    int_k = [int_k, value(intVar)];
end
%For biped
%Delta_k(1:params.nx,:)= Z_k(1:params.nx,:);

%For cartpole
%Delta_k(1:params.nx,1) = 0;
%Delta_k(7,isnan(Delta_k(7,:))) = 0; %since u is not part of objective nor of constraints it returns NaN in cartpole

assert(all(size(Delta_k) == [params.dim, N - 1]))
end

