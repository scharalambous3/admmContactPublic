function [Delta_k, int_k] = projectionSubproblem(Z_k, P_k, params, ~)
%UNTITLED4 Summary of this function goes here
% Detailed explanation goes here
N = params.N;
Delta_k = [];
int_k = [];
ops = sdpsettings('solver','mosek','cachesolvers',1,'verbose',0);

G_k = params.projG_k;
%For cartpole
%G_k(end, end) = 0;


for i = 1:(N - 1)
    delta = sdpvar(params.dim, 1);
    intVar = binvar(params.orthDim, 1);
 
    %consensus cost
    obj = (Z_k(:, i) - delta + P_k(:, i))' * G_k * (Z_k(:, i) - delta + P_k(:, i));

    constr = [];

    %TODO: constraints on integer variables
    
    constr = [constr, params.Adelta_delta * delta + params.Adelta_int * intVar >= params.bdelta];

    %For biped
    constr = [constr, [params.Ax zeros(size(params.Ax, 1), params.nu)] * delta >= params.bx];
    constr = [constr, [zeros(size(params.Au, 1), params.nx) params.Au] * delta >= params.bu];

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

