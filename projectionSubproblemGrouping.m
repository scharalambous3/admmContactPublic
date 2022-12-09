function [Delta_k, int_k] = projectionSubproblemGrouping(Z_k, P_k, params)
%UNTITLED4 Summary of this function goes here
% Detailed explanation goes here
N = params.N;
Delta_k = [];
int_k = [];
ops = sdpsettings('solver','mosek','cachesolvers',1,'verbose',0);

G_k = params.projG_k;
%For cartpole
%G_k(end, end) = 0;


for i = 1:(N - 1)/params.groupingN
    delta = sdpvar(params.dim, params.groupingN);
    intVar = binvar(params.orthDim, 1);

    obj = 0;
    constr = [];
    for j = 1:params.groupingN
        ndx = (i - 1) * params.groupingN + j;
        %consensus cost
        obj = obj + (Z_k(:, ndx) - delta(:, j) + P_k(:, ndx))' * G_k * (Z_k(:, ndx) - delta(:, j) + P_k(:, ndx));
    
    
        %TODO: constraints on integer variables
        
        constr = [constr, params.Adelta_delta * delta(:, j) + params.Adelta_int * intVar >= params.bdelta];
    
        %For biped
        constr = [constr, [params.Ax zeros(size(params.Ax, 1), params.nu)] * delta(:, j) >= params.bx];
        constr = [constr, [zeros(size(params.Au, 1), params.nx) params.Au] * delta(:, j) >= params.bu];

    end
    optimize(constr,obj,ops);
    Delta_k = [Delta_k, value(delta)];
    int_k = [int_k, repelem(value(intVar), 1, params.groupingN)];
end
%For biped
%Delta_k(1:params.nx,:)= Z_k(1:params.nx,:);

%For cartpole
%Delta_k(1:params.nx,1) = 0;
%Delta_k(7,isnan(Delta_k(7,:))) = 0; %since u is not part of objective nor of constraints it returns NaN in cartpole

assert(all(size(Delta_k) == [params.dim, N - 1]))
end

