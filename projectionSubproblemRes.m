function [Delta_k, int_k] = projectionSubproblemAlt(Z_k, P_k, params, ~)
%UNTITLED4 Summary of this function goes here
% Detailed explanation goes here
N = params.N;
Delta_k = [];
int_k = [];
ops = sdpsettings('solver','mosek','cachesolvers',1,'verbose',0);

G_k = params.projG_k;
%For cartpole
%G_k(end, end) = 0;

xc = params.X0;
res = zeros(params.dim, 1);


for i = 1:(N - 1)
    delta = sdpvar(params.dim, 1);
    intVar = binvar(params.orthDim, 1);
    constr = [];
    %For biped
    %if ((i ==1) || (i == (N - 1)))
    %    constr = [constr, intVar(1) == 1, intVar(2) == 1];
    %end
    %scalingFactor1 = (100 * (xc(5) - Z_k(5, i)))^4;
    %scalingFactor2 = (100 * (xc(6) - Z_k(6, i)))^4;
    %G_k(7:8, 7:8) = params.projG_k(7:8, 7:8) *  diag([scalingFactor1, scalingFactor2]);
    %G_k(17:18, 17:18) = params.projG_k(17:18, 17:18) *  diag([scalingFactor1, scalingFactor2]);

    %consensus cost
    obj = (Z_k(:, i) - delta + P_k(:, i) + res)' * G_k * (Z_k(:, i) - delta + P_k(:, i) + res);

    if (i > 1)
        %obj = obj + (intVar - int_k(:, end))' * params.RInt * (intVar - int_k(:, end));
    end

    %For walking
    %constr = [constr, intVar(1) + intVar(2) >= 1];
    
    constr = [constr, params.Adelta_delta * delta + params.Adelta_int * intVar >= params.bdelta];

    %For biped
    constr = [constr, [params.Ax zeros(size(params.Ax, 1), params.nu)] * delta >= params.bx];
    constr = [constr, [zeros(size(params.Au, 1), params.nx) params.Au] * delta >= params.bu];

    optimize(constr,obj,ops);
    Delta_k = [Delta_k, value(delta)];
    int_k = [int_k, value(intVar)];

    uc = value(delta(params.nx+1:end));% - P_k(params.nx+1:end, i);
    xc = discreteDyn(xc, uc, params);
    
    res = res + Z_k(:, i) - value(delta);
end
%For biped
%Delta_k(1:params.nx,:)= Z_k(1:params.nx,:);

%For cartpole
%Delta_k(1:params.nx,1) = 0;
%Delta_k(7,isnan(Delta_k(7,:))) = 0; %since u is not part of objective nor of constraints it returns NaN in cartpole

assert(all(size(Delta_k) == [params.dim, N - 1]))
end

