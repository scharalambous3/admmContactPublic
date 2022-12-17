function [Delta_k, int_k] = projectionSubproblemGrouping(Z_k, P_k, params, prevInt_k)
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

for i = 1:(N - 1)/params.groupingN
    delta = sdpvar(params.dim, params.groupingN);
    intVar = binvar(params.orthDim, 1);

    obj = 0;
    constr = [];

    %For biped
    if ((i ==1) || (i == (N - 1)/params.groupingN))
        %constr = [constr, intVar(4) == 1, intVar(8) == 1];
        %Alt
        %constr = [constr, intVar(1) == 1, intVar(2) == 1];
    end
    for j = 1:params.groupingN
        ndx = (i - 1) * params.groupingN + j;

        if (~isempty(prevInt_k))
            %obj = obj + (intVar - prevInt_k(:, ndx))' * params.RInt * (intVar - prevInt_k(:, ndx));
        end

        scalingFactor1 = exp(abs(xc(5) - Z_k(5, ndx)));
        scalingFactor2 = exp(abs(xc(6) - Z_k(6, ndx)));
        G_k(7:8, 7:8) = params.projG_k(7:8, 7:8) *  diag([scalingFactor1, scalingFactor2]);
        G_k(17:18, 17:18) = params.projG_k(17:18, 17:18) *  diag([scalingFactor1, scalingFactor2]);

        %consensus cost
        obj = obj + (Z_k(:, ndx) - delta(:, j) + P_k(:, ndx))' * G_k * (Z_k(:, ndx) - delta(:, j) + P_k(:, ndx));
    
        constr = [constr, params.Adelta_delta * delta(:, j) + params.Adelta_int * intVar >= params.bdelta];
    
        %For biped
        constr = [constr, [params.Ax zeros(size(params.Ax, 1), params.nu)] * delta(:, j) >= params.bx];
        constr = [constr, [zeros(size(params.Au, 1), params.nx) params.Au] * delta(:, j) >= params.bu];
        %For walking
        %constr = [constr, intVar(4) + intVar(8) >= 1];
        %Alt
        %constr = [constr, intVar(1) + intVar(2) >= 1];


    end
    optimize(constr,obj,ops);
    Delta_k = [Delta_k, value(delta)];
    int_k = [int_k, repelem(value(intVar), 1, params.groupingN)];
    for k = 1:params.groupingN
        ndx = (i - 1) * params.groupingN + k;
        uc = value(delta(params.nx+1:end, k)) - P_k(params.nx+1:end, ndx);
        xc = discreteDyn(xc, uc, params);
    end
end
%For biped
Delta_k(1:params.nx, 1)= Z_k(1:params.nx, 1);

%For cartpole
%Delta_k(1:params.nx,1) = 0;
%Delta_k(7,isnan(Delta_k(7,:))) = 0; %since u is not part of objective nor of constraints it returns NaN in cartpole

assert(all(size(Delta_k) == [params.dim, N - 1]))
end

