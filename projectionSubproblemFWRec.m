function [Delta_k, int_k] = projectionSubproblemFWRec(Z_k, P_k, params, prevInt_k)
%UNTITLED4 Summary of this function goes here
% Detailed explanation goes here
N = params.N;
Delta_k = [];
int_k = [];
ops = sdpsettings('solver','mosek','cachesolvers',1,'verbose',0);

%TODO Temporary hack since Delta_k P_k and Z_k lack xN
Z_k = [Z_k, [params.xDes(:, end); zeros(params.nu, 1)]];
P_k = [P_k, [params.xDes(:, end); zeros(params.nu, 1)]];

T = 1000 * diag([1000, 50000, 1000, 50000, 50000, 50000, 1000, 1000,]); %c, cdot, rx, rdotx

G_k = params.projG_k;
xc = params.X0;

for i = 1:(N - 1)
    delta = sdpvar(params.dim, 1);
    intVar = binvar(params.orthDim, 1);
 
    %consensus cost
    %Just so that NaNs do not appear
    obj = 0.000001 * (Z_k(:, i) - delta + P_k(:, i))' * G_k * (Z_k(:, i) - delta + P_k(:, i));
    obj = obj + (Z_k(1:params.nx, i+1) - discreteDyn(xc, delta(params.nx+1:end), params) + P_k(1:params.nx, i+1))' *...
        T * (Z_k(1:params.nx, i+1) - discreteDyn(xc, delta(params.nx+1:end), params) + P_k(1:params.nx, i+1));

    if (~isempty(prevInt_k))
       obj = obj + (intVar - prevInt_k(:, i))' * params.RInt * (intVar - prevInt_k(:, i));
    end

    constr = [delta(1:params.nx) == xc];

    %TODO: constraints on integer variables
    
    constr = [constr, params.Adelta_delta * delta + params.Adelta_int * intVar >= params.bdelta];
    tempConstr = constr;
    feasible = false;
    relaxEps = 1e-6;
    while ~feasible
        if ~isempty(params.Ax)
            constr = [tempConstr, [params.Ax zeros(size(params.Ax, 1), params.nu)] * delta >= params.bx - relaxEps];
            constr = [tempConstr, params.Ax * discreteDyn(xc, delta(params.nx+1:end), params) >= params.bx - relaxEps];
        end
        if ~isempty(params.Au)
            constr = [tempConstr, [zeros(size(params.Au, 1), params.nx) params.Au] * delta >= params.bu - relaxEps];
        end
        relaxEps = relaxEps * 10;
        opt = optimize(constr,obj,ops);
        feasible = ~opt.problem;
    end


    Delta_k = [Delta_k, value(delta)];
    int_k = [int_k, value(intVar)];

    uc = value(delta(params.nx+1:end));
    xc = discreteDyn(xc, uc, params);
end
%For biped
%Delta_k(1:params.nx,:)= Z_k(1:params.nx,:);

%For cartpole
%Delta_k(1:params.nx,1) = 0;
%Delta_k(7,isnan(Delta_k(7,:))) = 0; %since u is not part of objective nor of constraints it returns NaN in cartpole


assert(all(size(Delta_k) == [params.dim, N - 1]))
end

