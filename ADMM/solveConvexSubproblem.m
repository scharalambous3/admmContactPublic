function [Z_k, objValue] = solveConvexSubproblem(Delta_k, P_k, G_k, x_k, params, epsilon)
%solveConvexSubproblem Solve convex subproblem (QP since linear
%complementarity constraints are not included here)

Q=params.Q; Qf=params.Qf; R=params.R; N=params.N;

ops = params.convexSubproblemSettings;
X = sdpvar(params.nx,N);
U = sdpvar(params.nu,N - 1);

if (params.NLInitialization)
    assign(X, params.X0Init);
    assign(U, params.U0Init);
end

obj =0;
constr = [X(:,1) == x_k];

for i = 1:(N - 1)
    %quadratic state input cost
    obj = obj + (X(:,i) - params.xDes(:, i))' * Q * (X(:,i) - params.xDes(:, i)) + U(:, i)' * R * U(:, i);

    %consensus cost
    obj = obj + ([X(:,i); U(:, i)] - Delta_k(:, i) + P_k(:, i))' * G_k * ([X(:,i); U(:, i)] - Delta_k(:, i) + P_k(:, i));

    %dynamics
    constr = [constr, X(:, i + 1) == discreteDyn(X(:, i), U(:, i), params)];
    
    %Linear ineq constraints on X
    if ~isempty(params.Ax)
    constr = [constr, params.Ax * X(:,i) >= params.bx];
    end
    %Linear ineq constraints on U
    if ~isempty(params.Au)
    constr = [constr, params.Au * U(:,i) >= params.bu];
    end

    %relaxed orthogonality constraint
    %constr = [constr, (params.Aorth(:, 1:params.nx) * X(:,i) + params.Aorth(:, params.nx+1:end) * U(:,i) + params.aorth)' *...
    %          (params.Borth(:, 1:params.nx) * X(:,i) + params.Borth(:, params.nx+1:end) * U(:,i) + params.borth) <=epsilon];

end
% terminal cost
obj = obj + (X(:, N) - params.xDes(:, end))' * Qf * (X(:, N) - params.xDes(:, end));

%Terminal state constraint
if ~isempty(params.AxTerminal)
    constr = [constr, params.AxTerminal * X(:,N) == params.bxTerminal];
end
opt = optimize(constr,obj,ops);

Z_k = [value(X(:,1:N-1)); value(U)];
objValue = value(obj);
end