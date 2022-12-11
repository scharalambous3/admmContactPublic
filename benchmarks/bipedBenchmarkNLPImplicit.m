%%
clear
close all
clc
yalmip('clear')
%%
params = getBipedParamsImplicit();

%For jumping
params.Ax(1:4,:) = [];
params.bx(1:4) = [];

ops = sdpsettings('solver','snopt','cachesolvers',1,'verbose',1,'usex0', 1, 'snopt.Iterations_limit', 10000);
%%
[X1, U1] = solveNLP(ops, params, params.X0Init, params.U0Init, 100);
%%
[X2, U2] = solveNLP(ops, params, X1, U1, 50);
%%
[X3, U3] = solveNLP(ops, params, X2, U2, 10);
%%
[X4, U4] = solveNLP(ops, params, X3, U3, 1);
%%
[X5, U5] = solveNLP(ops, params, X4, U4, 0.1);
%%
%xTraj = value(X);
%uTraj = value(U);
Z_k = [value(X5(:,1:params.N-1)); value(U5)];
[xTrajRollout, uTrajRollout, cost] = getRollout(Z_k, params.X0, params);
%%
tRange=[0:params.dt:params.horizon-params.dt];
figure(1)
plot(tRange, xTraj);
legend("x1", "x2", "x3", "x4", "x5", "x6")

figure(2)
mat = 1-value(intVar);  % Your sample matrix
[r, c] = size(mat);                          % Get the matrix size
imagesc((1:c)+0.5, (1:r)+0.5, mat);          % Plot the image
colormap(gray);                              % Use a gray colormap
set(gca, 'XTick', 1:(c+1), 'YTick', 1:(r+1), ...  % Change some axes properties
         'XLim', [1 c+1], 'YLim', [1 r+1], ...
         'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');

%%
if (params.animation)
    hax = figure(3);
    
    %gif('A.gif','DelayTime',1/5)
    
    while (true)
        for k=1:length(tRange)-1
            drawFloatingBase(hax, X(1:2,k),U(1:4,k), params.xDes(1), [X(5,k);0.0], [X(6,k);0.0]);
        end
    
    end
    
    
end

function rotMatrix = rotz(theta)
rotMatrix = [cos(theta), - sin(theta), 0;...
             sin(theta), cos(theta), 0;...
             0, 0, 1];

end

%%
function [X, U] = solveNLP(ops, params, X0Init, U0Init, epsilon)
Q=params.Q; Qf=params.Qf; R=params.R; N=params.N;
x_k = params.X0;

X = sdpvar(params.nx,N);
U = sdpvar(params.nu,N - 1);

assign(X, X0Init);
assign(U, U0Init);

obj =0;
constr = [X(:,1) == x_k];

for i = 1:(N - 1)
    %quadratic state input cost
    obj = obj + (X(:,i) - params.xDes(:, i))' * Q * (X(:,i) - params.xDes(:, i)) + U(:, i)' * R * U(:, i);

    %dynamics
    constr = [constr, X(:, i + 1) == discreteDyn(X(:, i), U(:, i), params)];
    constr = [constr, params.Aorth * [X(:,i); U(:, i)] + params.aorth >= 0];
    constr = [constr, params.Borth * [X(:,i); U(:, i)] + params.borth >= 0];

    constr = [constr,...
        (params.Aorth * [X(:,i); U(:, i)] + params.aorth)' * (params.Borth * [X(:,i); U(:, i)] + params.borth) <= epsilon];


    %Linear ineq constraints on X
    constr = [constr, params.Ax * X(:,i) >= params.bx];
    
    %Linear ineq constraints on U
    if ~isempty(params.Au)
    constr = [constr, params.Au * U(:,i) >= params.bu];
    end

end
% terminal cost
obj = obj + (X(:, N) - params.xDes(:, end))' * Qf * (X(:, N) - params.xDes(:, end));
if ~isempty(params.AxTerminal)
    constr = [constr, params.AxTerminal * X(:,N) == params.bxTerminal];
end

opt = optimize(constr,obj,ops);

objValue = value(obj);
X=value(X);
U=value(U);
end