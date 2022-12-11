%%
clear
close all
clc
yalmip('clear')
%%
params = getBipedParamsImplicit();
x_k = params.X0;


Q=params.Q; Qf=params.Qf; R=params.R; N=params.N;

ops = sdpsettings('solver','mosek','cachesolvers',1,'allownonconvex',0,'verbose',1);
X = sdpvar(params.nx,N);
U = sdpvar(params.nu,N - 1);
intVar = binvar(params.orthDim, N - 1);

%assign(intVar, params.initZ)

obj =0;
constr = [X(:,1) == x_k];
Q = Q;
R = R;
Qf = Qf;
for i = 1:(N - 1)
    %quadratic state input cost
    obj = obj + (X(:,i) - params.xDes)' * Q * (X(:,i) - params.xDes) + U(:, i)' * R * U(:, i);

    %dynamics
    constr = [constr, X(:, i + 1) == discreteDyn(X(:, i), U(:, i), params)];
    constr = [constr, params.Adelta_delta * [X(:,i); U(:, i)] + params.Adelta_int * intVar(:, i) >= params.bdelta];

    %Linear ineq constraints on X
    if ~isempty(params.Ax)
    constr = [constr, params.Ax * X(:,i) >= params.bx];
    end
    %Linear ineq constraints on U
    if ~isempty(params.Au)
    constr = [constr, params.Au * U(:,i) >= params.bu];
    end

end
% terminal cost
obj = obj + (X(:, N) - params.xDes)' * Qf * (X(:, N) - params.xDes);

%Terminal state constraint
if ~isempty(params.AxTerminal)
    constr = [constr, params.AxTerminal * X(:,N) == params.bxTerminal];
end

opt = optimize(constr,obj,ops);

objValue = value(obj);

%%
xTraj = value(X);
uTraj = value(U);
Z_k = [value(X(:,1:N-1)); value(U)];
[xTrajRollout, uTrajRollout, cost] = getRollout(Z_k, x_k, params);
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

% subplot(3,3,3)
% plot(1:size(dynViol, 2), vecnorm(dynViol(1:3:end,:),2, 1));
% xlabel('Iterations')  
% ylabel('cdotx Violation')
% ylim([0, inf])
% 
% subplot(3,3,4)
% plot(1:size(dynViol, 2), vecnorm(dynViol(2:3:end,:),2, 1));
% xlabel('Iterations')
% ylabel('cdoty Violation')
% ylim([0, inf])
% 
% subplot(3,3,5)
% plot(1:size(dynViol, 2), vecnorm(dynViol(3:3:end,:),2, 1));
% xlabel('Iterations')
% ylabel('l Violation')
% ylim([0, inf])
% 
% subplot(3,3,6)
% plot(1:size(dynViol, 2), vecnorm(dynViol(4:3:end,:),2, 1));
% xlabel('Iterations')
% ylabel('rci1 Violation')
% ylim([0, inf])
% 
% subplot(3,3,7)
% plot(1:size(dynViol, 2), vecnorm(dynViol(5:3:end,:),2, 1));
% xlabel('Iterations')
% ylabel('rci2 Violation')
% ylim([0, inf])

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