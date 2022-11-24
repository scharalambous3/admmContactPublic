%%
clear
close all
clc
yalmip('clear')
%%
params = getBipedParams();
x_k = params.X0;
dim = params.dim;
Delta0 = params.Delta0;
P0 = params.P0;
G0 = (params.rho/2) * eye(params.dim, dim);

[Z,Delta, intVars, orthViol, objValueVec, primalResidualVec, dualResidualVec] = solveADMM(Delta0, P0, G0, x_k, params);

%%
% X=Z(1:params.nx,:);
% U=Z(1:params.nu,:);

[X, U, cost] = getRollout(Z, x_k, params);
%%
tRange=[0:params.dt:params.horizon-params.dt];
figure(1)
plot(tRange, X);
legend("x1", "x2", "x3", "x4", "x5", "x6")

figure(2)
title('Performance metrics')
subplot(2,2,1)
plot(1:size(orthViol, 2), orthViol);
xlabel('Iterations')
ylabel('Orthogonality violation')
ylim([0, inf])

subplot(2,2,2)
plot(1:size(objValueVec, 2), objValueVec);
xlabel('Iterations')
ylabel('Objective value')

subplot(2,2,3)
plot(1:size(primalResidualVec, 2), primalResidualVec);
xlabel('Iterations')
ylabel('z delta violation')
ylim([0, inf])
title('Prima Residual')

subplot(2,2,4)
plot(1:size(dualResidualVec, 2), dualResidualVec);
xlabel('Iterations')
ylabel('Dual residual')
ylim([0, inf])
title('Dual Residual')

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