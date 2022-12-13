%%
clear
close all
clc
yalmip('clear')
%%
params = getBipedParamsImplicit();
x_k = params.X0;
dim = params.dim;
Delta0 = params.Delta0;
P0 = params.P0;
G0 = params.G0;

[Z,Delta, intVars, orthViol, objValueVec, primalResidualArr, dualResidualArr] = solveADMM(Delta0, P0, G0, x_k, params);

%%
% X=Z(1:params.nx,:);
% U=Z(1:params.nu,:);

[X, U, cost] = getRollout(Z, x_k, params);
%%
tRange=[0:params.dt:params.horizon];
figure(1)
plot(tRange, X);
%legend("x1", "x2", "x3", "x4", "x5", "x6")

figure(2)
hax = gca;
plotPerf(hax, orthViol,objValueVec, primalResidualArr, dualResidualArr, intVars(params.separationIndices, :), X)
plotStates(intVars(params.separationIndices, :), Z, params);
plotProj(intVars(params.separationIndices, :),  Z, Delta, params);
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