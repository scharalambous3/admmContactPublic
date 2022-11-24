%%
clear
close all
clc
yalmip('clear')
%%
params = getBipedParams();
N = params.N; finalTime = params.finalTime; xDes = params.xDes;


%%
% initZ = zeros(params.orthDim, N - 1);
% stanceIndices = 0.5/params.dt;
% reps = (N-1)/stanceIndices;
% trotVec = repelem([1,0], stanceIndices);
% trotDefault = repmat(trotVec,1, reps/2);
% initZ(1, :) = trotDefault(1:N-1);
% initZ(2, :) = 1 - initZ(1, :);
% initZ(:, 1:5) = 1;
% initZ(:, end-5:end) = 1;
% 
% X0 = [[linspace(0, xDes(1), N); 0.5 * ones(1,N)];...
%     [(xDes(1)/params.horizon) * ones(1, N); zeros(1,N)];...
%     [linspace(0, xDes(1), N) - 0.25; linspace(0, xDes(1), N) + 0.25]];
% F0 = [initZ(1, :);...
%      9.8 * params.m * initZ(1, :);...
%      initZ(2, :);...
%      9.8 * params.m * initZ(2, :);...
%      1-initZ];
%%
x_k = params.X0;

dim = params.dim;
Delta0=zeros(params.dim , N-1);
%Delta0 = [X0(:,1:end-1); F0];
P0 = zeros(params.dim, (N-1));
G0 = (params.rho/2) * eye(params.dim, dim);

[Z,Delta, intVars, orthViol, objValueVec, zdeltaViolVec] = solveADMM(Delta0, P0, G0, x_k, params);

%%
X=Z(1:params.nx,:);
U=Z(1:params.nu,:);
%%
tRange=[0:params.dt:params.horizon-params.dt];
figure(1)
plot(tRange(1:end-1), X);
legend("x1", "x2", "x3", "x4", "x5", "x6")

figure(2)
title('Performance metrics')
subplot(2,2,1)
plot(1:size(orthViol, 2), orthViol);
xlabel('Iterations')
ylabel('orthViol')
ylim([0, inf])

subplot(2,2,2)
plot(1:size(objValueVec, 2), objValueVec);
xlabel('Iterations')
ylabel('Objective value')

subplot(2,2,3)
plot(1:size(zdeltaViolVec, 2), zdeltaViolVec);
xlabel('Iterations')
ylabel('z delta violation')
ylim([0, inf])

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
            drawFloatingBase(hax, c(:,k),f(:,k), xDes(1), [rci(1,k);0.0], [rci(2,k);0.0]);
        end
    
    end
    
    
end

function rotMatrix = rotz(theta)
rotMatrix = [cos(theta), - sin(theta), 0;...
             sin(theta), cos(theta), 0;...
             0, 0, 1];

end