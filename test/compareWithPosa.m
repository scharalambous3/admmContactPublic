clear;
close all;
%Compare projection function
clc;
params = getTestParams();
load('deltaProj.mat');
deltaExpected = reshape(deltaProjOutput(1:70), 7, 10);
G_k = diag([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.0]);
P_k = zeros(7, 10);

load('z.mat');
Z_kPosa = reshape(z(1:70), 7, 10);
[Delta_k, ~] = projectionSubproblem(Z_kPosa,P_k,G_k, params);


assert(norm(Delta_k-deltaExpected) < 1e-3);
%%
Delta0=zeros(params.dim , params.N-1);
X0 =[0.3; 0; 0.3; 0];

[Z_k, ~] = solveConvexSubproblem(Delta0, P_k, G_k, X0, params);
assert(norm(Z_k-Z_kPosa) < 0.01);

%%
load('zFinalPosa.mat');
load('zFinalMine.mat');
zFinalPosa = reshape(zFinalPosa(1:70), 7, 10);
assert(norm(zFinalPosa-zFinalMine) < 0.5);
fprintf('Norm difference between final Z %d \n', norm(zFinalPosa-zFinalMine)) 
