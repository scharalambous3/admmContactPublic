clear;
clc;
load('projectionIssue.mat');
%load('projectionIssueModified.mat');
[Z_k, objValue] = solveConvexSubproblem(Delta_k, P_k, G_k, x_k, params);

%%
[Delta_k, int_k] = projectionSubproblem(Z_k, P_k, G_k, params);
