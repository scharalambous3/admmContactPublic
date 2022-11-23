function [xTraj, uTraj, cost] = getRollout(Z_k, x_k, params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
cost = 0;

nx = params.nx; nu = params.nu; Q = params.Q; Qf = params.Qf; R = params.R; N = params.N;
assert(nx+nu==size(Z_k, 1)); assert(size(Q, 1) == nx);

xc = x_k;
xTraj = [];
uTraj = [];
for i = 1:(N - 1)
    uc = Z_k(nx + 1:end, i);

    xTraj = [xTraj, xc];
    uTraj = [uTraj, uc];

    cost = cost +  (xc - params.xDes)' * Q * (xc - params.xDes) +  uc' * R * uc;

    %simulate forward
    xc = discreteDyn(xc, uc, params);
end
xTraj = [xTraj, xc];
cost = cost + (xc - params.xDes)' * Qf * (xc - params.xDes);
end