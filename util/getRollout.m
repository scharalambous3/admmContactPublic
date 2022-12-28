function [xTraj, uTraj, cost] = getRollout(Z_k, x_k, params)
%getRollout Returns state and input trajectories from rollout of solution
%from convex subproblem

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

    cost = cost +  (xc - params.xDes(:, i))' * Q * (xc - params.xDes(:, i)) +  uc' * R * uc;

    %simulate forward
    xc = discreteDyn(xc, uc, params);
end
xTraj = [xTraj, xc];
cost = cost + (xc - params.xDes(:, end))' * Qf * (xc - params.xDes(:, end));
end