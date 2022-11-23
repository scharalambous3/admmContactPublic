function [params] = getBipedParams()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
params.dt = 0.05;
params.N = 61;
params.horizon = params.N * params.dt;

params.finalTime = 3.0;
params.simSteps = params.finalTime/params.dt;

params.nx = 6;
params.nu = 6;
params.dim = params.nx + params.nu;
params.animation = false;

params.X0 = [0; 0.5; 0; 0; -0.25; 0.25];
params.xDes= [params.finalTime * 0.5; 0.5; 0; 0; 0; 0];

params.g = [0; -9.8];

params.m = 50.0;
params.I = (params.m/12) * (0.75^2 + 0.5^2);

%Since dynamics are linear now, Ts doesnt matter. Once I have NL dynamics
%Ts can be smaller than dt for better integration
params.Ts = params.dt;
A =zeros(params.nx, params.nx);
A(1:2, 3:4) = eye(2); %cdot
B = zeros(params.nx, params.nu);
B(3:4, 1:2) = (1/params.m) * eye(2); %f1
B(3:4, 3:4) = (1/params.m) * eye(2); %f2
B(5:6, 5:6) = eye(2); %rdocit

d = zeros(params.nx,1);
d(3:4, 1) = params.g;

params.A = eye(params.nx) + params.Ts * A;
params.B = params.Ts*B;
params.d = params.Ts*d;

params.M=1000;

params.Q = diag([500000, 100000, 100, 100, 0, 0]);
params.R = diag([0.1, 0.01, 0.1, 0.01, 1, 1]);
%params.Qf = idare(params.A, params.B, params.Q, params.R,[],[]);
params.Qf = params.Q;

%ADMM
params.epsDyn = 1e-16;
params.rho = 0.2;
params.rhoScale = 2;
params.maxIters=10;


%rci box constraint Ax * X >= bx
%     X(1, i) - 0.25 - 0.1 <= X(5, i), X(5, i) <= X(1, i) - 0.25 + 0.1
%     X(1, i) + 0.25 - 0.1 <= X(6, i), X(6, i) <= X(1, i) + 0.25 + 0.1
params.Ax = [-1, 0, 0, 0, 1, 0;...
                1, 0, 0, 0,-1, 0;...
                1, 0, 0, 0, 0, -1;...
               -1, 0, 0, 0, 0, 1];
params.bx = [-0.35; 0.15; -0.35; 0.15];

params.mu=0.9;
%unilateral force constraint Au U >= bu.  
params.Au = zeros(6, params.nu);
params.Au(1:2,[2, 4]) = eye(2); % f2>=0, f4>=0
params.Au(3:6,1:4) = blkdiag([-1, params.mu; 1,params.mu],[-1, params.mu; 1,params.mu]); % friction cone
params.bu = zeros(6, 1);

params.orthDim = 2;

params.Aorth = zeros(params.orthDim, params.dim);
params.Aorth(:, [8, 10]) = eye(params.orthDim);
params.aorth = zeros(params.orthDim, 1);
params.Borth = zeros(params.orthDim, params.dim);
params.Borth(:, [11, 12]) = eye(params.orthDim);
params.borth = zeros(params.orthDim, 1);

%Constraint Adelta_delta * Delta + Adelta_int * int >= bdelta
params.Adelta_delta = [params.Borth; params.Aorth; -params.Aorth; params.Borth; - params.Borth];

params.Adelta_int = [zeros(params.orthDim);...  %%Borth * Z + borth >= 0
    params.M * eye(params.orthDim);... % - M*i <= Aorth * Z + aorth
    params.M * eye(params.orthDim);... % Aorth * Z  + aorth <= + M*i 
    -params.M * eye(params.orthDim);... % -M*(1-i) <= Borth * Z + borth
    -params.M * eye(params.orthDim)]; % Borth * Z + borth <= M*(1-i)

params.bdelta = [- params.borth; - params.aorth; params.aorth; - params.M * ones(params.orthDim, 1) - params.borth; - params.M * ones(params.orthDim, 1) + params.borth];

%Initialization
initZ = zeros(params.orthDim, N - 1);
stanceIndices = 0.5/params.dt;
reps = (N-1)/stanceIndices;
trotVec = repelem([1,0], stanceIndices);
trotDefault = repmat(trotVec,1, reps/2);
initZ(1, :) = trotDefault(1:N-1);
initZ(2, :) = 1 - initZ(1, :);
initZ(:, 1:5) = 1;
initZ(:, end-5:end) = 1;

X0 = [[linspace(0, xDes(1), N); 0.5 * ones(1,N)];...
    [(xDes(1)/params.horizon) * ones(1, N); zeros(1,N)];...
    [linspace(0, xDes(1), N) - 0.25; linspace(0, xDes(1), N) + 0.25]];
F0 = [initZ(1, :);...
     9.8 * params.m * initZ(1, :);...
     initZ(2, :);...
     9.8 * params.m * initZ(2, :);...
     1-initZ];
params.Delta0 = [X0(:,1:end-1); F0];
params.P0= zeros(params.dim, (params.N-1));
end