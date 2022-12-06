function [params] = getBipedParams()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
params.dt = 0.005;
params.N = 101;
params.horizon = params.N * params.dt;
params.NLInitialization=0;
%params.convexSubproblemSettings = sdpsettings('solver','snopt','cachesolvers',1,'allownonconvex',1, 'usex0', params.NLInitialization);%, 'snopt.Iterations_limit', 500);%, 'osqp.time_limit', 0.01);
params.convexSubproblemSettings = sdpsettings('solver','mosek','cachesolvers',1,'allownonconvex',0);%, 'snopt.Iterations_limit', 500);%, 'osqp.time_limit', 0.01);
params.finalTime = 3.0;
params.simSteps = params.finalTime/params.dt;

params.nx = 12;
lambdaDim = 8;
uDim = 4;
params.nu = lambdaDim + uDim;
params.dim = params.nx + params.nu;
params.animation = false;
params.liveGraphs = true;

params.X0 = [0; 0.5; 0; 0; -0.25; 0; 0.25; 0; zeros(4, 1)];
params.xDes= [0 * params.finalTime * 0.5; 5; 0; 0;...
    0 * params.finalTime * 0.5 - 0.25; 5; 0 * params.finalTime * 0.5 - 0.25; 5; zeros(4, 1)];

%ADMM
params.epsDyn = 1e-16;
%params.rho = 0.2/1000;
params.rho = 1;
params.rhoScale = 1.2;
params.maxIters=30;
params.epsilon0 = 100;
%The state, fx1 and fx2 do not enter the orthogonality constraint
%In the projection subproblem I need nonzero entries in G to retun non-NaNs
params.projG_k = blkdiag(eye(params.nx), eye(params.nu));
%The elements of omega and delta corresponding to these are irrelevant by
%setting the respective weights in G to 0
params.G0 = (params.rho/2) * blkdiag(diag([zeros(1, 5), 1, 0, 1, zeros(1, 4)]),... %since in orthogonality only x6,x8,x9,x10,x11,x12 appear
                                     eye(params.nu));

params.g = [0; -9.8];
initialization = false;
m=10.0;
params.m = m;
params.I = (params.m/12) * (0.75^2 + 0.5^2);

params.Q = diag([50000, 50000, 50000, 50000, 0, 0, 0, 0, 100, 100, 100, 100]);
params.R = diag([zeros(1, 1), 1, 1, 1, zeros(1, 1), 1, 1, 1 , 0.001, 0.001, 0.001, 0.001]);
%params.Qf = idare(params.A, params.B, params.Q, params.R,[],[]);
params.Qf = params.Q;

%Since dynamics are linear now, Ts doesnt matter. Once I have NL dynamics
%Ts can be smaller than dt for better integration
params.Ts = params.dt;
h = params.Ts;
A =zeros(params.nx, params.nx);
A(1:2, 3:4) = eye(2); %cdot
A(5:8, 9:12) = eye(4); %rdot
B=[0  -(h^2)/m   (h^2)/m    0     0   -(h^2)/m   (h^2)/m     0         0         0         0     0;...
   0     0          0   (h^2)/m   0       0        0      (h^2)/m      0         0         0     0;...
   0   -h/m        h/m      0     0      -h/m     h/m        0         0         0         0     0;...
   0     0          0      h/m    0       0        0        h/m        0         0         0     0;...
   0     0          0       0     0       0        0         0        h^2        0         0     0;...
   0     0          0       0     0       0        0         0         0        h^2        0     0;...
   0     0          0       0     0       0        0         0         0         0        h^2    0;...
   0     0          0       0     0       0        0         0         0         0         0    h^2;...
   0     0          0       0     0       0        0         0         h         0         0     0;...
   0     0          0       0     0       0        0         0         0         h         0     0;...
   0     0          0       0     0       0        0         0         0         0         h     0;...
   0     0          0       0     0       0        0         0         0         0         0     h];

d = zeros(params.nx,1);
d(1:2, 1) = (h^2) * params.g;
d(3:4, 1) = h * params.g;

params.A = eye(params.nx) + params.Ts * A;
params.B = B;
params.d = d;

params.M=1000;

%rci box constraint Ax * X >= bx
%     cx - 0.25 - 0.1 <= rc1x, rc1x <= cx - 0.25 + 0.1
%     cx + 0.25 - 0.1 <= rc2x, rc2x <= cx + 0.25 + 0.1
%     cy >= cyDes - 0.15
%     cy <= cyDes + 0.15 --> - (cyDes + 0.15) <= - X(2, i)
% params.Ax = [-1, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0;...
%              1, 0, 0, 0,-1,  0,  0, 0, 0, 0, 0, 0;...
%              1, 0, 0, 0, 0,  0, -1, 0, 0, 0, 0, 0;...
%             -1, 0, 0, 0, 0,  0,  1, 0, 0, 0, 0, 0;...
%              0, 1, 0, 0, 0,  0,  0, 0, 0, 0, 0, 0;...
%              0, -1, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0];
% params.bx = [-0.35; 0.15; -0.35; 0.15; params.xDes(2) - 0.15; - (params.xDes(2) + 0.15)];
%     cy - 0.5 - 0.1 <= rc1y, rc1y <= cy - 0.5 + 0.1
%     cy - 0.5 - 0.1 <= rc2y, rc2y <= cy - 0.5 + 0.1
params.Ax = [-1, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0;...
             1, 0, 0, 0,-1,  0,  0, 0, 0, 0, 0, 0;...
             1, 0, 0, 0, 0,  0, -1, 0, 0, 0, 0, 0;...
            -1, 0, 0, 0, 0,  0,  1, 0, 0, 0, 0, 0;...
             0, -1, 0, 0, 0, 1,  0, 0, 0, 0, 0, 0;...
             0, 1, 0, 0,  0,-1,  0, 0, 0, 0, 0, 0;...
             0, -1, 0, 0, 0, 0,  0, 1, 0, 0, 0, 0;...
             0, 1, 0, 0,  0, 0,  0, -1, 0, 0, 0, 0;...
             0, 0, 0, 0,  0, 1,  0, 0, 0, 0, 0, 0;...
             0, 0, 0, 0,  0, 0,  0, 1, 0, 0, 0, 0];
params.bx = [-0.35; 0.15; -0.35; 0.15; -0.6; 0.4; -0.6; 0.4; 0; 0];

mu=0.9;
params.mu = mu;
%unilateral force constraint Au * U >= bu.  
params.Au = zeros(2, params.nu);
params.Au(1:2,[4, 8]) = eye(2); % f2>=0, f4>=0
params.bu = zeros(2, 1);

params.orthDim = lambdaDim;

params.Aorth = [zeros(params.orthDim, params.nx), eye(lambdaDim), zeros(params.orthDim, uDim)];
params.aorth = zeros(params.orthDim, 1);
E=[0     0     0     0     0     0     0     0     0     0     0     0;...
   0     0     0     0     0     0     0     0     1     0     0     0;...
   0     0     0     0     0     0     0     0    -1     0     0     0;...
   0     0     0     0     0     1     0     0     0     h     0     0;...
   0     0     0     0     0     0     0     0     0     0     0     0;...
   0     0     0     0     0     0     0     0     0     0     1     0;...
   0     0     0     0     0     0     0     0     0     0    -1     0;...
   0     0     0     0     0     0     0     1     0     0     0     h];

FH=[0    -1    -1    mu     0     0     0     0     0     0     0     0;...
    1    0     0     0      0     0     0     0     h     0     0     0;...
    1    0     0     0      0     0     0     0    -h     0     0     0;...
    0    0     0     0      0     0     0     0     0     h^2   0     0;...
    0    0     0     0      0    -1    -1     mu     0     0     0     0;...
    0    0     0     0      1     0     0     0     0     0     h     0;...
    0    0     0     0      1     0     0     0     0     0    -h     0;...
    0    0     0     0      0     0     0     0     0     0     0     h^2];
params.Borth = [E, FH]; % q = Ex + F lambda + H u
%params.Borth(:, [11, 12]) = eye(params.orthDim);
params.borth = zeros(params.orthDim, 1);

%Constraint Adelta_delta * Delta + Adelta_int * int >= bdelta
params.Adelta_delta = [params.Aorth; params.Borth; params.Aorth; -params.Aorth; params.Borth; - params.Borth];

params.Adelta_int = [zeros(params.orthDim);...  %%Aorth * Z + aorth >= 0
    zeros(params.orthDim);...  %%Borth * Z + borth >= 0
    params.M * eye(params.orthDim);... % - M*i <= Aorth * Z + aorth
    params.M * eye(params.orthDim);... % Aorth * Z  + aorth <= + M*i 
    -params.M * eye(params.orthDim);... % -M*(1-i) <= Borth * Z + borth
    -params.M * eye(params.orthDim)]; % Borth * Z + borth <= M*(1-i)

params.bdelta = [- params.aorth; - params.borth; - params.aorth; params.aorth; - params.M * ones(params.orthDim, 1) - params.borth; - params.M * ones(params.orthDim, 1) + params.borth];

%Initialization
initZ = zeros(params.orthDim, params.N - 1);
stanceIndices = 0.25/params.dt;
reps = (params.N-1)/stanceIndices;
trotVec = repelem([1,0], stanceIndices);
trotDefault = repmat(trotVec,1, reps/2);
initZ(1, :) = trotDefault(1:params.N-1);
initZ(2, :) = 1 - initZ(1, :);
initZ(:, 1:5) = 1;
initZ(:, end-5:end) = 1;

params.initZ = initZ;

params.X0Init = [[linspace(0, params.xDes(1), params.N); 0.5 * ones(1,params.N)];...
    [(params.xDes(1)/params.horizon) * ones(1, params.N); zeros(1,params.N)];...
    [linspace(0, params.xDes(1), params.N) - 0.25; linspace(0, params.xDes(1), params.N) + 0.25]];
params.U0Init = [initZ(1, :);...
     9.8 * params.m * initZ(1, :);...
     initZ(2, :);...
     9.8 * params.m * initZ(2, :);...
     (1-initZ)];
if (initialization)
    params.Delta0 = [params.X0Init(:,1:end-1); params.U0Init];
else
    params.Delta0=zeros(params.dim , params.N-1);
end
params.P0= zeros(params.dim, (params.N-1));
end