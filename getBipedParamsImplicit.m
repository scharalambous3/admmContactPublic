function [params] = getBipedParamsImplicit()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
params.dt = 0.05;
params.N = 51;
params.groupingN  = 0.25/params.dt;
params.horizon = (params.N-1) * params.dt;
params.NLInitialization=0;
%params.convexSubproblemSettings = sdpsettings('solver','snopt','cachesolvers',1,'allownonconvex',1, 'usex0', params.NLInitialization);%, 'snopt.Iterations_limit', 500);%, 'osqp.time_limit', 0.01);
params.convexSubproblemSettings = sdpsettings('solver','mosek','cachesolvers',1,'allownonconvex',0);%, 'snopt.Iterations_limit', 500);%, 'osqp.time_limit', 0.01);
params.finalTime = 3.0;
params.simSteps = params.finalTime/params.dt;

%params.separationIndices = [4,8];
params.separationIndices = [2, 4];

params.nx = 12;
lambdaDim = 8;
uDim = 4;
params.nu = lambdaDim + uDim;
params.dim = params.nx + params.nu;
params.animation = false;
params.liveGraphs = true;

params.X0 = [0; 0.5; 0; 0; -0.25; 0; 0.25; 0; zeros(4, 1)];
%For jumping
% xDesMiddle = [0; 5; 0; 0;...
%     - 0.25; 4.5; 0.25; 4.5; zeros(4, 1)];
% xDesFinal = [0; 0.5; 0; 0; -0.25; 0; 0.25; 0; zeros(4, 1)];
% params.xDes = repmat(xDesFinal, 1, params.N);
% midIndex = floor(params.N/2);
% params.xDes(2, :) = [linspace(params.X0(2), xDesMiddle(2), midIndex), linspace(xDesMiddle(2), xDesFinal(2), params.N - midIndex) ];
% params.xDes(6, :) = [linspace(params.X0(6), xDesMiddle(6), midIndex), linspace(xDesMiddle(6), xDesFinal(6), params.N - midIndex) ];
% params.xDes(8, :) = [linspace(params.X0(8), xDesMiddle(8), midIndex), linspace(xDesMiddle(8), xDesFinal(8), params.N - midIndex) ];

xDesFinal= [params.finalTime * 0.5; 0.5; 0; 0;...
   0; 0; 0; 0; zeros(4, 1)];
params.xDes = repmat(xDesFinal, 1, params.N);
% Sinusoidal cy reference
tRange=[0:params.dt:params.horizon];
%params.xDes(2,:) = xDesFinal(2) + 0.1 * sin(tRange/params.horizon * 4 * pi);

%ADMM
params.epsDyn = 1e-16;
%params.rho = 0.2/1000;
params.rho = 1;
params.rhoScale = 1.5;
params.maxIters= 20;
params.epsilon0 = 100;
%The state, fx1 and fx2 do not enter the orthogonality constraint
%In the projection subproblem I need nonzero entries in G to retun non-NaNs
%Note in orthogonality only x6,x8,x9,x10,x11,x12 appear                                    
params.projG_k = blkdiag(diag([ones(1, 4), 1000 * ones(1, 4), 100 * ones(1, 4), ones(1, 8), ones(1,4)])); 
%The elements of omega and delta corresponding to these are irrelevant by
%setting the respective weights in G to 0

params.G0 = (params.rho/2) * blkdiag(diag([ones(1, 4), 100 * ones(1, 4), 10 * ones(1, 4), ones(1, 8), ones(1,4)]));


params.g = [0; -9.8];
initialization = false;
m=10.0;
params.m = m;
params.I = (params.m/12) * (0.75^2 + 0.5^2);

params.Q = diag([50000, 50000, 1000, 1000, 0, 0, 0, 0, 1000, 1000, 1000, 1000]);
%params.R = diag([zeros(1, 1), 1, 1, 0.01, zeros(1, 1), 1, 1, 0.01 , 0.05, 0.05, 0.05, 0.05]);
params.R =  diag([1, 1, 1, 0.1, 1, 1, 1, 0.1, 5, 5, 5, 5]);
%params.Qf = idare(params.A, params.B, params.Q, params.R,[],[]);
params.Qf = diag([50000, 50000, 1000, 1000, 0, 0, 0, 0, 1000, 1000, 1000, 1000]);
params.RInt = diag([0,0,0,1,0,0,0,1]);

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
%     cy - 0.5 - 0.1 <= rc2y, rc2y <= cy - 0.5 + 0.1 --> cy - 0.6 <= rcy <= cy - 0.4
params.Ax = [-1, 0, 0, 0, 1, 0,  0, 0, 0, 0, 0, 0;...%box constraints for rcix
             1, 0, 0, 0,-1,  0,  0, 0, 0, 0, 0, 0;...
             1, 0, 0, 0, 0,  0, -1, 0, 0, 0, 0, 0;...
            -1, 0, 0, 0, 0,  0,  1, 0, 0, 0, 0, 0;...
             0, -1, 0, 0, 0, 1,  0, 0, 0, 0, 0, 0;...%box constraints for rciy
             0, 1, 0, 0,  0,-1,  0, 0, 0, 0, 0, 0;...
             0, -1, 0, 0, 0, 0,  0, 1, 0, 0, 0, 0;...
             0, 1, 0, 0,  0, 0,  0, -1, 0, 0, 0, 0;...
             0, 0, 0, 0,  0, 1,  0, 0, 0, 0, 0, 0;...%rci >= 0
             0, 0, 0, 0,  0, 0,  0, 1, 0, 0, 0, 0];
params.bx = [-0.35; 0.15; -0.35; 0.15; -0.6; 0.4; -0.6; 0.4; 0; 0];

mu=0.7;
params.mu = mu;
%unilateral force constraint Au * U >= bu.  
params.Au = [];
params.bu = [];
unilateralMat = zeros(6, params.nu);
unilateralMat(:,[2, 3, 4, 6, 7, 8]) = eye(6); % f1N>=0, f2N>=0, and fT+ and fT- >= 0
params.Au = [params.Au; unilateralMat];
params.bu = [params.bu; zeros(6, 1)];

frConesMat = zeros(2, params.nu);
frConesMat(1,[2,3,4]) = [-1, -1, mu]; %friction cone for contact 1
frConesMat(2,[6,7,8]) = [-1, -1, mu]; %friction cone for contact 2
params.Au = [params.Au; frConesMat];
params.bu = [params.bu; zeros(2, 1)];

forceLimitsMat = zeros(2, params.nu);
forceLimitsMat(:,[4, 8]) = -eye(2); % force limits
params.Au = [params.Au; forceLimitsMat];
params.bu = [params.bu;  -75 * ones(2, 1)];

% rddotLimitsMat = zeros(8, params.nu);
% rddotLimitsMat(1:4,[9:12]) = eye(4); % rddot limits
% rddotLimitsMat(5:8,[9:12]) = -eye(4); % rddot limits
% params.Au = [params.Au; rddotLimitsMat];
% params.bu = [params.bu;  -5 * ones(8, 1)];


params.AxTerminal = zeros(6,params.nx);
params.AxTerminal(:, [6,8,9,10,11,12]) = eye(6);
params.bxTerminal = zeros(6, 1);

%params.orthDim = lambdaDim;
params.orthDim = 4;
%Aorth is lambda
% params.Aorth = [zeros(params.orthDim, params.nx), eye(lambdaDim), zeros(params.orthDim, uDim)];
params.Aorth = zeros(params.orthDim, params.dim);
params.Aorth([1,3],[params.nx + 4, params.nx + 8]) = eye(2);
params.Aorth([2,4],[params.nx + 4, params.nx + 8]) = eye(2);
params.aorth = zeros(params.orthDim, 1);
%Borth is q
% E=[0     0     0     0     0     0     0     0     0     0     0     0;...
%    0     0     0     0     0     0     0     0     1     0     0     0;...
%    0     0     0     0     0     0     0     0    -1     0     0     0;...
%    0     0     0     0     0     1     0     0     0     h     0     0;...
%    0     0     0     0     0     0     0     0     0     0     0     0;...
%    0     0     0     0     0     0     0     0     0     0     1     0;...
%    0     0     0     0     0     0     0     0     0     0    -1     0;...
%    0     0     0     0     0     0     0     1     0     0     0     h];
% 
% FH=[0    -1    -1    mu     0     0     0     0     0     0     0     0;...
%     1    0     0     0      0     0     0     0     h     0     0     0;...
%     1    0     0     0      0     0     0     0    -h     0     0     0;...
%     0    0     0     0      0     0     0     0     0     h^2   0     0;...
%     0    0     0     0      0    -1    -1     mu     0     0     0     0;...
%     0    0     0     0      1     0     0     0     0     0     h     0;...
%     0    0     0     0      1     0     0     0     0     0    -h     0;...
%     0    0     0     0      0     0     0     0     0     0     0     h^2];

E=[0     0     0     0     0     0     0     0     1     0     0     0;... %rdot1x
   0     0     0     0     0     1     0     0     0     h     0     0;...%r1y + h rdot1y
   0     0     0     0     0     0     0     0     0     0     1     0;...%rdot2x
   0     0     0     0     0     0     0     1     0     0     0     h];%r2y + h rdot2y

FH=[0    0     0     0      0     0     0     0     h     0     0     0;...%h rddot1x
    0    0     0     0      0     0     0     0     0     h^2   0     0;...%h^2 rddot1y
    0    0     0     0      0     0     0     0     0     0     h     0;...%h rddot2x
    0    0     0     0      0     0     0     0     0     0     0     h^2];%h^2 rddot2y

params.Borth = [E, FH]; % q = Ex + F lambda + H u
%params.Borth(:, [11, 12]) = eye(params.orthDim);
params.borth = zeros(params.orthDim, 1);

%Constraint Adelta_delta * Delta + Adelta_int * int >= bdelta
%params.Adelta_delta = [params.Aorth; params.Borth; params.Aorth; -params.Aorth; params.Borth; - params.Borth];
params.Adelta_delta = [params.Aorth; params.Aorth; -params.Aorth; params.Borth; - params.Borth];

params.Adelta_int = [zeros(params.orthDim);...  %%Aorth * Z + aorth >= 0
    %zeros(params.orthDim);...  %%Borth * Z + borth >= 0
    params.M * eye(params.orthDim);... % - M*i <= Aorth * Z + aorth
    params.M * eye(params.orthDim);... % Aorth * Z  + aorth <= + M*i 
    -params.M * eye(params.orthDim);... % -M*(1-i) <= Borth * Z + borth
    -params.M * eye(params.orthDim)]; % Borth * Z + borth <= M*(1-i)

%params.bdelta = [- params.aorth; - params.borth; - params.aorth; params.aorth; - params.M * ones(params.orthDim, 1) - params.borth; - params.M * ones(params.orthDim, 1) + params.borth];
params.bdelta = [- params.aorth; - params.aorth; params.aorth; - params.M * ones(params.orthDim, 1) - params.borth; - params.M * ones(params.orthDim, 1) + params.borth];

%Initialization TODO: Update for implicit formulation
% initZ = zeros(params.orthDim, params.N - 1);
% stanceIndices = 0.25/params.dt;
% reps = (params.N-1)/stanceIndices;
% trotVec = repelem([1,0], stanceIndices);
% trotDefault = repmat(trotVec,1, reps/2);
% initZ(1, :) = trotDefault(1:params.N-1);
% initZ(2, :) = 1 - initZ(1, :);
% initZ(:, 1:5) = 1;
% initZ(:, end-5:end) = 1;
%params.initZ = initZ;

%For jumping
initZ = zeros(2, params.N - 1);
initZ(:, 1:5) = 1;
initZ(:, end-5:end) = 1;
params.X0Init = params.xDes;
params.U0Init = zeros(params.nu, params.N - 1);
params.U0Init(4, :) = 10 * 9.8 * params.m * initZ(1, :);
params.U0Init(8, :) = 10 * 9.8 * params.m * initZ(2, :);

% params.X0Init = [[linspace(0, params.xDes(1, end), params.N); 0.5 * ones(1,params.N)];...
%     [(params.xDes(1, end)/params.horizon) * ones(1, params.N); zeros(1,params.N)];...
%     [linspace(0, params.xDes(1, end), params.N) - 0.25; linspace(0, params.xDes(1, end), params.N) + 0.25]];
%TODO: Update for implicit formulation
% params.U0Init = [initZ(1, :);...
%      9.8 * params.m * initZ(1, :);...
%      initZ(2, :);...
%      9.8 * params.m * initZ(2, :);...
%      (1-initZ)];
if (initialization)
    params.Delta0 = [params.X0Init(:,1:end-1); params.U0Init];
else
    params.Delta0=zeros(params.dim , params.N-1);
end
params.P0= zeros(params.dim, (params.N-1));

params.fn1Ndx = params.nx + 4;
params.fn2Ndx = params.nx + 8;
params.phi1Ndx = 6;
params.phi2Ndx = 8;
params.phiDot1Ndx = 10;
params.phiDot2Ndx = 12;
params.phiDDot1Ndx = params.nx + 10;
params.phiDDot2Ndx = params.nx + 12;
params.vT1Ndx = 9;
params.vT2Ndx = 11;
end