function [params] = getGaitingParams()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
params.dt = 0.1;
params.N = 10;
params.horizon = params.N * params.dt;

params.convexSubproblemSettings = sdpsettings('solver','osqp','cachesolvers',1,'allownonconvex',0, 'osqp.time_limit', 0.01);

params.finalTime = 3.0;
params.simSteps = params.finalTime/params.dt;

params.nx = 6;
params.nu = 10; % Input = [6 dimensional lambda; 4 dimensional u]
params.dim = params.nx + params.nu;
params.animation = false;
params.liveGraphs = true;

params.X0 =  [-8;0;3;0;4;0];

params.xDes= zeros(params.nx, 1);
params.NLInitialization=0;
%Since dynamics are linear now, Ts doesnt matter. Once I have NL dynamics
%Ts can be smaller than dt for better integration
g = 9.81;
mu = 1;
h = params.dt;
A = [[1,h,0,0,0,0];[0, 1, 0, 0, 0, 0]; [0, 0, 1, h, 0, 0 ]; [0, 0, 0, 1, 0, 0]; [0, 0, 0, 0, 1, h]; [0, 0, 0, 0, 0, 1]];
B = [[0,0,0,0]; [0, 0, 0, 0]; [h*h, 0, 0, 0]; [h, 0, 0, 0]; [0, h*h, 0, 0]; [0, h, 0, 0]];
D = [[0, h*h, -h*h, 0, h*h, -h*h]; [0, h, -h, 0, h, -h]; [0, -h*h, h*h, 0, 0, 0]; [0, -h, h, 0, 0, 0]; [0, 0, 0, 0, -h*h, h*h]; [0, 0, 0, 0, -h, h]];
E = [[0, 0, 0, 0, 0, 0]; [0, 1, 0, -1, 0, 0]; [0, -1, 0, 1, 0, 0]; [0, 0, 0, 0, 0, 0]; [0, 1, 0, 0, 0, -1]; [0, -1, 0, 0, 0, 1]];
F = [[0, -1, -1, 0, 0, 0]; [1, 2*h, -2*h, 0, h, -h]; [1, -2*h, 2*h, 0, -h, h]; [0, 0, 0, 0, -1,-1]; [0, h, -h, 1, 2*h, -2*h]; [0, -h, h, 1, -2*h, 2*h]];
c = [[0];[-h*g]; [h*g]; [0]; [-h*g]; [h*g]];
d = [[-g*h*h];[-g*h];[0];[0];[0];[0]];
H = [[0, 0, mu, 0]; [-h, 0, 0, 0]; [h, 0, 0, 0]; [0, 0, 0, mu]; [0, -h, 0, 0]; [0, h, 0, 0]];


params.A = A;
params.B = [D, B]; % Since input in my case is [lambda; u]
params.d = d;

params.M=1000;

params.Q = [[5000, 0, 0, 0, 0, 0]; [0, 10, 0, 0, 0, 0]; [0, 0, 10, 0, 0 ,0]; [0, 0, 0, 10, 0, 0]; [0, 0, 0, 0, 10, 0]; [0, 0, 0, 0, 0, 10]];
params.R = diag([zeros(1, 6), ones(1,4) ]);
%params.Qf = idare(params.A, params.B, params.Q, params.R,[],[]);
params.Qf = params.Q;



%ADMM
params.epsilon0=inf;
params.epsDyn = 1e-16;
params.rho = 1;
params.rhoScale = 1.2;
params.maxIters=10;
params.projG_k = blkdiag(1000*eye(params.nx), eye(params.nu));
params.G0 = (params.rho/2) * eye(params.dim, params.dim);

% Constraints Ax * X >= bx and Au * U >= bu
%u3 and u4 >= 0
params.Au = zeros(2, params.nu);
params.Au(:, end-1:end) = eye(2);
params.bu = zeros(2, 1);

params.Ax = zeros(4, params.nx);
params.Ax(1:2,[3,5]) = eye(2);
params.Ax(3:4,[3,5]) = - eye(2);
params.bx = [1;3;-3;-5];

params.AxTerminal = [];
params.bxTerminal = [];

params.orthDim = 6;

params.Aorth = zeros(params.orthDim, params.dim);
params.Aorth(:, params.nx+1:params.nx+params.orthDim) = eye(params.orthDim);
params.aorth = zeros(params.orthDim, 1);
params.Borth = [E, F, H];
params.borth = c;


%Constraint Adelta_delta * Delta + Adelta_int * int >= bdelta
params.Adelta_delta = [params.Aorth; params.Borth; params.Aorth; -params.Aorth; params.Borth; - params.Borth];

params.Adelta_int = [zeros(params.orthDim);... %Aorth * Z + aorth >= 0
    zeros(params.orthDim);...  %%Borth * Z + borth >= 0
    params.M * eye(params.orthDim);... % - M*i <= Aorth * Z + aorth redundant in cartpole
    params.M * eye(params.orthDim);... % Aorth * Z  + aorth <= + M*i 
    -params.M * eye(params.orthDim);... % -M*(1-i) <= Borth * Z + borth
    -params.M * eye(params.orthDim)]; % Borth * Z + borth <= M*(1-i)

params.bdelta = [- params.aorth; - params.borth; - params.aorth; params.aorth; - params.M * ones(params.orthDim, 1) - params.borth; - params.M * ones(params.orthDim, 1) + params.borth];

%Initialization
params.Delta0=zeros(params.dim, params.N-1);
params.P0 = zeros(params.dim, (params.N-1));
end