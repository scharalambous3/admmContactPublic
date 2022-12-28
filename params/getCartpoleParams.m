function [params] = getCartpoleParams()
%getCartpoleParams Params for cart-pole with soft walls

params.dt = 0.01;
params.N = 11;
params.horizon = (params.N - 1)* params.dt;
params.epsilon0=1000000;
params.finalTime = 1.0;
params.simSteps = params.finalTime/params.dt;
params.convexSubproblemSettings = sdpsettings('solver','osqp','cachesolvers',1,'allownonconvex',0,'osqp.time_limit', 0.1);%, 'snopt.Iterations_limit', 500);%, 'osqp.time_limit', 0.01);
params.NLInitialization=0;
params.separationIndices=[1,2];
params.nx = 4;
params.nu = 3;
params.dim = params.nx + params.nu;
params.animation = false;
params.liveGraphs = true;

xDesFinal= [0; 0; 0; 0];
params.xDes = repmat(xDesFinal, 1, params.N);
params.X0 =[0.3; 3; 0.3; 0];

g = 9.81;
mp = 0.411; 
mc = 0.978; 
len_p = 0.6; 
len_com = 0.4267; 
d1 = 0.35;
d2 = -0.35;
ks= 100;
Ts = 0.01;
A = [[0, 0, 1, 0]; [0, 0, 0, 1]; [0, g*mp/mc, 0, 0]; [0, g*(mc+mp)/(len_com*mc), 0, 0]];
B = [0;0;1/mc;1/(len_com*mc)];
D = [[0,0]; [0,0]; [(-1/mc) + (len_p/(mc*len_com)), (1/mc) - (len_p/(mc*len_com)) ]; [(-1 / (mc*len_com) ) + (len_p*(mc+mp)) / (mc*mp*len_com*len_com)  , -((-1 / (mc*len_com) ) + (len_p*(mc+mp)) / (mc*mp*len_com*len_com))]];
E = [[-1, len_p, 0, 0]; [1, -len_p, 0, 0 ]];
F = 1/ks * eye(2);
d = zeros(4,1);
H = zeros(2,1);
c = [d1; -d2];
params.A = eye(params.nx) + Ts * A;
Bu = Ts*B;
Dlambda = Ts*D;
params.d = Ts*d;

params.B = [Dlambda, Bu];

params.M=1000;

params.Q = 2* [[10, 0, 0, 0]; [0, 3, 0, 0]; [0, 0, 1, 0]; [0, 0, 0, 1]];
params.R = 2* diag([0, 0, 1]);
params.Qf = idare(params.A, Bu, params.Q, 1,[],[]); % 1 is for the u part only (ie w/o lambda part)


%ADMM
params.epsDyn = 1e-16;
params.rho = 0.2; %to meet 0.2/2 = 0.1 of Posa's code
params.rhoScale = 2;
params.maxIters=10;

params.projG_k = blkdiag(1000 * eye(params.nx), eye(2), eye(1));
params.G0 = 2 * (params.rho/2) * eye(params.dim, params.dim);

params.Ax = [];
params.bx = [];

params.AxTerminal = [];
params.bxTerminal = [];

%unilateral force constraint Au U >= bu.  f2>=0, f4>=0
params.Au = [];
params.bu = [];

params.orthDim = 2;

params.Aorth = zeros(params.orthDim, params.dim);
params.Aorth(:, 5:6) = eye(2);
params.aorth = zeros(2, 1);
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