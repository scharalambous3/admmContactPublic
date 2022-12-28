function [params] = getBipedParams()
%getBipedParams Get params for 2D floating base system with 2 contacts,
%linear dynamics, quadratic state-input cost
%Separation and sliding complementarity constraints 

params.dt = 0.05;
params.N = 51;
params.horizon = (params.N - 1) * params.dt;
params.NLInitialization=0;
%params.convexSubproblemSettings = sdpsettings('solver','snopt','cachesolvers',1,'allownonconvex',1, 'usex0', params.NLInitialization);%, 'snopt.Iterations_limit', 500);%, 'osqp.time_limit', 0.01);
params.convexSubproblemSettings = sdpsettings('solver','mosek','cachesolvers',1,'allownonconvex',0);%, 'snopt.Iterations_limit', 500);%, 'osqp.time_limit', 0.01);
params.finalTime = 3.0;
params.simSteps = params.finalTime/params.dt;

params.nx = 6;
params.nu = 6;
params.dim = params.nx + params.nu;
params.animation = false;
params.liveGraphs = true;

params.X0 = [0; 0.5; 0; 0; -0.25; 0.25];
xDesFinal= [params.finalTime * 0.5; 0.5; 0; 0; 0; 0];

params.xDes = repmat(xDesFinal, 1, params.N);
%ADMM
params.epsDyn = 1e-16;
%params.rho = 0.2/1000;
params.rho = 1;
params.rhoScale = 1.2;
params.maxIters=30;
params.epsilon0 = 100;
%The state, fx1 and fx2 do not enter the orthogonality constraint
%In the projection subproblem I need nonzero entries in G to retun non-NaNs
params.projG_k = blkdiag(eye(params.nx), eye(4), 1000*eye(2));
%The elements of omega and delta corresponding to these are irrelevant by
%setting the respective weights in G to 0
params.G0 = (params.rho/2) * blkdiag(zeros(params.nx), diag([0, 1, 0, 1, 100000, 100000]));

params.g = [0; -9.8];
initialization = false;
params.m = 50.0;
params.I = (params.m/12) * (0.75^2 + 0.5^2);

params.Q = diag([50000, 50000, 50000, 50000, 0, 0]);
params.R = diag([0.001, 0.0001, 0.001, 0.0001, 1, 1]);
%params.Qf = idare(params.A, params.B, params.Q, params.R,[],[]);
params.Qf = params.Q;

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

%rci box constraint Ax * X >= bx
%     X(1, i) - 0.25 - 0.1 <= X(5, i), X(5, i) <= X(1, i) - 0.25 + 0.1
%     X(1, i) + 0.25 - 0.1 <= X(6, i), X(6, i) <= X(1, i) + 0.25 + 0.1
%     X(2, i) >= cyDes - 0.15
%     X(2, i) <= cyDes + 0.15 --> - (cyDes + 0.15) <= - X(2, i)
params.Ax = [-1, 0, 0, 0, 1, 0;...
             1, 0, 0, 0,-1, 0;...
             1, 0, 0, 0, 0, -1;...
            -1, 0, 0, 0, 0, 1;...
             0, 1, 0, 0, 0, 0;...
             0, -1, 0, 0, 0, 0];
params.bx = [-0.35; 0.15; -0.35; 0.15; params.xDes(2, 1) - 0.15; - (params.xDes(2, 1) + 0.15)];

params.mu=0.9;
%unilateral force constraint Au U >= bu.  
params.Au = zeros(6, params.nu);
params.Au(1:2,[2, 4]) = eye(2); % f2>=0, f4>=0
params.Au(3:6,1:4) = blkdiag([-1, params.mu; 1,params.mu],[-1, params.mu; 1,params.mu]); % friction cone
params.bu = zeros(6, 1);

params.AxTerminal = [];
params.bxTerminal = [];

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

params.X0Init = [[linspace(0, params.xDes(1, end), params.N); 0.5 * ones(1,params.N)];...
    [(params.xDes(1, end)/params.horizon) * ones(1, params.N); zeros(1,params.N)];...
    [linspace(0, params.xDes(1, end), params.N) - 0.25; linspace(0, params.xDes(1), params.N) + 0.25]];
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