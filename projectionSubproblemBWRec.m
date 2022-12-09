function [Delta_k, int_k] = projectionSubproblem(Z_k, P_k, params)
%UNTITLED4 Summary of this function goes here
% Detailed explanation goes here
N = params.N;
Delta_k = [];
int_k = [];
ops = sdpsettings('solver','mosek','cachesolvers',1,'verbose',0);

%delta = [1cx, 2cy, 3cdotx, 4cdoty, 5rc1, 6rc2, 7fx1, 8fy1, 9fx2, 10fy2, 11rdotc1, 12rdotc2]
%delta = [delta_c∈R2;delta_cdot∈R2;delta_rci∈R2;delta_f∈R4;delta_rdotci∈R2]

G_k = params.projG_k;
%For cartpole
%G_k(end, end) = 0;

rdot1_kplusone = zeros(2,1);
r1y_kplusone = 0;
rdot2_kplusone = zeros(2,1);
r2y_kplusone = 0;
for i = (N - 1):-1:1
    delta = sdpvar(params.dim, 1);
    intVar = binvar(params.orthDim, 1);
 
    %consensus cost
    obj = (Z_k(:, i) - delta + P_k(:, i))' * G_k * (Z_k(:, i) - delta + P_k(:, i));

    constr = [];
    % orthogonality constraints:  f>= 0, rdotci orthogonal to f

    %TODO: constraints on integer variables
    %constr = [constr, intVar(1) + intVar(2) >= 1];
    v1_kplusone = [1 0] * rdot1_kplusone;
    v2_kplusone = [1 0] * rdot2_kplusone;
    phi1_kplusone = r1y_kplusone;
    phi2_kplusone = r2y_kplusone;
    borth = [0; v1_kplusone'; -v1_kplusone'; phi1_kplusone; 0; v2_kplusone'; -v2_kplusone'; phi2_kplusone];
    bdelta = [- params.aorth; - borth; - params.aorth; params.aorth; - params.M * ones(params.orthDim, 1) - borth; - params.M * ones(params.orthDim, 1) + borth];

    % zero force during swing
    constr = [constr, params.Adelta_delta * delta + params.Adelta_int * intVar >= bdelta];

    constr = [constr,...
        [r1y_kplusone;r2y_kplusone;rdot1_kplusone;rdot2_kplusone] == [params.A([6,8,9,10,11,12], :) params.B([6,8,9,10,11,12], :)] * delta + params.d([6,8,9,10,11,12])];
    
    optimize(constr,obj,ops);
    Delta_k = [value(delta), Delta_k];
    int_k = [value(intVar), int_k];

    r1y_kplusone = value(delta(6));
    r2y_kplusone = value(delta(8));
    rdot1_kplusone = value(delta(9:10));
    rdot2_kplusone = value(delta(11:12));
end
%For biped
%Delta_k(1:params.nx,:)= Z_k(1:params.nx,:);

%For cartpole
%Delta_k(1:params.nx,1) = 0;
%Delta_k(7,isnan(Delta_k(7,:))) = 0; %since u is not part of objective nor of constraints it returns NaN in cartpole

assert(all(size(Delta_k) == [params.dim, N - 1]))
end

