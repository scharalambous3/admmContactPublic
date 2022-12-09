function [Z_k,Delta_k, int_k, violVec, objValueVec, primalResidualVec, dualResidualVec] = solveADMM(Delta0, P0, G0, x_k, params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Delta_k = Delta0;
Z_k = [];
P_k = P0;
G_k = G0;
epsilon = params.epsilon0;
violVec = [];
objValueVec =[];

primalResidualVec=[];
dualResidualVec=[];

prevZ_k = [];
int_k = [];

PVec = [];
viol = inf;
for iter = 1:params.maxIters
    if (viol <= params.epsDyn)
       fprintf('Exited at iteration %i. Orthogonality violation norm: %d. Rho: %d \n', iter-1, violNorm, rho) 
       break; 
    end
    % For cartpole
    %G_k(params.dim, params.dim) = 0;

    %For cartpole (?)
    %P_k(1:params.nx,1) = 0; %Dual should have no effect on initial state
    
    [Z_k, objValue] = solveConvexSubproblem(Delta_k, P_k, G_k, x_k, params, epsilon);

    viol = getOrthogonalityViolation(Z_k, params);

    prevInt_k = int_k;
    [Delta_k, int_k] = projectionSubproblemGrouping(Z_k, P_k, params, prevInt_k);

    [xTraj, uTraj, rolloutObjValue] = getRollout(Z_k, x_k, params);
    objValueVec =[objValueVec, rolloutObjValue];

    primalResidualVec = [primalResidualVec, norm(vecnorm(Z_k - Delta_k,2,1))];
%    primalResidualVec = [primalResidualVec, norm(Z_k(:) - Delta_k(:))];
    %dual residual = rho * I * (-I) * (z_k+1 - z+k)
    if (iter > 1)
        dualResidualVec = [dualResidualVec, norm(params.rho * (Z_k(:) - prevZ_k(:)))];
    end
    prevZ_k = Z_k;
    P_k = P_k + Z_k - Delta_k;
    P_k = P_k / params.rhoScale;
    G_k = G_k * params.rhoScale;
    fprintf('Iteration %i. Orthogonality violation norm: %d. Rho: %d \n', iter, viol, norm(G_k)) 
    violVec = [violVec, viol];
    epsilon = epsilon/params.rhoScale;

    if (params.liveGraphs)
        figure(2)
        hax = gca;
        plotPerf(hax, violVec,objValueVec, primalResidualVec, dualResidualVec, int_k, xTraj)
    end
end


end

