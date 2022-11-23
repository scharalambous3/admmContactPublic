%%
clear
close all
clc
yalmip('clear')
%%
params = getBipedParams();

x_k = params.X0;

dim = params.dim;
Delta0 = params.Delta0;
P0 = params.P0;
G0 = (params.rho/2) * eye(params.dim, dim);

X=[];
U=[];
for k = 1:(params.simSteps - 1)
    [Z_k,Delta_k, intVars, orthViol, objValueVec, zdeltaViolVec] = solveADMM(Delta0, P0, G0, x_k, params);
    u_k = Z_k((params.nx + 1):end, 1);
    x_k = discreteDyn(x_k, u_k, params);
    X=[X, x_k];
    U=[U, u_k];
end


%%
tRange=[0:params.dt:size(X, 2) * params.dt];
figure(1)
plot(tRange(1:end-1), X);
stateDim=[1:1:params.nx];
legend (sprintfc('%.1f',stateDim))

%%
if (params.animation)
    hax = figure(3);
    
    %gif('A.gif','DelayTime',1/5)
    
    while (true)
        for k=1:length(tRange)-1
            drawFloatingBase(hax, c(:,k),f(:,k), xDes(1), [rci(1,k);0.0], [rci(2,k);0.0]);
        end
    
    end
    
    
end

function rotMatrix = rotz(theta)
rotMatrix = [cos(theta), - sin(theta), 0;...
             sin(theta), cos(theta), 0;...
             0, 0, 1];

end