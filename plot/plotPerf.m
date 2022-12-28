function [] = plotPerf(hax, orthViol,objValueVec, primalResidualArr, dualResidualArr, intVars, xTraj)
%plotPerf Plots metrics (orthogonality violation, objective value of given
%state trajectory, primal and dual residuals), trajectory of state
%variables and boolean assignment of the integer decision variables

figure(2)
subplot(5,2,1)
semilogy(1:size(orthViol, 2), orthViol);
xlabel('Iterations')
ylabel('Orthogonality violation')
ylim([0, inf])

subplot(5,2,2)
semilogy(1:size(objValueVec, 2), objValueVec);
xlabel('Iterations')
ylabel('Objective value')

subplot(5,2,3)
semilogy(1:size(primalResidualArr, 2), vecnorm(primalResidualArr, 2 , 1));
xlabel('Iterations')
ylabel('z delta violation')
ylim([0, inf])
title('Primal Residual')

subplot(5,2,4)
semilogy(1:size(dualResidualArr, 2), vecnorm(dualResidualArr, 2, 1));
xlabel('Iterations')
ylabel('Dual residual')
ylim([0, inf])
title('Dual Residual')

subplot(5,2,5)
semilogy(1:size(primalResidualArr, 2), primalResidualArr);
xlabel('Iterations')
ylabel('z delta violation')
ylim([0, inf])
title('Primal Residual')

subplot(5,2,6)
semilogy(1:size(dualResidualArr, 2), dualResidualArr);
xlabel('Iterations')
ylabel('Dual residual')
ylim([0, inf])
title('Dual Residual')

subplot(5,2,7:8)
mat = 1-intVars;  % Your sample matrix
[r, c] = size(mat);                          % Get the matrix size
imagesc((1:c)+0.5, (1:r)+0.5, mat);          % Plot the image
colormap(gray);                              % Use a gray colormap
set(gca, 'XTick', 1:(c+1), 'YTick', 1:(r+1), ...  % Change some axes properties
         'XLim', [1 c+1], 'YLim', [1 r+1], ...
         'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');

subplot(5,2,9:10)
tRange=[1:length(xTraj)];
plot(tRange, xTraj);
xlabel('time')
ylabel('x')
end