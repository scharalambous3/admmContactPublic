function [] = plotPerf(hax, orthViol,objValueVec, primalResidualVec, dualResidualVec, intVars, xTraj)
%UNTITLED2 Summary of this function goes here
figure(2)
subplot(4,2,1)
semilogy(1:size(orthViol, 2), orthViol);
xlabel('Iterations')
ylabel('Orthogonality violation')
ylim([0, inf])

subplot(4,2,2)
semilogy(1:size(objValueVec, 2), objValueVec);
xlabel('Iterations')
ylabel('Objective value')

subplot(4,2,3)
semilogy(1:size(primalResidualVec, 2), primalResidualVec);
xlabel('Iterations')
ylabel('z delta violation')
ylim([0, inf])
title('Prima Residual')

subplot(4,2,4)
semilogy(1:size(dualResidualVec, 2), dualResidualVec);
xlabel('Iterations')
ylabel('Dual residual')
ylim([0, inf])
title('Dual Residual')

subplot(4,2,5:6)
mat = 1-intVars;  % Your sample matrix
[r, c] = size(mat);                          % Get the matrix size
imagesc((1:c)+0.5, (1:r)+0.5, mat);          % Plot the image
colormap(gray);                              % Use a gray colormap
set(gca, 'XTick', 1:(c+1), 'YTick', 1:(r+1), ...  % Change some axes properties
         'XLim', [1 c+1], 'YLim', [1 r+1], ...
         'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');

subplot(4,2,7:8)
tRange=[1:length(xTraj)];
plot(tRange, xTraj);
xlabel('time')
ylabel('x')
end