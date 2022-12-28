function [] = plotSeparationAlt(intVars, Z, Delta, params)
%UNTITLED2 Summary of this function goes here

h = params.dt;
timeTraj = 0:params.dt:params.horizon;

figure(4)
subplot(2,4,1)
plot(timeTraj(1:end-1), Z(params.fn1Ndx, :));
hold on;
plot(timeTraj(1:end-1), 100 * Z(params.vT1Ndx, :));
%plot(timeTraj, 100 * Delta(6, :),'x');
for k = 1 : size(Delta, 2)
    Delta_k = Delta(:, k);
    plot(timeTraj(k:k+1), 100 * [Delta_k(params.vT1Ndx), Delta_k(params.vT1Ndx) + h * Delta_k(params.vDotT1Ndx)],'--*')
end
plotAreas(timeTraj, intVars(1,:));
hold off
xlabel('time')
title('Contact 1 CC')
legend('fN1 from Z_k','100 * vT1 from Z_k', '100 * vT1(k) from Δ_k');


subplot(2,4,2)
plot(timeTraj(1:end-1), Z(params.fn2Ndx, :));
hold on;
plot(timeTraj(1:end-1), 100 * Z(params.vT2Ndx, :));
%plot(timeTraj, 100 * Delta(6, :),'x');
for k = 1 : size(Delta, 2)
    Delta_k = Delta(:, k);
    plot(timeTraj(k:k+1), 100 * [Delta_k(params.vT2Ndx), Delta_k(params.vT2Ndx) + h * Delta_k(params.vDotT2Ndx)],'--*')
end
plotAreas(timeTraj, intVars(2,:));
hold off
xlabel('time')
title('Contact 2 CC')
legend('fN2 from Z_k','100 * vT2 from Z_k', '100 * vT2(k) from Δ_k');

subplot(2,4,3)
plot(timeTraj(1:end-1), Z(params.fn1Ndx, :));
hold on;
plot(timeTraj(1:end-1), Delta(params.fn1Ndx, :),'--*');
plotAreas(timeTraj, intVars(1,:));
hold off
xlabel('time')
title('fN1 Z_k vs Δ_k')
legend('fN1 from Z_k','fN1 from Δ_k');

subplot(2,4,4)
plot(timeTraj(1:end-1), Z(params.fn2Ndx, :));
hold on;
plot(timeTraj(1:end-1), Delta(params.fn2Ndx, :),'--*');
plotAreas(timeTraj, intVars(2,:));
hold off
xlabel('time')
title('fN2 Z_k vs Δ_k')
legend('fN2 from Z_k','fN2 from Δ_k');

subplot(2,4,5)
plot(timeTraj(1:end-1), Z(params.vDotT1Ndx, :));
hold on;
plot(timeTraj(1:end-1), Delta(params.vDotT1Ndx, :));
plotAreas(timeTraj, intVars(1,:));
hold off
xlabel('time')
title('vDotT1 Z_k vs Δ_k')
legend('vDotT1 from Z_k','vDotT1 from Δ_k');

subplot(2,4,6)
plot(timeTraj(1:end-1), Z(params.vDotT2Ndx, :));
hold on;
plot(timeTraj(1:end-1), Delta(params.vDotT2Ndx, :));
plotAreas(timeTraj, intVars(2,:));
hold off
xlabel('time')
title('vDotT2 Z_k vs Δ_k')
legend('vDotT2 from Z_k','vDotT2 from Δ_k');

end
function [] = plotAreas(timeTraj, intVars)
YL = get(gca, 'YLim');
ind = 1:length(intVars);
for i = ind(logical(intVars))
    x_points = [timeTraj(i), timeTraj(i), timeTraj(i+1), timeTraj(i+1)];  
    y_points = [YL(1), YL(2), YL(2), YL(1)];
    color = [0, 0, 1];
    a = fill(x_points, y_points, color);
    a.FaceAlpha = 0.1;
    a.EdgeColor = "none";
end
end