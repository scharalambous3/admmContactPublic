function [] = plotSeparation(intVars, Z, Delta, params)
%UNTITLED2 Summary of this function goes here
if ~isfield(params,'phi1Ndx')
    plotSeparationAlt(intVars, Z, Delta, params)
    return;
end
h = params.dt;
timeTraj = 0:params.dt:params.horizon;

figure(4)
subplot(2,4,1)
plot(timeTraj(1:end-1), Z(params.fn1Ndx, :));
hold on;
plot(timeTraj(1:end-1), 100 * Z(params.phi1Ndx, :));
%plot(timeTraj, 100 * Delta(6, :),'x');
for k = 1 : size(Delta, 2)
    Delta_k = Delta(:, k);
    plot(timeTraj(k:k+1), 100 * [Delta_k(params.phi1Ndx), Delta_k(params.phi1Ndx) + h * Delta_k(params.phiDot1Ndx) + h^2 * Delta_k(params.phiDDot1Ndx)],'--*')
end
plotAreas(timeTraj, intVars(1,:));
hold off
xlabel('time')
title('Contact 1 Separation Constraint')
legend('fN1 from Z_k','100 * phi1 from Z_k', '100 * phi1(k) from Δ_k');


subplot(2,4,2)
plot(timeTraj(1:end-1), Z(params.fn2Ndx, :));
hold on;
plot(timeTraj(1:end-1), 100 * Z(params.phi2Ndx, :));
%plot(timeTraj, 100 * Delta(6, :),'x');
for k = 1 : size(Delta, 2)
    Delta_k = Delta(:, k);
    plot(timeTraj(k:k+1), 100 * [Delta_k(params.phi2Ndx), Delta_k(params.phi2Ndx) + h * Delta_k(params.phiDot2Ndx) + h^2 * Delta_k(params.phiDDot2Ndx)],'--*')
end
plotAreas(timeTraj, intVars(2,:));
hold off
xlabel('time')
title('Contact 2 Separation Constraint')
legend('fN2 from Z_k','100 * phi2 from Z_k', '100 * phi2(k) from Δ_k');

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
plot(timeTraj(1:end-1), Z(params.phiDot1Ndx, :));
hold on;
plot(timeTraj(1:end-1), Delta(params.phiDot1Ndx, :));
plotAreas(timeTraj, intVars(1,:));
hold off
xlabel('time')
title('phidot1 Z_k vs Δ_k')
legend('phidot1 from Z_k','phidot1 from Δ_k');

subplot(2,4,6)
plot(timeTraj(1:end-1), Z(params.phiDot2Ndx, :));
hold on;
plot(timeTraj(1:end-1), Delta(params.phiDot2Ndx, :));
plotAreas(timeTraj, intVars(2,:));
hold off
xlabel('time')
title('phidot2 Z_k vs Δ_k')
legend('phidot2 from Z_k','phidot2 from Δ_k');


subplot(2,4,7)
plot(timeTraj(1:end-1), Z(params.phiDDot1Ndx, :));
hold on;
plot(timeTraj(1:end-1), Delta(params.phiDDot1Ndx, :));
plotAreas(timeTraj, intVars(1,:));
hold off
xlabel('time')
title('phiddot1 Z_k vs Δ_k')
legend('phiddot1 from Z_k','phiddot1 from Δ_k');

subplot(2,4,8)
plot(timeTraj(1:end-1), Z(params.phiDDot2Ndx, :));
hold on;
plot(timeTraj(1:end-1), Delta(params.phiDDot2Ndx, :));
plotAreas(timeTraj, intVars(2,:));
hold off
xlabel('time')
title('phiddot2 Z_k vs Δ_k')
legend('phiddot2 from Z_k','phiddot2 from Δ_k');

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