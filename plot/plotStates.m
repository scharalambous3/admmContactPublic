function [] = plotStates(intVars, Z, params)
%plotStates Plots state variables involved in complementarity constraints

timeTraj = 0:params.dt:params.horizon;
figure(3)
subplot(2,4,1)
plot(timeTraj(1:end-1), Z(1, :));
hold on;
plot(timeTraj(1:end-1), Z(2, :));
yline(params.xDes(1, 1), '--');
hold off
xlabel('time')
title('CoM position')
legend('cx','cy');

if ~isfield(params,'fn1Ndx')
    return;
end
subplot(2,4,4)
plot(timeTraj(1:end-1), Z(params.fn1Ndx, :));
hold on;
plot(timeTraj(1:end-1), 100 * Z(params.vT1Ndx, :));
plotAreas(timeTraj, intVars(1,:));
hold off
xlabel('time')
title('Contact 1 Sliding')
legend('fN1','100 * vT1');

subplot(2,4,5)
plot(timeTraj(1:end-1), Z(params.fn2Ndx, :));
hold on;
plot(timeTraj(1:end-1), 100 * Z(params.vT2Ndx, :));
plotAreas(timeTraj, intVars(2,:));
hold off
xlabel('time')
title('Contact 2 Sliding')
legend('fN2','100 * vT2');

if ~isfield(params,'phi1Ndx')
    return;
end

subplot(2,4,2)
plot(timeTraj(1:end-1), Z(params.fn1Ndx, :));
hold on;
plot(timeTraj(1:end-1), 100 * Z(params.phi1Ndx, :));
plotAreas(timeTraj, intVars(1,:));
hold off
xlabel('time')
title('Contact 1 Separation Constraint')
legend('fN1','100 * phi1');

subplot(2,4,3)
plot(timeTraj(1:end-1), Z(params.fn2Ndx, :));
hold on;
plot(timeTraj(1:end-1), 100 * Z(params.phi2Ndx, :));
plotAreas(timeTraj, intVars(2,:));
hold off
xlabel('time')
title('Contact 2 Separation Constraint')
legend('fN2','100 * phi2');
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