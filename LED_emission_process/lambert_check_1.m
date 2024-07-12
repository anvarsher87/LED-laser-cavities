% Clear all variables and close all figures
clear;
close all;

% Initialize variables
u = [];
v = [];
w = [];

count = [];
bin_theta = [];
nbins = 45;
delta_bin = 90 / nbins;

% Surface properties
tang1 = [1, 0, 0];
tang2 = [0, 1, 0];
norm = [0, 0, 1];

% Initialize data
for bin = 0:nbins-1
    count(bin + 1) = 0;
    bin_theta(bin + 1) = bin * delta_bin;
end

% Sample random points
for it = 1:1000000
    sin_theta = sqrt(rand())*sin(pi/2);
    cos_theta = sqrt(1 - sin_theta^2);

    % Random in-plane angle
    psi = rand() * 2 * pi;

    % Three vector components
    a = sin_theta * cos(psi);
    b = sin_theta * sin(psi);
    c = cos_theta;

    % Multiply by corresponding directions
    v1 = a * tang1;
    v2 = b * tang2;
    v3 = c * norm;

    % Add up to get velocity, vel = v1 + v2 + v3
    vel = v1 + v2 + v3;

    % Update histogram
    theta = acos(cos_theta) * 180 / pi;
    bin = floor(theta / delta_bin) + 1;
    count(bin) = count(bin) + 1;

    % Add every 1000th particle to the visualization list
    if mod(it, 1000) == 0
        u = [u, vel(1)];
        v = [v, vel(2)];
        w = [w, vel(3)];
    end
end

% Plot results
figure;
subplot(2, 1, 1);
scatter3(u, v, w, 'r.');
xlabel('u');
ylabel('v');
zlabel('w');

% Divide by bin area
for i = 1:nbins
    t1 = (i - 1) * delta_bin * pi / 180;
    t2 = i * delta_bin * pi / 180;
    A = 2 * pi * ((1 - cos(t2)) - (1 - cos(t1)));
    count(i) = count(i) / A;
end

% Normalize data
c0 = 0.5 * (count(1) + count(2));
count = count / c0;
cos_theta_plot = cosd(bin_theta);

% Add xy plot
subplot(2, 1, 2);
plot(bin_theta, count, 'r');
hold on;
plot(bin_theta, cos_theta_plot, 'k');
xlabel('angle');
ylabel('normalized count');
hold off;

figure
subplot(1,2,1)
polarplot(deg2rad(bin_theta), count, 'r');
thetalim([0 90])
set(gca, 'ThetaDir', 'counterclockwise', 'ThetaZeroLocation', 'top')
subplot(1,2,2)
plot(bin_theta, count, 'r');
xlabel('angle');
ylabel('normalized count');
hold off;


