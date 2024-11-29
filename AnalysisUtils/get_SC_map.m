function get_SC_map()
% Define parameters
Bx = 1.1;
By = 1.8;
A = 0.9;

% Create a grid of polar coordinates
R = linspace(0, 60, 100); % radial distance from 0 to 30 degrees for a wider range
deg = [-90, -60, -30, 0 30 60 90];
theta = deg2rad(deg); % angle from 0 to 2*pi radians
% [R, theta] = meshgrid(R, theta); % create a mesh grid

X = [];
Y = [];
for i = 1:length(R)
    for j = 1:length(theta)
        % Convert polar coordinates to cartesian coordinates using the given equations
        X(i, j) = Bx * log(sqrt(R(i)^2 + 2 * A * (R(i) * cos(theta(j)))  + A ^ 2) / A);
        Y(i, j) = By * atan((R(i) .* sin(theta(j))) / (R(i) * cos(theta(j)) + A));
    end
end

plot(X, Y, 'k')

% Annotate the plot with theta values
for j = 1:length(theta)
    text(X(end, 1)+ 0.2, Y(end, j), sprintf('%d°', deg(j)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
end

% Define parameters
Bx = 1.1;
By = 1.8;
A = 0.9;

% Create a set of fixed R values
R = [0.5, 1, 2, 5, 10, 20, 40, 60]; % example R values
theta = linspace(-pi/2, pi/2, 100); % angle from -90 to 90 degrees

X = [];
Y = [];
for i = 1:length(R)
    for j = 1:length(theta)
        % Convert polar coordinates to cartesian coordinates using the given equations
        X(i, j) = Bx * log(sqrt(R(i)^2 + 2 * A * (R(i) * cos(theta(j)))  + A ^ 2) / A);
        Y(i, j) = By * atan((R(i) * sin(theta(j))) / (R(i) * cos(theta(j)) + A));
    end
end

% Plot the lines for different Rs
hold on;
for i = 1:length(R)
    plot(X(i, :), Y(i, :), 'k');
end

% Annotate the plot with R values
for i = 1:length(R)
    text(X(i, end)-0.2/i, Y(i, end)+0.1*i/10, [num2str(R(i)), '°'], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end

hold off;
title('SC Map Chen et al. 2019');
xlabel('X (mm)');
ylabel('Y (mm)');
axis equal;
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
box off
ylim([-3, 3])
xlim([-0.5, 5])

end