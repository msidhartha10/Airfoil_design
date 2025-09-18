% NACA 5 and 6 Series Airfoil Plot Generator with Cosine Spacing and Excel Export
% Requirements: MATLAB, xlswrite (built-in), or writematrix (newer MATLAB)

clear; clc;

% ---- User Parameters ----
naca_type = 5; % 5 for NACA 5-series, 6 for NACA 6-series
naca_code = '23012'; % Example: '23012' for 5-series, '641212' for 6-series
num_points = 100;    % Number of coordinate points (cosine spaced)
chord = 1.0;         % Chord length
output_excel = 'airfoil_coords.xlsx'; % Output Excel file name
% ------------------------

% Generate cosine spaced x-coordinates
beta = linspace(0, pi, num_points);
x = (1 - cos(beta))/2 * chord;

if naca_type == 5
    [xu, yu, xl, yl, xc, yc] = naca5digit(naca_code, x, chord);
elseif naca_type == 6
    [xu, yu, xl, yl, xc, yc] = naca6digit(naca_code, x, chord);
else
    error('Only 5 or 6 series supported.');
end

% Plot
figure;
plot(xu, yu, 'b', xl, yl, 'r', xc, yc, 'k--', 'LineWidth', 1.5);
legend('Upper Surface', 'Lower Surface', 'Camber Line');
axis equal; grid on;
xlabel('x'); ylabel('y');
title(['NACA ', naca_code, ' airfoil']);

% Prepare data for Excel
T = table(xu', yu', xl', yl', xc', yc', ...
    'VariableNames', {'Upper_X','Upper_Y','Lower_X','Lower_Y','Camber_X','Camber_Y'});
writetable(T, output_excel);
disp(['Exported coordinates to: ', output_excel]);

% --- NACA 5-digit function ---
function [xu, yu, xl, yl, xc, yc] = naca5digit(code, x, c)
    % Parse code
    m = str2double(code(1))/20; % Design lift coefficient
    p = str2double(code(2))/20; % Position of max camber
    t = str2double(code(4:5))/100; % Thickness

    % Camber line coefficients (see Abbott & von Doenhoff)
    if str2double(code(3)) == 0
        % Standard camber
        k1 = 15.957*m;
    else
        % Reflex camber
        k1 = 15.793*m;
    end

    % Camber line
    yc = zeros(size(x));
    dyc_dx = zeros(size(x));
    for i = 1:length(x)
        xi = x(i)/c;
        if xi < p
            yc(i) = (k1/6)*(xi^3 - 3*p*xi^2 + p^2*(3 - p)*xi);
            dyc_dx(i) = (k1/6)*(3*xi^2 - 6*p*xi + p^2*(3-p));
        else
            yc(i) = (k1*p^3/6)*(1 - xi);
            dyc_dx(i) = -(k1*p^3/6);
        end
    end

    theta = atan(dyc_dx);

    % Thickness distribution
    yt = 5*t*c*(0.2969*sqrt(x/c) - 0.1260*(x/c) - 0.3516*(x/c).^2 ...
        + 0.2843*(x/c).^3 - 0.1015*(x/c).^4);

    xu = x - yt.*sin(theta);
    yu = yc + yt.*cos(theta);
    xl = x + yt.*sin(theta);
    yl = yc - yt.*cos(theta);
    xc = x;
end

% --- NACA 6-digit function ---
function [xu, yu, xl, yl, xc, yc] = naca6digit(code, x, c)
    % Only basic 6-series supported: NACA MPXXYY, e.g., 641212
    m = str2double(code(1)); % 6-series
    p = str2double(code(2)); % Position of min pressure
    a = str2double(code(3)); % Design lift coefficient
    t = str2double(code(end-1:end))/100; % Thickness

    % Assume symmetric for now (camber line = 0)
    yc = zeros(size(x));
    dyc_dx = zeros(size(x));

    % Thickness distribution (same as 4/5-series)
    yt = 5*t*c*(0.2969*sqrt(x/c) - 0.1260*(x/c) - 0.3516*(x/c).^2 ...
        + 0.2843*(x/c).^3 - 0.1015*(x/c).^4);

    xu = x - yt;
    yu = yc + yt;
    xl = x + yt;
    yl = yc - yt;
    xc = x;
end
