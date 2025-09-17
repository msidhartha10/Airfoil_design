clc;
clear vars;
close all;

c = 180;                           % Chord length in mm

% 299 total points, so 149 upper, 149 lower, 1 at LE (shared)
[xu, yu, xl, yl, x, yc] = generate_naca4_airfoil('5420',c, 299, true); %[output:6e8f5d5a] %[output:4b36cc68]

%%
function [xu, yu, xl, yl, x, yc] = generate_naca4_airfoil(code,c, nTotalPoints, closedTE)

%% Inputs:
% code         = '2412' (NACA 4-digit airfoil)
% nTotalPoints = total number of surface points (e.g., 299)
% closedTE     = true (closed trailing edge) or false (open)

%% Parameters from NACA code
m = str2double(code(1)) / 100;     % Max camber
p = str2double(code(2)) / 10;      % Max camber location
t = str2double(code(3:4)) / 100;   % Thickness


%% Cosine-spaced x values (0 to 1 chord, cosine distribution)
n = nTotalPoints;
beta = linspace(0, pi, n);
x = (1 - cos(beta)) / 2;

%% Thickness distribution
a0 = 0.2969;
a1 = -0.1260;
a2 = -0.3516;
a3 = 0.2843;
a4 = closedTE * (-0.1036) + (~closedTE) * (-0.1015);  % closed/open TE

yt = 5 * t * (a0 * sqrt(x) + a1 * x + a2 * x.^2 + a3 * x.^3 + a4 * x.^4);

%% Camber line and slope
yc = zeros(size(x));
dyc_dx = zeros(size(x));

for i = 1:length(x)
    if x(i) < p
        yc(i) = m/p^2 * (2*p*x(i) - x(i)^2);
        dyc_dx(i) = 2*m/p^2 * (p - x(i));
    else
        yc(i) = m/(1 - p)^2 * ((1 - 2*p) + 2*p*x(i) - x(i)^2);
        dyc_dx(i) = 2*m/(1 - p)^2 * (p - x(i));
    end
end

theta = atan(dyc_dx);

%% Upper and lower surface
xu = x - yt .* sin(theta);
yu = yc + yt .* cos(theta);
xl = x + yt .* sin(theta);
yl = yc - yt .* cos(theta);

%% Scale by chord length (180 mm)
xu = xu * c;
xl = xl * c;
x = x * c;
yu = yu * c;
yl = yl * c;
yc = yc * c;

%% Plot
figure;
hold on;
axis equal;
grid on;
plot(xu, yu, 'r-', 'LineWidth', 0.8);
plot(xl, yl, 'b-', 'LineWidth', 0.8);
plot(x, yc, 'k--');
legend('Upper Surface','Lower Surface','Camber Line');
title(['NACA ', code, ' Airfoil (Chord = ', num2str(c), ' mm)']);

%% Export to Excel
T = array2table([xu', yu', xl', yl', x', yc'], ...
    'VariableNames', {'X_upper', 'Y_upper', 'X_lower', 'Y_lower', 'X_Chamberline', 'Y_chamberline'});

filename = ['NACA', code, '_airfoil_', num2str(c), 'mm.xlsx'];
writetable(T, filename);
disp(['Excel file "', filename, '" has been successfully exported.']);

end
