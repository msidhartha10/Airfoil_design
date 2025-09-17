% NACA 5- and 6-Series Airfoil Plot Generator with Cosine Spacing
% and Excel export. Supports Open or Closed Trailing Edge.
%
% Notes:
% - 5-digit camber line implementation is an approximate formula tuned for
%   the common p=0.15 case (e.g., NACA 23xxx). It will work best for those.
%   Reflexed camber is supported approximately.
% - 6-series cambered mean line is not included (complex pressure distribution
%   method). 6-series is generated as symmetric (zero camber).
% - Thickness distribution uses classic NACA 4-digit polynomial with a
%   selectable trailing edge: open (-0.1015) or closed (-0.1036).
%
% By: Copilot
% MATLAB R2019b+ recommended (uses writetable). For older versions, switch to xlswrite.

clear; clc;

% ------------ User Parameters ------------
% Choose series and code (examples):
%   - 5-series: '23012'  (approx cambered 5-digit, best for p=0.15 i.e. second digit=3)
%   - 6-series: '641212' (treated as symmetric thickness-only in this script)
naca_series = 5;                 % 5 or 6
naca_code   = '23012';           % e.g., '23012' or '641212'
chord       = 1.0;               % chord length
num_points  = 201;               % number of cosine-spaced points (>= 2)
output_dir  = '.';               % output folder for Excel
prompt_for_trailing_edge = true; % set false to force trailing edge type via variable below
trailing_edge_type = 'open';     % 'open' or 'closed' (used if prompt_for_trailing_edge=false)
% ----------------------------------------

% Ask user for trailing edge type (as requested)
if prompt_for_trailing_edge
    te_in = input('Trailing edge type (open/closed) [open]: ', 's');
    te_in = strtrim(lower(te_in));
    if isempty(te_in)
        trailing_edge_type = 'open';
    elseif any(strcmp(te_in, {'open','closed'}))
        trailing_edge_type = te_in;
    else
        warning('Invalid trailing edge type. Defaulting to OPEN.');
        trailing_edge_type = 'open';
    end
else
    % validate provided trailing_edge_type
    if ~any(strcmpi(trailing_edge_type, {'open','closed'}))
        warning('Invalid trailing_edge_type provided. Defaulting to OPEN.');
        trailing_edge_type = 'open';
    end
end
fprintf('Using a %s trailing edge.\n', upper(trailing_edge_type));

% Cosine spacing along chord (0 at LE to c at TE)
beta = linspace(0, pi, num_points);
x = (1 - cos(beta)) * (chord/2);

% Generate coordinates
switch naca_series
    case 5
        [xu, yu, xl, yl, xc, yc] = airfoil_5digit(naca_code, x, chord, trailing_edge_type);
    case 6
        [xu, yu, xl, yl, xc, yc] = airfoil_6series(naca_code, x, chord, trailing_edge_type);
    otherwise
        error('Only NACA 5 or 6 series are supported.');
end

% Plot
figure('Color','w');
plot(xu, yu, 'b', 'LineWidth', 1.5); hold on;
plot(xl, yl, 'r', 'LineWidth', 1.5);
plot(xc, yc, 'k--', 'LineWidth', 1.2);
axis equal; grid on; box on;
xlabel('x'); ylabel('y');
title(sprintf('NACA %s (%s TE) - %s points', naca_code, trailing_edge_type, num2str(num_points)));
legend('Upper Surface','Lower Surface','Camber Line','Location','best');

% Prepare and export to Excel
out_table = table(xu(:), yu(:), xl(:), yl(:), xc(:), yc(:), ...
    'VariableNames', {'Upper_X','Upper_Y','Lower_X','Lower_Y','Camber_X','Camber_Y'});

outfile = fullfile(output_dir, sprintf('NACA_%s_%s_TE.xlsx', naca_code, trailing_edge_type));
try
    writetable(out_table, outfile);
    fprintf('Exported coordinates to: %s\n', outfile);
catch ME
    warning('Failed to write Excel file via writetable. Error: %s', ME.message);
    % Fallback: write CSVs
    out_csv = fullfile(output_dir, sprintf('NACA_%s_%s_TE.csv', naca_code, trailing_edge_type));
    writetable(out_table, out_csv);
    fprintf('Exported coordinates to CSV (fallback): %s\n', out_csv);
end

% --------------- Functions ---------------

function [xu, yu, xl, yl, xc, yc] = airfoil_5digit(code, x, c, te_type)
    % Parse 5-digit code LPQTT
    % L (1st digit): design lift Cl = 0.15 * L
    % P (2nd digit): pos of max camber at p = P/20
    % Q (3rd digit): 0 = standard camber, 1 = reflex
    % TT (last two): thickness in % chord
    if length(code) ~= 5 || any(~isstrprop(code, 'digit'))
        error('Invalid NACA 5-digit code. Example: 23012');
    end
    L = str2double(code(1));
    P = str2double(code(2));
    Q = str2double(code(3));
    TT = str2double(code(4:5));

    Cl_design = 0.15 * L;
    p = P / 20;      % location of max camber
    t = TT / 100;    % thickness

    % Approximate 5-digit mean camber line (works best for p ~= 0.15)
    % This uses the common k1 constants scaled by Cl_design/0.3 at p=0.15,
    % and slightly adjusted for reflex. For other p values this is an
    % approximation.
    [yc, dyc_dx] = meanline_5digit_approx(x/c, p, Cl_design, Q);

    theta = atan(dyc_dx);
    yt = thickness_distribution(x, c, t, te_type);

    xu = x - yt .* sin(theta);
    yu = yc .* c + yt .* cos(theta);
    xl = x + yt .* sin(theta);
    yl = yc .* c - yt .* cos(theta);

    xc = x;
    % yc currently non-dimensional; scale camber line to chord
    yc = yc .* c;
end

function [xu, yu, xl, yl, xc, yc] = airfoil_6series(code, x, c, te_type)
    % Basic 6-series: thickness only (symmetric), camber not modeled here.
    % Accepts codes like '641212' but ignores camber-specific content.
    if any(~isstrprop(code, 'digit'))
        error('Invalid NACA 6-series code. Example: 641212');
    end
    t = str2double(code(end-1:end)) / 100;

    yc = zeros(size(x));         % symmetric mean line
    dyc_dx = zeros(size(x));     % zero slope
    theta = zeros(size(x));      % zero

    yt = thickness_distribution(x, c, t, te_type);

    xu = x - yt .* sin(theta);
    yu = yc + yt .* cos(theta);
    xl = x + yt .* sin(theta);
    yl = yc - yt .* cos(theta);

    xc = x;
    % yc is already dimensional (zeros)
end

function yt = thickness_distribution(x, c, t, te_type)
    % NACA thickness distribution with selectable trailing edge
    % te_type: 'open' => a5 = -0.1015 (non-zero TE thickness)
    %          'closed' => a5 = -0.1036 (zero TE thickness)
    a0 = 0.2969;
    a1 = -0.1260;
    a2 = -0.3516;
    a3 = 0.2843;

    if strcmpi(te_type, 'closed')
        a4 = -0.1036; % closed trailing edge
    else
        a4 = -0.1015; % open trailing edge (classic NACA)
    end

    xc = x ./ c;
    yt = 5 * t * c .* (a0*sqrt(max(xc,0)) + a1*xc + a2*xc.^2 + a3*xc.^3 + a4*xc.^4);
end

function [yc, dyc_dx] = meanline_5digit_approx(xc, p, Cl_design, reflexFlag)
    % Approximate NACA 5-digit mean line.
    % Uses a piecewise cubic with k1 scaled around the common p=0.15 case.
    %
    % Inputs:
    %   xc         - nondimensional x (x/c)
    %   p          - location of max camber (P/20)
    %   Cl_design  - design lift coefficient (0.15*L)
    %   reflexFlag - 0 = standard, 1 = reflex
    %
    % Output:
    %   yc, dyc_dx (both nondimensional)

    % Choose a baseline K1 for p=0.15 case, scale with Cl.
    % For standard (non-reflex), many references list K1_base ≈ 15.957 for Cl=0.3.
    % For reflex, K1_base ≈ 15.793 for Cl=0.3.
    Cl_ref = 0.30; % reference design Cl for the base constants
    if reflexFlag == 0
        K1_base = 15.957;
    else
        K1_base = 15.793;
    end

    % Heuristic scaling for other p: adjust by (0.15/p)^2 to keep shape reasonable
    % This is a practical approximation; for high fidelity use tabulated 5-digit constants.
    p = max(min(p, 0.35), 0.05); % clamp p to a reasonable range
    scale_p = (0.15 / p)^2;

    k1 = K1_base * (Cl_design / Cl_ref) * scale_p;

    % Piecewise camber line (Abbott & von Doenhoff style)
    yc = zeros(size(xc));
    dyc_dx = zeros(size(xc));
    for i = 1:numel(xc)
        xci = xc(i);
        if xci < p
            yc(i) = (k1/6) * (xci^3 - 3*p*xci^2 + p^2*(3 - p)*xci);
            dyc_dx(i) = (k1/6) * (3*xci^2 - 6*p*xci + p^2*(3-p));
        else
            if reflexFlag == 0
                yc(i) = (k1 * p^3 / 6) * (1 - xci);
                dyc_dx(i) = -(k1 * p^3 / 6);
            else
                % Simple reflex tail-up tweak: subtract a small cubic near TE
                % This approximates reflex behavior (non-exact).
                yc(i) = (k1 * p^3 / 6) * (1 - xci) - 0.01 * (xci - p).^3;
                dyc_dx(i) = -(k1 * p^3 / 6) - 0.03 * (xci - p).^2;
            end
        end
    end
end
