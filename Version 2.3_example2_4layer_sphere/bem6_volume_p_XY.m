%   This script accurately computes and displays electric potential sampled
%   on a cross-section (transverse plane) via the FMM method with accurate
%   neighbor integration
%
%   Copyright SNM/WAW 2017-2020

%%  Load/prepare data
planeABCD = [0 0 1 -Z*1e-3]; %Analytical equation of the observation point plane for neighbor search acceleration
label       = 't';

%% Define observation points in the cross-section
Ms = 200;
x = linspace(xmin, xmax, Ms);
y = linspace(ymin, ymax, Ms);
[X0, Y0]  = meshgrid(x, y);
clear pointsXY;
pointsXY(:, 1) = reshape(X0, 1, Ms^2);
pointsXY(:, 2) = reshape(Y0, 1, Ms^2);  
pointsXY(:, 3) = Z*ones(1, Ms^2);

%% Find the potential at each observation point
tic
pointsXY       = 1e-3*pointsXY;     % Convert back to m
Ppri           = zeros(Ms*Ms, 1);
Psec           = zeros(Ms*Ms, 1);
[~, Ppri]      = bemf3_inc_field_electric(strdipolePplus,...
                                        strdipolePminus,...
                                        strdipolesig, strdipoleCurrent, ...                                        
                                        P, t, pointsXY, Area, normals, 0, 0); 
R = 2;  %   precise integration
Psec           = bemf5_volume_field_potential(pointsXY, c, P, t, Center, Area, normals, R, planeABCD);
Ptotal         = Ppri + Psec;   
disp([newline 'Potential calculated in ' num2str(toc) ' s']);

%% Plot the potential in the cross-section
figure;
% Potential contour plot
temp      = Ptotal;
th1 = +1e-7;             %   in V
th2 = -1e-7;             %   in V
levels      = 10;
bemf2_graphics_vol_field(temp, th1, th2, levels, x, y);
xlabel('Distance x, mm');
ylabel('Distance y, mm');
title(strcat('Potential (V), ', label, '-in the transverse plane'));

% Tissue boundaries
color   = prism(length(tissue)); color(3, :) = [1 0 1];
for m = countXY
    edges           = EofXY{m};             %   this is for the contour
    points          = [];
    points(:, 1)    = +PofXY{m}(:, 1);       %   this is for the contour  
    points(:, 2)    = +PofXY{m}(:, 2);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end

% Dipoles
hold on
bemf1_graphics_dipole(1e3*strdipolePplus, 1e3*strdipolePminus, 1e3*strdipoleCurrent, 3);

% General settings 
axis 'equal';  axis 'tight';     
colormap parula; colorbar;
axis([xmin xmax ymin ymax]);
grid on; set(gcf,'Color','White');