%   This script accurately computes and displays electric potential sampled
%   on a cross-section (coronal plane) via the FMM method with accurate
%   neighbor integration
%
%   Copyright SNM/WAW 2018-2020

planeABCD = [0 1 0 -Y*1e-3]; % Analytical equation of observation point plane for neighbor search acceleration
label       = 't';

%%  Define observation points in the cross-section
Ms = 200;
x = linspace(xmin, xmax, Ms);
z = linspace(zmin, zmax, Ms);
[X0, Z0]  = meshgrid(x, z);
clear pointsXZ;
pointsXZ(:, 1) = reshape(X0, 1, Ms^2);
pointsXZ(:, 2) = Y*ones(1, Ms^2);
pointsXZ(:, 3) = reshape(Z0, 1, Ms^2);  

%% Find potential at each observation point in the cross-section       
tic
pointsXZ       = 1e-3*pointsXZ;     % Convert back to m
Ppri           = zeros(Ms*Ms, 1);
Psec           = zeros(Ms*Ms, 1);
[~, Ppri]      = bemf3_inc_field_electric(strdipolePplus,...
                                        strdipolePminus,...
                                        strdipolesig, strdipoleCurrent, ...                                        
                                        P, t, pointsXZ, Area, normals, 0, 0);   
R = 4;  %   precise integration            
Psec           = bemf5_volume_field_potential(pointsXZ, c, P, t, Center, Area, normals, R, planeABCD);
Ptotal         = Ppri + Psec;   
disp([newline 'Potential calculated in ' num2str(toc) ' s']);

%% Plot the potential in the cross-section
figure;
% Potential contour plot
temp      = Ptotal;
th1 = +3*1e-6;           %   in V/m
th2 = -3*1e-6;           %   in V/m
levels      = 20;
bemf2_graphics_vol_field(temp, th1, th2, levels, x, z);
xlabel('Distance x, mm');
ylabel('Distance z, mm');
title(['Potential (V), ', label, '-in the coronal plane']);
 
% Tissue boundaries
color   = prism(length(tissue)); color(4, :) = [0 1 1];
for m = countXZ
    edges           = EofXZ{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofXZ{m}(:, 1);       %   this is for the contour  
    points(:, 2)    = +PofXZ{m}(:, 3);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end

% Dipoles
hold on
bemf1_graphics_dipole(1e3*strdipolePplus, 1e3*strdipolePminus, strdipoleCurrent, 1);

% General settings 
axis 'equal';  axis 'tight';     
colormap parula; colorbar;
axis([xmin xmax zmin zmax]);
grid on; set(gcf,'Color','White');