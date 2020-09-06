%   This script accurately computes and displays magnetic fields sampled on
%   a cross-section (sagittal plane) via the FMM method with accurate
%   neighbor integration
%
%   Copyright SNM/WAW 2018-2020

%%  Load/prepare data
planeABCD = [1 0 0 -X*1e-3]; % Analytical equation of observation point plane for neighbor search acceleration

%%  Post processing parameters
component   = 4;        %   field component to be plotted (1, 2, 3 or x, y, z, or 4 - total) 
temp        = ['x' 'y' 'z' 't'];
label       = temp(component);

%%  Define observation points in the cross-section
Ms = 200;
y = linspace(ymin, ymax, Ms);
z = linspace(zmin, zmax, Ms);
[Y0, Z0]  = meshgrid(y, z);
clear pointsYZ;
pointsYZ(:, 1) = X*ones(1, Ms^2);
pointsYZ(:, 2) = reshape(Y0, 1, Ms^2);
pointsYZ(:, 3) = reshape(Z0, 1, Ms^2);

%% Find the B-field at each observation point in the cross-section
tic
pointsYZ       = 1e-3*pointsYZ;     % Convert back to m
Bpri           = zeros(Ms*Ms, 3);
Bsec           = zeros(Ms*Ms, 3);
difference     = condin - condout;
Bpri           = bemf3_inc_field_magnetic(strdipolemvector, strdipolemcenter, strdipolemstrength, pointsYZ, mu0);                                     
R = 5;        %   precise integration            
Bsec           = bemf5_volume_field_magnetic(pointsYZ, Ptot, P, t, Center, Area, normals, difference, mu0, R, 1, planeABCD);
Btotal         = Bpri + Bsec;     
disp([newline 'B-field calculated in ' num2str(toc) ' s']);

%%  Plot the B-field in the cross-section
figure;
% B-field contour plot
if component == 4
    temp      = abs(sqrt(dot(Btotal, Btotal, 2)));
else
    temp      = abs(Btotal(:, component));
end
th1 = 1e-13;                %   in T
th2 = 0;                    %   in T
levels      = 10;
bemf2_graphics_vol_field(temp, th1, th2, levels, y, z);
xlabel('Distance y, mm');
ylabel('Distance z, mm');
title(strcat('B-field (T), ', label, '-component in the sagittal plane'));

% Tissue boundaries
color   = prism(length(tissue)); color(4, :) = [0 1 1];
for m = countYZ
    edges           = EofYZ{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofYZ{m}(:, 2);       %   this is for the contour  
    points(:, 2)    = +PofYZ{m}(:, 3);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end

% Dipole
hold on
bemf1_graphics_dipole(1e3*strdipolePplus, 1e3*strdipolePminus, strdipoleCurrent, 2);

% General settings 
axis 'equal';  axis 'tight';     
colormap parula; colorbar;
axis([ymin ymax zmin zmax]);
grid on; set(gcf,'Color','White');