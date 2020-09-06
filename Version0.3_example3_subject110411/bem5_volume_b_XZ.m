%   This script accurately computes and displays magnetic fields sampled on
%   a surface (coronal plane) via the FMM method with accurate
%   neighbor triangle integration. 
%
%   Copyright SNM 2018-2020

%%  Load/prepare data
planeABCD = [0 1 0 -Y*1e-3]; %Analytical equation of the observation plane for neighbor triangle search acceleration

%  Post processing parameters
component   = 4;        %   field component to be plotted (1, 2, 3 or x, y, z, or 4 - total) 
temp        = ['x' 'y' 'z' 't'];
label       = temp(component);

%%  Define observation points in the cross-section (MsxMs observation points)  
Ms = 200;
x = linspace(xmin, xmax, Ms);
z = linspace(zmin, zmax, Ms);
[X0, Z0]  = meshgrid(x, z);
clear pointsXZ;
pointsXZ(:, 1) = reshape(X0, 1, Ms^2);
pointsXZ(:, 2) = Y*ones(1, Ms^2);
pointsXZ(:, 3) = reshape(Z0, 1, Ms^2);  

%%  Find the B-field at each observation point in the cross-section
tic
pointsXZ       = 1e-3*pointsXZ;     % Convert back to m
Bpri           = zeros(Ms*Ms, 3);
Bsec           = zeros(Ms*Ms, 3);
difference     = condin - condout;
Bpri           = bemf3_inc_field_magnetic(strdipolemvector, strdipolemcenter, strdipolemstrength, pointsXZ, mu0);                                     
R = 5;        %   precise integration            
Bsec           = bemf5_volume_field_magnetic(pointsXZ, Ptot, P, t, Center, Area, normals, difference, mu0, R, 1, planeABCD);
Btotal         = Bpri + Bsec;     
fieldPlaneTime = toc  

%%  Plot the B-field in the cross-section
figure;
% Contour plot
if component == 4
    temp      = abs(sqrt(dot(Btotal, Btotal, 2)));
else
    temp      = abs(Btotal(:, component));
end
th1 = +1e-11;         %   in T
th2 = 0;              %   in T
levels      = 10;
bemf2_graphics_vol_field(temp, th1, th2, levels, x, z);
xlabel('Distance x, mm');
ylabel('Distance z, mm');
title(strcat('B-field (T), ', label, '-component in the coronal plane'));
 
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
bemf1_graphics_dipole(1e3*strdipolePplus, 1e3*strdipolePminus, strdipoleCurrent, 2);

%   General settings 
axis 'equal';  axis 'tight';     
colormap parula; colorbar;
axis([xmin xmax zmin zmax]);
grid on; set(gcf,'Color','White');