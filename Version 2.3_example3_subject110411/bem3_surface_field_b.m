%   This script computes and plots the incident magnetic field magnitude
%   (or any of the components) for any brain compartment surface/interface
%   (plots the surface field + optionally coil geometry) 
%
%   Copyright SNM/WAW 2017-2020

%%   Compute the B-field for all surfaces/interfaces
tissue_to_plot = 'Skin';
objectnumber = find(strcmp(tissue, tissue_to_plot));    
Points = Center(Indicator==objectnumber, :); 
Normals = normals(Indicator==objectnumber, :); 

d = 10e-3; % Magnetometer distance from skin surface
obsPtsMag = Points + d*Normals; % Observation points for magnetic field

tic
h    = waitbar(0.5, 'Please wait - computing magnetic field');
difference      = condin - condout;
Bpri            = bemf3_inc_field_magnetic(strdipolemvector, strdipolemcenter, strdipolemstrength, obsPtsMag, mu0);   
Bsec            = bemf5_volume_field_magnetic(obsPtsMag, Ptot, P, t, Center, Area, normals, difference, mu0, 0, 0);
Btotal          = Bpri + Bsec;
temp            = abs(sqrt(dot(Btotal, Btotal, 2)));
close(h);
BTime = toc

%%  Digitize figure
figure;
step = 10;
temp = round(step*temp/max(temp)).*(max(temp))/step;
bemf2_graphics_surf_field(P, t, temp, Indicator, objectnumber);
title(strcat('Solution: Magnetic field in T at 10 mm for:', tissue{objectnumber}));

%%  Plot centerline
hold on;
plot3(1e-3*pointsline(:, 1), 1e-3*pointsline(:, 2), 1e-3*pointsline(:, 3), '-r', 'lineWidth', 3);
view(148, 38); axis off; camzoom(1)