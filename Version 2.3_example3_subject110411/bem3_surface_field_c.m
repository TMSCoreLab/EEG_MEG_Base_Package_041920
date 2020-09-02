%   This script plots the induced surface charge density (proportional to
%   the normal electric field) for any brain compartment surface (plots the
%   density)
%
%   Copyright SNM/WAW 2018-2020

%%   Graphics
tissue_to_plot = 'CSF';
objectnumber    = find(strcmp(tissue, tissue_to_plot));
temp            = eps0*c(Indicator==objectnumber);  % the real charge density is eps0*c

%%  Digitize figure
figure;
step = 5;
temp = round(step*temp/max(temp)).*(max(temp))/step;
bemf2_graphics_surf_field(P, t, temp, Indicator, objectnumber);
title(strcat('Solution: Surface charge density in C/m^2 for:', tissue{objectnumber}));

%%  Plot centerline
hold on;
plot3(1e-3*pointsline(:, 1), 1e-3*pointsline(:, 2), 1e-3*pointsline(:, 3), '-r', 'lineWidth', 3);
view(148, 38); axis off; camzoom(1)