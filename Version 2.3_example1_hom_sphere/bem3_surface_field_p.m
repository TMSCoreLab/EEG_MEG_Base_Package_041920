%   This script plots the electric potential of the primary, secondary, or
%   the full field for any brain compartment surface
%
%   Copyright SNM/WAW 2017-2020

%%   Graphics
tissue_to_plot = 'Skin';
objectnumber    = find(strcmp(tissue, tissue_to_plot));
temp            = Ptot(Indicator==objectnumber);

figure;
bemf2_graphics_surf_field(P, t, temp, Indicator, objectnumber);
title(['Solution: Electric potential in V for: ', tissue{objectnumber}]);