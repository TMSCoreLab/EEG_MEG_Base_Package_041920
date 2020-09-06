%   This script compares analytical and numerical solutions for the potential
%   A finite-length axial dipole within a homogeneous sphere is considered
%
%   Copyright SNM 2018-2020

%   Numerical solution
PotNum = Ptot;  %   Total plotential

tissue_to_plot  = 'Skin';
objectnumber    = find(strcmp(tissue, tissue_to_plot));
index           = Indicator==objectnumber;
PotNum          = PotNum(index);

%   Analytical solution
tic
rad     = 0.092;
if exist('output_analytical_solution.mat', 'file')
    load('output_analytical_solution.mat')
else
    PotAnl  = a_p_1layer_finite(I0, strdipolePplus, strdipolePminus, strdipolesig(1), rad, Center(index, :));
end
disp([newline 'Analytical solution loaded or evaluated in ' num2str(toc) ' s']);

f1 = figure;
bemf2_graphics_surf_field(P, t, PotNum, Indicator, objectnumber);
title('numerical')

f2 = figure;
title('analytical')
bemf2_graphics_surf_field(P, t, PotAnl, Indicator, objectnumber);

Error_2norm = norm(PotNum - PotAnl)/norm(PotAnl)
Error_RDM   = norm(PotNum/norm(PotNum) - PotAnl/norm(PotAnl))