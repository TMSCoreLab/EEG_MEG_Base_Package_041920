%   This script compares analytical and numerical solutions for the
%   magnetic field. An infinitesimally-short dipole within a four-layer
%   sphere (or any multilayer sphere) is considered
%
%   Copyright SNM 2018-2020

obs     = load('meshsphereo_97.mat');       %   observation surface, mm
Points  = 1e-3*meshtricenter(obs.P, obs.t); %   observation surface, m

%   Numerical solution
difference     = condin - condout;
Bpri           = bemf3_inc_field_magnetic(strdipolemvector, strdipolemcenter, strdipolemstrength, Points, mu0);                                     
R = 5;         %   precise integration            
Bsec           = bemf5_volume_field_magnetic(Points, Ptot, P, t, Center, Area, normals, difference, mu0, R, 0);
Btotal         = Bpri + Bsec;     
Bnum           = Btotal;
temp           = sqrt(dot(Bnum, Bnum, 2));

f1 = figure;
bemf2_graphics_surf_field(obs.P, obs.t, temp, ones(1, size(obs.t, 1)), 1);
title('numerical')

%   Analytical solution
dvector     =      strdipolePplus - strdipolePminus;
dcenter     = 0.5*(strdipolePplus + strdipolePminus);
moment      = repmat(strdipoleCurrent(1:end/2), 1, 3).*dvector;
Banl        = zeros(size(Bnum));
for m = 1:size(strdipolePplus, 1)
    Banl        = Banl + a_m_infinite(moment(m, :), dcenter(m, :), Points);
end
temp        = sqrt(dot(Banl, Banl, 2));

f2 = figure;
bemf2_graphics_surf_field(obs.P, obs.t, temp, ones(1, size(obs.t, 1)), 1);
title('analytical')

Error_2norm = norm(Bnum - Banl)/norm(Banl)
Error_RDM   = norm(Bnum/norm(Bnum) - Banl/norm(Banl))

