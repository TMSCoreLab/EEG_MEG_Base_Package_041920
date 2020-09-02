function Bpri = bemf3_inc_field_magnetic(strdipolemvector, strdipolemcenter, strdipolemstrength, Points, mu0)  
%   Computes incident magnetic field from the dipole distribution via the FMM
%   in terms of the pseudo electric potential
%
%   Copyright SNM 2018-2020
    
    %  Compute dipole positions, directions, and moments    
    M  = size(strdipolemvector, 1);  % number of subdipoles
    N  = size(Points, 1);             % number of points
    
    %   Define centers, scalar moments and unit directions of magn. dipoles
    segmoment   = sqrt(dot(strdipolemvector, strdipolemvector, 2));
    moments     = segmoment.*strdipolemstrength;
    segpoints   = strdipolemcenter; 
    directions  = strdipolemvector./repmat(segmoment, 1, 3);
    
    %%  General FMM: input data    
    srcinfo.sources = segpoints';   %   source points [3, N] - dipole centers
    srcinfo.nd      = 3;            %   three effective pseudo dipole sets    
    targ            = Points';      %   target points
    prec            = 1e-2;         %   precision-->OK for surfaces
    pg              = 0;            %   nothing is evaluated at sources
    pgt             = 1;            %   only potential is evaluated at targets
     
    %%  Three sets (M1, M2,M3) of effective pseudo dipole moments
    nx(:, 1) = +0*directions(:, 1).*moments;
    nx(:, 2) = -1*directions(:, 3).*moments;
    nx(:, 3) = +1*directions(:, 2).*moments;
    ny(:, 1) = +1*directions(:, 3).*moments;
    ny(:, 2) = +0*directions(:, 2).*moments;
    ny(:, 3) = -1*directions(:, 1).*moments;
    nz(:, 1) = -1*directions(:, 2).*moments;
    nz(:, 2) = +1*directions(:, 1).*moments;
    nz(:, 3) = +0*directions(:, 3).*moments;
    
    srcinfo.dipoles(1, :, :)     = conj(nx.');  %   first set of dipole moments  
    srcinfo.dipoles(2, :, :)     = conj(ny.');  %   second set of dipole moments
    srcinfo.dipoles(3, :, :)     = conj(nz.');  %   third set of dipole moments
    U = lfmm3d(prec, srcinfo, pg, targ, pgt);   %   FMM for three sets
    Bpri = mu0/(4*pi)*U.pottarg.';              %   potentials for three sets
end

