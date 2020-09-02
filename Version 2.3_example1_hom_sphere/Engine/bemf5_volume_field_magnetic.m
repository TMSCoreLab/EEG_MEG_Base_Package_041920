function B = bemf5_volume_field_magnetic(Points, Pot, P, t, Center, Area, normals, difference, mu0, R, flag, planeABCD)
%   Computes secondary magnetic field from the volumetric current
%   distribution via the FMM in terms of the equivalent surface currents
%   and pseudo electric potential anywhere in space. Performs accurate
%   computations of all neighbor surface integrals
%
%   Copyright SNM 2018-2020

    if(nargin < 12)
        planeABCD = [];
    end

    %  Compute equivalent dipole positions, directions, and moments
    M  = size(Center, 1);  %    number of equivalent dipoles
    
    moments     = -difference.*Pot.*Area;   %   dipole moments
    directions  = normals;                  %   dipole directions
   
    %%  General FMM: input data    
    srcinfo.sources = Center';      %   source points [3, N] - face centers
    srcinfo.nd      = 3;            %   three effective pseudo dipole sets    
    targ            = Points';      %   target points
    prec            = 1e-2;         %   precision-->this value is OK for surfaces
    pg              = 0;            %   nothing is evaluated at sources
    pgt             = 1;            %   potential is evaluated at targets
     
    %%  Three sets (M1, M2,M3) of effective dipole moments
    nx(:, 1) = +0*directions(:, 1).*moments;
    nx(:, 2) = -1*directions(:, 3).*moments;
    nx(:, 3) = +1*directions(:, 2).*moments;
    ny(:, 1) = +1*directions(:, 3).*moments;
    ny(:, 2) = +0*directions(:, 2).*moments;
    ny(:, 3) = -1*directions(:, 1).*moments;
    nz(:, 1) = -1*directions(:, 2).*moments;
    nz(:, 2) = +1*directions(:, 1).*moments;
    nz(:, 3) = +0*directions(:, 3).*moments;
    
    srcinfo.dipoles(1, :, :)     = conj(nx.');                      %   first set of dipole moments  
    srcinfo.dipoles(2, :, :)     = conj(ny.');                      %   second set of dipole moments
    srcinfo.dipoles(3, :, :)     = conj(nz.');                      %   third set of dipole moments
    U                   = lfmm3d(prec, srcinfo, pg, targ, pgt);     %   FMM for three sets
    B                   = mu0/(4*pi)*U.pottarg.';                   %   effective potentials for three sets
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if flag == 0    % Only the center-point approximation is used
        return;
    end   
    %   Undo the effect of the m-th triangle equiv. surface current on neighbor obs. points and
    %   add precise integration instead   
    Size             = mean(sqrt(Area));    %   average triangle size
    const = mu0/(4*pi);     
    if(isempty(planeABCD))
        eligibleTriangles = 1:size(t, 1);
    else
        d1 = abs(planeABCD(1)*Center(:,1) + planeABCD(2)*Center(:,2) + planeABCD(3)*Center(:,3) + planeABCD(4));
        d2 = norm(planeABCD(1:3));
        d = d1./d2;
        eligibleTriangles = find(d <= R*Size);
    end
    ineighborlocal   = rangesearch(Points, Center(eligibleTriangles, :), R*Size, 'NSMethod', 'kdtree'); % over triangles: M by X  
    %ineighborlocal   = rangesearch(Points, Center, R*Size, 'NSMethod', 'kdtree'); % over triangles: M by X  
    %for m =1:M
    for j = 1:length(eligibleTriangles)
        %index       = ineighborlocal{m};  % index into points that are close to triangle m   
        index       = ineighborlocal{j};
        m = eligibleTriangles(j);
        if ~isempty(index)
            temp        = repmat(Center(m, :), length(index), 1) - Points(index, :);   %   these are distances to the observation points
            DIST        = sqrt(dot(temp, temp, 2));                                    %   single column                
            tempn       = repmat(normals(m, :), length(index), 1);
            temp        = cross(tempn, temp, 2);                                       %   cross-product with the normal vector              
            I           = Area(m)*temp./repmat(DIST.^3, 1, 3);                         %   center-point integral, standard format 
            factor1     = +const*difference(m)*Pot(m)*I;
            B(index, :) = B(index, :) - factor1;         
            r1          = P(t(m, 1), :);    %   row
            r2          = P(t(m, 2), :);    %   row
            r3          = P(t(m, 3), :);    %   row           
            I           = potint2(r1, r2, r3, normals(m, :), Points(index, :));     %   analytical precise integration MATLAB
            factor2     = +const*difference(m)*Pot(m)*cross(tempn, I, 2);
            B(index, :)= B(index, :) + factor2;
        end    
    end
end