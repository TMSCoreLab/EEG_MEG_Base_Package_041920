function [Epri, Ppri] = bemf3_inc_field_electric(strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, ...                                        
                                         P, t, Center, Area, normals, R, flag)
%   Computes potential and electric field from the dipole distribution via the FMM
%   Includes accurate neighbor triangle integrals for facets located close
%   to the dipoles
%
%   Copyright SNM 2018-2020

    %   Define source (pole) positions and FMM pseudo charges         
    Positions   = [strdipolePplus; strdipolePminus];
    PseudoQ     = strdipoleCurrent./strdipolesig;
    %   FMM 2019
    srcinfo.nd      = 1;                    %   one vector of charges  
    srcinfo.sources = Positions';           %   source points
    targ            = Center';              %   target points
    prec            = 1e-2;                 %   precision->OK for surfaces    
    pg      = 0;                            %   nothing is evaluated at sources
    pgt     = 2;                            %   potential/field are evaluated at targets
    srcinfo.charges(1, :)    = PseudoQ.';   %   pseudo charges    
    U                        = lfmm3d(prec, srcinfo, pg, targ, pgt);
    Ppri                     = +1/(4*pi)*U.pottarg.';
    Epri(:, 1)               = -1/(4*pi)*U.gradtarg(1, :);
    Epri(:, 2)               = -1/(4*pi)*U.gradtarg(2, :);
    Epri(:, 3)               = -1/(4*pi)*U.gradtarg(3, :); 
    
    if flag == 0    % Only the center-point approximation is used 
        return;
    end    
    %   Replace the center-point approximation by precise integration when
    %   triangles are close to the dipole sources. For every triangle,
    %   variable ineighborlocal returns index into the closest source
    %   positions    
    Size             = mean(sqrt(Area));    %   average triangle size
    PositionsCenter  = (Positions(1:end/2, :) + Positions(end/2+1:end, :))/2;   %   This is dipole center position 
                                                                                %   Critical for the loop; otherwise one pole may be ignored 
    N = size(PositionsCenter, 1);                                               %   Total number of dipoles
    ineighborlocal   = rangesearch(PositionsCenter, Center, R*Size, 'NSMethod', 'kdtree'); 
    % Loop over triangles: M by X  
    M = size(Center, 1);
    CurrentOverSigma = strdipoleCurrent./strdipolesig;
    for m = 1:M
        inde       = ineighborlocal{m};               %   index into dipole centers that are close to triangle m   
        if ~isempty(inde)
            index        = [inde N+inde];             %   index into poles that are close to triangle m  
            VectorCurrent   = repmat(CurrentOverSigma(index), 1, 3); 
            temp            = Positions(index, :) - repmat(Center(m, :), length(index), 1);  %   these are distances to the observation point
            DIST            = sqrt(dot(temp, temp, 2));                                      %   single column   
            I               = VectorCurrent.*temp./repmat(DIST.^3, 1, 3);                    %   integral, standard format        
            Epri(m, :)      = Epri(m, :) - (-1/(4*pi)*sum(I, 1));                            %   for the m-th triangle, undo the effect of all neighbor sources  
            r1 = P(t(m, 1), :);
            r2 = P(t(m, 2), :);
            r3 = P(t(m, 3), :);       
            Int = potint2(r1, r2, r3, normals(m, :), Positions(index, :));
            Int = VectorCurrent./(4*pi).*Int/Area(m);
            Epri(m, :) = Epri(m, :) + sum(Int, 1);
            I              =  CurrentOverSigma(index)./DIST;                        %   integral, standard format        
            Ppri(m)        = Ppri(m) - (1/(4*pi)*sum(I, 1));                        %   for the m-th triangle, undo the effect of all neighbor sources          
            [Int, ~]       = potint(r1, r2, r3, normals(m, :), Positions(index, :));
            Int            = CurrentOverSigma(index)./(4*pi).*Int/Area(m);      
            Ppri(m) = Ppri(m) + sum(Int, 1);
        end    
    end      
end