function [Epri, Ppri] = bemf3_inc_field_electric_plain(strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, Points)
%   Computes potential and electric field from the dipole distribution via the FMM
%   at observation points (Points)
%
%   Copyright SNM 2018-2020

    %   Define source (pole) positions and FMM pseudo charges         
    Positions   = [strdipolePplus; strdipolePminus];
    PseudoQ     = strdipoleCurrent./strdipolesig;
    %   FMM 2019
    srcinfo.nd      = 1;                    %   one vector of charges  
    srcinfo.sources = Positions';           %   source points
    targ            = Points';              %   target points
    prec            = 1e-2;                 %   precision->OK for surfaces    
    pg      = 0;                            %   nothing is evaluated at sources
    pgt     = 2;                            %   potential/field are evaluated at targets
    srcinfo.charges(1, :)    = PseudoQ.';   %   pseudo charges    
    U                        = lfmm3d(prec, srcinfo, pg, targ, pgt);
    Ppri                     = +1/(4*pi)*U.pottarg.';
    Epri(:, 1)               = -1/(4*pi)*U.gradtarg(1, :);
    Epri(:, 2)               = -1/(4*pi)*U.gradtarg(2, :);
    Epri(:, 3)               = -1/(4*pi)*U.gradtarg(3, :); 
end