function Potential = bemf4_surface_field_potential_accurate(c, Center, Area, PC)
%   This function computes CONTINUOUS electric field and the full potential
%   on a surface facet due to charges on ALL OTHER facets including
%   accurate neighbor integrals. Self-terms causing discontinuity may not
%   be included for electric field
%   To obtain the true field/potential, divide the result(s) by eps0;

%   Copyright SNM 2018-2020    
   
    %  FMM 2019   
    %   Potentials of surface charges
    %   FMM plus correction
    tic
    const           = 1/(4*pi);
    eps             = 1e-3;
    pg              = 1;            %   only potential here
    srcinfo.sources = Center';
    srcinfo.charges = (c.*Area)';
    U               = lfmm3d(eps, srcinfo, pg);
    Potential       = +U.pot';    
    %   Near-field correction   
    Potential = const*Potential;
    Potential = Potential + PC*c;    
    toc     
end
