function B = a_m_infinite(p, r0, R)
%   Computes magnetic field outside for a current dipole inside
%   Inputs:
%   p - vector current dipole moment (A times m)
%   r0 - dipole position, 1x3, m
%   R - Mx3 array of positions outside the sphere where the field 
%   is plotted
%   Output(s):
%   B - magnetic field outside the sphere
%   Source: Sarvas J. Basic mathematical and electromagnetic concepts of
%   the biomagnetic inverse problem. Phys. Med. Biol. 1987; 32(1):11-22. 
%
%   SNM 2018-2020

    mu0         = 1.25663706e-006;
    M = size(R, 1); B = zeros(M, 3);
    for m = 1:M
        r = R(m, :);  %   Single observation point
        a   = r  - r0;
        as  = sqrt(dot(a, a));
        rs  = sqrt(dot(r, r));
        F   = as*(rs*as + dot(r, a));
        G2F = (as^2/rs + 2*as + 2*rs + dot(r, a)/as)*r - ...
              (as + 2*rs + dot(r, a)/as)*r0;        
        B(m, :) = +mu0/(4*pi*F^2)*(F*cross(p, r0) ...
                  - dot(cross(p, r0), r)*G2F);       
    end
end

