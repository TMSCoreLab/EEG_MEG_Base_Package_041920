function Pot = a_p_1layer_finite(I0, Pplus, Pminus, sigma, rad, R)
%   Computes surface electric potential for an axial finite-length
%   electric dipole within a homogeneous sphere with conductivity sigma
%   Inputs:
%   I0     - current amplitude
%   PPlus  - position of +I0 source
%   PMinus - position of -I0 source
%   sigma  - conductivity
%   rad    - radius of the sphere
%   R      - Mx3 array of positions on the sphere where the potential is plotted
%   Output(s):
%   Pot - surface electric potential
%   Source: Ernest Frank, Electric Potential Produced by Two Point Current
%   Sources in a Homogeneous Conducting Sphere, Journal of Applied Physics
%   23, 1225 (1952); doi: 10.1063/1.1702037
%
%   SNM 2018-2020

    M = size(R, 1); Pot = zeros(M, 1);
    N = 50;
    a = Pminus(3);
    b = Pplus(3);  
    parpool(16);
    parfor m = 1:M
        r = R(m, :);  %   Single observation point
        costheta = r(3)/rad;
        Pot(m) = 0;
        LP     = legendreP(1:N, costheta);
        for n = 1:N
            Pot(m) = Pot(m) + ...
                (((2*n+1)/n)/rad^(n+1))*(b^n - a^n)*LP(n); 
        end
        Pot(m) = (I0/(4*pi*sigma))*Pot(m);                 
    end
    delete(gcp('nocreate'));
end

