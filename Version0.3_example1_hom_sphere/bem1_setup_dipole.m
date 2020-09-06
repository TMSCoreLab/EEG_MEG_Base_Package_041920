%   This script creates the base dipole (single)
%
%   Copyright SNM/WAW 2018-2020

%%   Single-dipole example
I0 = 1e-6;              %   source current, A

%  Vertical finite-length dipole (2 mm length)
strdipolePplus      = [0.000 0.000 0.076];    %   in m
strdipolePminus     = [0.000 0.000 0.074];    %   in m 

Ctr = mean([strdipolePplus; strdipolePminus]);
R = 0.002;
dlength = norm(strdipolePplus - strdipolePminus)

strdipolesig           = [cond(1) cond(1)]';
strdipoleCurrent       = [+I0 -I0]';
M                      = size(strdipolePplus, 1);

%%   Magnetic dipole subdivision (optional)
D = 1;                        %   number of smaller subdipoles
strdipolemvector   = zeros(D*M, 3);
strdipolemcenter   = zeros(D*M, 3);
strdipolemstrength = zeros(D*M, 1);
for m = 1:M
    temp = (1/D)*(strdipolePplus(m, :) - strdipolePminus(m, :));
    for d = 1:D 
        arg = d+D*(m-1);
        strdipolemvector(arg, :)     = temp;
        strdipolemcenter(arg, :)     = strdipolePminus(m, :) + (d-1/2)*temp;
        strdipolemstrength(arg, :)   = strdipoleCurrent(m);                  
    end
end