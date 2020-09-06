%   This script creates the cortical dipole layer (multiple dipoles)
%
%   Copyright SNM/WAW 2018-2020

%%   Multiple dipole example
I0 = 1e-6;              %   source current, A
d  = 1e-3;              %   finite-dipole length, m
s  = 2.5e-3;            %   spacing from GM, m       
R  = 0.005;             %   radius of the enclosing sphere in m   
Ctr= [0 0 75.5e-3];     %   position of enclosing sphere in m
GM              = load('meshsphere6_78.mat');   %   in mm
GM.P            = 1e-3*GM.P;                    %   now in m
GM.Center       = meshtricenter(GM.P, GM.t);    %   base for the dipole layer
%   Indexes into GM triangles strictly within the sphere
indexg1         = find( (GM.P(GM.t(:, 1), 1)-Ctr(1)).^2 + (GM.P(GM.t(:, 1), 2)-Ctr(2)).^2 + (GM.P(GM.t(:, 1), 3)-Ctr(3)).^2 < R^2);
indexg2         = find( (GM.P(GM.t(:, 2), 1)-Ctr(1)).^2 + (GM.P(GM.t(:, 2), 2)-Ctr(2)).^2 + (GM.P(GM.t(:, 2), 3)-Ctr(3)).^2 < R^2);
indexg3         = find( (GM.P(GM.t(:, 3), 1)-Ctr(1)).^2 + (GM.P(GM.t(:, 3), 2)-Ctr(2)).^2 + (GM.P(GM.t(:, 3), 3)-Ctr(3)).^2 < R^2);
indexg          = intersect(intersect(indexg1, indexg2), indexg3); 
M               = length(indexg);
strdipolePplus              = GM.Center(indexg, :)  - (s + 0)*GM.normals(indexg, :);
strdipolePminus             = GM.Center(indexg, :)  - (s + d)*GM.normals(indexg, :);
strdipolesig                = repmat(cond(4), 2*M, 1);
clear strdipoleCurrent; 
strdipoleCurrent(1:M, 1)    = +I0;
strdipoleCurrent(M+1:2*M, 1)= -I0;

NoDipoles = size(strdipolePplus, 1)

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

%%  Plot and check correct position
WM          = load('meshsphere6_73.mat');    %   in mm
WM.P        = 1e-3*WM.P;                    %   now in m
WM.Center   = meshtricenter(WM.P, WM.t);

indexw1 = find( (WM.P(WM.t(:, 1), 1)-Ctr(1)).^2 + (WM.P(WM.t(:, 1), 2)-Ctr(2)).^2 + (WM.P(WM.t(:, 1), 3)-Ctr(3)).^2 < R^2);
indexw2 = find( (WM.P(WM.t(:, 2), 1)-Ctr(1)).^2 + (WM.P(WM.t(:, 2), 2)-Ctr(2)).^2 + (WM.P(WM.t(:, 2), 3)-Ctr(3)).^2 < R^2);
indexw3 = find( (WM.P(WM.t(:, 3), 1)-Ctr(1)).^2 + (WM.P(WM.t(:, 3), 2)-Ctr(2)).^2 + (WM.P(WM.t(:, 3), 3)-Ctr(3)).^2 < R^2);
indexw  = intersect(intersect(indexw1, indexw2), indexw3); 

% Plot dipole(s) between WM and GM
f1 = figure;
str.EdgeColor = 'k'; str.FaceColor = [0 1 1]; str.FaceAlpha = 1.0; 
bemf2_graphics_base(WM.P, WM.t(indexw, :), str);
str.EdgeColor = 'k'; str.FaceColor = [0.5 0.5 0.5]; str.FaceAlpha = 1.0; 
bemf2_graphics_base(GM.P, GM.t(indexg, :), str);
bemf1_graphics_dipole(strdipolePplus, strdipolePminus, strdipoleCurrent, 4) 
axis 'equal';  axis 'tight';   
daspect([1 1 1]);
set(gcf,'Color','White');
camlight; lighting phong;
xlabel('x'); ylabel('y'); zlabel('z');
view(0, 0);

% Plot dipole(s) above WM
f2 = figure;
tissue_to_plot = 'WM';
t0 = t(Indicator==find(strcmp(tissue, tissue_to_plot)), :);    % (change indicator if necessary: 1-skin, 2-skull, etc.)
str.EdgeColor = 'none'; str.FaceColor = [0 1 1]; str.FaceAlpha = 1.0; 
bemf2_graphics_base(P, t0, str);
bemf1_graphics_dipole(strdipolePplus, strdipolePminus, strdipoleCurrent, 4) 
axis 'equal';  axis 'tight';   
daspect([1 1 1]);
set(gcf,'Color','White');
camlight; lighting phong;
xlabel('x'); ylabel('y'); zlabel('z');
view(0, 0);