%   This script creates the base dipole (single)
%
%   Copyright SNM/WAW 2018-2020

%%   Single dipole example
I0 = 1e-6;              %   source current, A

% Vertical short dipole (0.04 mm length)
% strdipolePplus      = [0.000 0.000 0.07552];    %   in m
% strdipolePminus     = [0.000 0.000 0.07548];    %   in m

%  Horizontal short dipole (0.04 mm length)
strdipolePplus      = [+0.00002 0.000 0.0755];  %   in m
strdipolePminus     = [-0.00002 0.000 0.0755];  %   in m

%   Rotate dipole if necessary
theta = 0.0;
strdipolePplus  = meshrotate2(strdipolePplus, [1 0 0], theta);
strdipolePminus = meshrotate2(strdipolePminus, [1 0 0], theta);

Ctr = mean([strdipolePplus; strdipolePminus]);
R = 0.002;
dlength = norm(strdipolePplus - strdipolePminus)

strdipolesig           = [cond(4) cond(4)]';
strdipoleCurrent       = [+I0 -I0]';
M                      = size(strdipolePplus, 1);

%%   Magnetic dipole subdivision (optional)
D = 10;                        %   number of smaller subdipoles
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
R           = 0.01;                         %   in m     
WM          = load('meshsphere6_73.mat');   %   in mm
WM.P        = 1e-3*WM.P;                    %   now in m
WM.Center   = meshtricenter(WM.P, WM.t);

GM          = load('meshsphere6_78.mat');   %   in mm
GM.P        = 1e-3*GM.P;                    %   now in m
GM.Center   = meshtricenter(GM.P, GM.t);

Ctr    = mean((strdipolePplus + strdipolePminus)/2, 1);

indexg1 = find( (GM.P(GM.t(:, 1), 1)-Ctr(1)).^2 + (GM.P(GM.t(:, 1), 2)-Ctr(2)).^2 + (GM.P(GM.t(:, 1), 3)-Ctr(3)).^2 < R^2);
indexg2 = find( (GM.P(GM.t(:, 2), 1)-Ctr(1)).^2 + (GM.P(GM.t(:, 2), 2)-Ctr(2)).^2 + (GM.P(GM.t(:, 2), 3)-Ctr(3)).^2 < R^2);
indexg3 = find( (GM.P(GM.t(:, 3), 1)-Ctr(1)).^2 + (GM.P(GM.t(:, 3), 2)-Ctr(2)).^2 + (GM.P(GM.t(:, 3), 3)-Ctr(3)).^2 < R^2);
indexg  = intersect(intersect(indexg1, indexg2), indexg3); 

indexw1 = find( (WM.P(WM.t(:, 1), 1)-Ctr(1)).^2 + (WM.P(WM.t(:, 1), 2)-Ctr(2)).^2 + (WM.P(WM.t(:, 1), 3)-Ctr(3)).^2 < R^2);
indexw2 = find( (WM.P(WM.t(:, 2), 1)-Ctr(1)).^2 + (WM.P(WM.t(:, 2), 2)-Ctr(2)).^2 + (WM.P(WM.t(:, 2), 3)-Ctr(3)).^2 < R^2);
indexw3 = find( (WM.P(WM.t(:, 3), 1)-Ctr(1)).^2 + (WM.P(WM.t(:, 3), 2)-Ctr(2)).^2 + (WM.P(WM.t(:, 3), 3)-Ctr(3)).^2 < R^2);
indexw  = intersect(intersect(indexw1, indexw2), indexw3); 

% Plot dipole between WM and GM
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

% Plot dipole above WM
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