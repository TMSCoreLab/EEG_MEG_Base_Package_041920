%   This script creates the cortical dipole layer (multiple dipoles)
%
%   Copyright SNM/WAW 2018-2020

%%   Setup extrernal conductivity for dipole sources
dipole_tissue = 'GM';   %   dipoles are always within the gray matter
dipole_tissue_idx = find(strcmp(tissue, dipole_tissue));

%%   Setup dipoles between GM and WM
d  = 0.3e-3;            %   finite-dipole length, m (change to 1e-3)
s  = 2e-3;              %   spacing from GM, m       
R  = 0.007;             %   radius of the enclosing sphere for the cluster in m   
Ctr= 1e-3*[31 0 56];    %   position of enclosing sphere in m
GM              = load('110411_gm.mat');        %   in mm
GM.P            = 1e-3*GM.P;                    %   now in m
GM.Center       = meshtricenter(GM.P, GM.t);    %   base for the dipole layer

%   Indexes into GM triangles strictly within the sphere
indexg1         = find( (GM.P(GM.t(:, 1), 1)-Ctr(1)).^2 + (GM.P(GM.t(:, 1), 2)-Ctr(2)).^2 + (GM.P(GM.t(:, 1), 3)-Ctr(3)).^2 < R^2);
indexg2         = find( (GM.P(GM.t(:, 2), 1)-Ctr(1)).^2 + (GM.P(GM.t(:, 2), 2)-Ctr(2)).^2 + (GM.P(GM.t(:, 2), 3)-Ctr(3)).^2 < R^2);
indexg3         = find( (GM.P(GM.t(:, 3), 1)-Ctr(1)).^2 + (GM.P(GM.t(:, 3), 2)-Ctr(2)).^2 + (GM.P(GM.t(:, 3), 3)-Ctr(3)).^2 < R^2);
indexg          = intersect(intersect(indexg1, indexg2), indexg3); 
M               = length(indexg);

%   Do barycentric subdivision
order = 2;      %   change to 3 if necessary
[coeff, weights, IndexF] = tri(order);
P1 = GM.P(GM.t(indexg, 1), :); %   first triangle node
P2 = GM.P(GM.t(indexg, 2), :); %   second triangle node
P3 = GM.P(GM.t(indexg, 3), :); %   third triangle node
K = length(weights);
strdipolePplus      = [];
strdipolePminus     = [];
for k = 1:K
    Set = coeff(1, k)*P1 + coeff(2, k)*P2 + coeff(3, k)*P3;
    strdipolePplus  = [strdipolePplus; Set-(s - d/2)*GM.normals(indexg, :)];
    strdipolePminus = [strdipolePminus;Set-(s + d/2)*GM.normals(indexg, :)];
end
strdipolesig      = repmat(cond(dipole_tissue_idx), 2*M*K, 1);
clear strdipoleCurrent; 

%%  Define I0 through Okada Murakami constant
AREA = mean(Area(indexg));
I0 = 1e-3*AREA/d/K   
strdipoleCurrent(1:M*K, 1)    = +I0;
strdipoleCurrent(M*K+1:2*M*K, 1)= -I0;

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

%%  Get centerline
Target      = [31 0 56];     %    in mm here
[pointsline] = targetctr('110411_skin.mat', Target);

%%  Plot and check correct position
WM          = load('110411_wm.mat');        %   in mm
WM.P        = 1e-3*WM.P;                    %   now in m
WM.Center   = meshtricenter(WM.P, WM.t);

indexw1 = find( (WM.P(WM.t(:, 1), 1)-Ctr(1)).^2 + (WM.P(WM.t(:, 1), 2)-Ctr(2)).^2 + (WM.P(WM.t(:, 1), 3)-Ctr(3)).^2 < R^2);
indexw2 = find( (WM.P(WM.t(:, 2), 1)-Ctr(1)).^2 + (WM.P(WM.t(:, 2), 2)-Ctr(2)).^2 + (WM.P(WM.t(:, 2), 3)-Ctr(3)).^2 < R^2);
indexw3 = find( (WM.P(WM.t(:, 3), 1)-Ctr(1)).^2 + (WM.P(WM.t(:, 3), 2)-Ctr(2)).^2 + (WM.P(WM.t(:, 3), 3)-Ctr(3)).^2 < R^2);
indexw  = intersect(intersect(indexw1, indexw2), indexw3); 

%% Plot dipole(s) between WM and GM
figure('color', 'w');
str.EdgeColor = 'k'; str.FaceColor = [0 1 1]; str.FaceAlpha = 1.0; 
bemf2_graphics_base(WM.P, WM.t(indexw, :), str);
str.EdgeColor = 'k'; str.FaceColor = [0.5 0.5 0.5]; str.FaceAlpha = 1.0; 
bemf2_graphics_base(GM.P, GM.t(indexg, :), str);
bemf1_graphics_dipole(strdipolePplus, strdipolePminus, strdipoleCurrent, 0);

S = load('sphere');
S.P = 1000*R*S.P;
S.P = S.P + repmat(Ctr, size(S.P, 1), 1);
str.EdgeColor = 'none'; str.FaceColor = [0 1 1]; str.FaceAlpha = 0.05; 
bemf2_graphics_base(S.P, S.t, str);

%%  Plot centerline
hold on;
index = 2000:3000;
plot3(1e-3*pointsline(index, 1), 1e-3*pointsline(index, 2), 1e-3*pointsline(index, 3), '-r', 'lineWidth', 4);

daspect([1 1 1])
camlight; lighting flat;
xlabel('x'); ylabel('y'); zlabel('z');
view(150, -10); axis off; camzoom(1.5);

%% Plot sphere above GM
figure('color', 'w');
tissue_to_plot = 'GM';
t0 = t(Indicator==find(strcmp(tissue, tissue_to_plot)), :);    % (change indicator if necessary: 1-skin, 2-skull, etc.)
str.EdgeColor = 'none'; str.FaceColor = [1 0.75 0.65]; str.FaceAlpha = 1.0; 
bemf2_graphics_base(P, t0, str);
str.EdgeColor = 'none'; str.FaceColor = [0 1 1]; str.FaceAlpha = 0.5; 
bemf2_graphics_base(S.P, S.t, str);

%%  Plot centerline
hold on;
plot3(1e-3*pointsline(:, 1), 1e-3*pointsline(:, 2), 1e-3*pointsline(:, 3), '-r', 'lineWidth', 4);

axis 'equal';  axis 'tight';   
daspect([1 1 1]);
set(gcf,'Color','White');
daspect([1 1 1]);
camlight('headlight'); lighting phong;   
xlabel('x'); ylabel('y'); zlabel('z');
view(33, 53); axis off; camzoom(3)