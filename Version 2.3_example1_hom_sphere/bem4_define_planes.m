% Prepare cross-section observation planes for subsequent visualizations

%% Parameters

X    = 0;
Y    = 0;
Z    = 75;
xmin = -100;
xmax = +100;
ymin = -100;
ymax = +100;
zmin = -100;
zmax = +100;



%% Model initialization and preparation
tic
%   Preallocate cell arrays
m_max = length(tissue);
tS = cell(m_max, 1);
nS = tS; %  Reuse this empty cell array for other initialization
eS = tS;
TriPS = tS;
TriMS = tS;
Eplot = tS;
PS = P * 1e3; % Convert to mm
for m = 1:m_max
    tS{m} = t(Indicator == m, :);
    nS{m} = normals(Indicator == m, :);
    [eS{m}, TriPS{m}, TriMS{m}] = mt(tS{m}); 
    %Eplot{m}  = Eoutside_m(Indicator == m);      %   Select the field here
end
disp([newline 'Initial surface data processed in ' num2str(toc) ' s']);
%SurfaceDataProcessTime = toc

%%  Plot cross-sectional plane locations with respect to the head model
figure;
% Plot the planes
patch(1e-3*[xmin xmin xmax xmax],1e-3*[ymin ymax ymax ymin], 1e-3*[Z Z Z Z], 'c', 'FaceAlpha', 0.35);
patch(1e-3*[xmin xmin xmax xmax],1e-3*[Y Y Y Y], 1e-3*[zmin zmax zmax zmin], 'c', 'FaceAlpha', 0.35);
patch(1e-3*[X X X X], 1e-3*[ymin ymin ymax ymax], 1e-3*[zmin zmax zmax zmin], 'c', 'FaceAlpha', 0.35);

% Plot the selected tissue of the head model
tissue_to_plot = 'Skin';
t0 = t(Indicator==find(strcmp(tissue, tissue_to_plot)), :);    % (change indicator if necessary: 1-skin, 2-skull, etc.)
str.EdgeColor = 'none'; str.FaceColor = [1 0.75 0.65]; str.FaceAlpha = 1.0; 
bemf2_graphics_base(P, t0, str);

% General settings
axis 'equal';  axis 'tight';   
daspect([1 1 1]);
set(gcf,'Color','White');
camlight; lighting phong;
view(157, 25);

%% Compute contours of intersected tissues in the XY cross-section
tissues = length(tissue);
PofXY = cell(tissues, 1);   %   intersection nodes for a tissue
EofXY = cell(tissues, 1);   %   edges formed by intersection nodes for a tissue
TofXY = cell(tissues, 1);   %   intersected triangles
NofXY = cell(tissues, 1);   %   normal vectors of intersected triangles
countXY = [];               %   number of every tissue present in the slice
for m = 1:tissues 
    [Pi, ti, polymask, flag] = meshplaneintXY(PS, tS{m}, eS{m}, TriPS{m}, TriMS{m}, Z);
    if flag % intersection found                
        countXY               = [countXY m];
        PofXY{m}            = Pi;               %   intersection nodes
        EofXY{m}            = polymask;         %   edges formed by intersection nodes
        TofXY{m}            = ti;               %   intersected triangles
        NofXY{m}            = nS{m}(ti, :);     %   normal vectors of intersected triangles        
    end
end


%% Compute contours of intersected tissues in the XZ cross-section
tissues = length(tissue);
PofXZ = cell(tissues, 1);   %   intersection nodes for a tissue
EofXZ = cell(tissues, 1);   %   edges formed by intersection nodes for a tissue
TofXZ = cell(tissues, 1);   %   intersected triangles
NofXZ = cell(tissues, 1);   %   normal vectors of intersected triangles
countXZ = [];               %   number of every tissue present in the slice
for m = 1:tissues 
    [Pi, ti, polymask, flag] = meshplaneintXZ(PS, tS{m}, eS{m}, TriPS{m}, TriMS{m}, Y);
    if flag % intersection found                
        countXZ               = [countXZ m];
        PofXZ{m}            = Pi;               %   intersection nodes
        EofXZ{m}            = polymask;         %   edges formed by intersection nodes
        TofXZ{m}            = ti;               %   intersected triangles
        NofXZ{m}            = nS{m}(ti, :);     %   normal vectors of intersected triangles        
    end
end

%% Compute contours of intersected tissues in the YZ cross-section
tissues = length(tissue);
PofYZ = cell(tissues, 1);   %   intersection nodes for a tissue
EofYZ = cell(tissues, 1);   %   edges formed by intersection nodes for a tissue
TofYZ = cell(tissues, 1);   %   intersected triangles
NofYZ = cell(tissues, 1);   %   normal vectors of intersected triangles
countYZ = [];               %   number of every tissue present in the slice
for m = 1:tissues 
    [Pi, ti, polymask, flag] = meshplaneintYZ(PS, tS{m}, eS{m}, TriPS{m}, TriMS{m}, X);
    if flag % intersection found                
        countYZ               = [countYZ m];
        PofYZ{m}            = Pi;               %   intersection nodes
        EofYZ{m}            = polymask;         %   edges formed by intersection nodes
        TofYZ{m}            = ti;               %   intersected triangles
        NofYZ{m}            = nS{m}(ti, :);     %   normal vectors of intersected triangles        
    end
end

%% Plot tissue contours in the XY cross-section
figure;

% Display the contours
color   = prism(length(tissue)); color(4, :) = [0 1 1];
for m = countXY
    edges           = EofXY{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofXY{m}(:, 1);       %   this is for the contour  
    points(:, 2)    = +PofXY{m}(:, 2);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end
patch([xmin xmin xmax xmax],[ymin ymax ymax ymin], 'c', 'FaceAlpha', 0.35);
% General settings
title(['Cross-section in the transverse plane at Z = ' num2str(Z) ' mm']);
xlabel('x, mm'); ylabel('y, mm');
axis 'equal';  axis 'tight'; 
set(gcf,'Color','White');

%% Plot tissue contours in the XZ cross-section
figure;

% Display the contours
color   = prism(length(tissue)); color(4, :) = [0 1 1];
for m = countXZ
    edges           = EofXZ{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofXZ{m}(:, 1);       %   this is for the contour  
    points(:, 2)    = +PofXZ{m}(:, 3);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end
patch([xmin xmin xmax xmax], [zmin zmax zmax zmin], 'c', 'FaceAlpha', 0.35);
% General settings
title(['Cross-section in the coronal plane at Y = ' num2str(Y) ' mm']);
xlabel('x, mm'); ylabel('z, mm');
axis 'equal';  axis 'tight'; 
set(gcf,'Color','White');

%% Plot tissue contours in the YZ cross-section
figure;

% Display the contours
color   = prism(length(tissue)); color(4, :) = [0 1 1];
for m = countYZ
    edges           = EofYZ{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofYZ{m}(:, 2);       %   this is for the contour  
    points(:, 2)    = +PofYZ{m}(:, 3);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end
patch([ymin ymin ymax ymax], [zmin zmax zmax zmin], 'c', 'FaceAlpha', 0.35);
% General settings
title(['Cross-section in the sagittal plane at X = ' num2str(X) ' mm']);
xlabel('y, mm'); ylabel('z, mm');
axis 'equal';  axis 'tight'; 
set(gcf,'Color','White');
