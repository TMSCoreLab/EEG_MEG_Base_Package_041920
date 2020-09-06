function [Epri, Ppri] = bemf3_inc_field_electric_gauss_selective(strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, P, t, Center, dipoleClusterCenter, gaussRadius)    
    %% Compute model subdivision parameters
                        %   number of integration points in the Gaussian quadrature  
                        %   for the outer potential integrals
                        %   Numbers 1, 4, 7, 13, 25 are permitted 
    %   Gaussian weights for analytical integration (for the outer integral)
    gauss = 7;
    if      gauss == 0  
        [coeffS, weightsS, IndexS]  = tri(subdivParam);
    elseif  gauss == 1  
        [coeffS, weightsS, IndexS]  = tri(1, 1);
    elseif  gauss == 4  
        [coeffS, weightsS, IndexS]  = tri(4, 3);
    elseif  gauss == 7
        [coeffS, weightsS, IndexS]  = tri(7, 5);
    elseif  gauss == 13 
        [coeffS, weightsS, IndexS]  = tri(13, 7);
    elseif  gauss == 25 
        [coeffS, weightsS, IndexS]  = tri(25, 10);
    else
        error('Invalid Gaussian subdivision parameter');
    end
    
    %% Find all triangles that require Gaussian subdivision based on distance from center of dipole cluster
    % Find indices of triangles that need to be subdivided
    
    % This commented-out method may be useful for multiple dipole clusters.
    % (as written, it only extracts the points close to the first cluster)
    %trianglesToSubdiv_temp = rangesearch(Center, dipoleClusterCenter, gaussRadius);
    % Convert to logical array
    %trianglesToSubdiv = logical(zeros(size(t, 1), 1));
    %trianglesToSubdiv(trianglesToSubdiv_temp{1}) = true;
    
    % This method is much faster for a single dipole cluster
    tempDist = vecnorm(Center - dipoleClusterCenter, 2, 2);
    trianglesToSubdiv = tempDist <= gaussRadius;
    
    
    %% Preallocate output variables
    Epri = zeros(size(t, 1), 3);
    Ppri = zeros(size(t, 1), 1);
    
    %% Calculate incident fields on triangles that do not require subdivision
    [Epri(~trianglesToSubdiv, :), Ppri(~trianglesToSubdiv, :)] = bemf3_inc_field_electric_plain(strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, Center(~trianglesToSubdiv, :));
    
    %% Subdivide the triangles that do require subdivision
    Center_subdiv = zeros(IndexS*sum(trianglesToSubdiv), 3);
    P1 = P(t(trianglesToSubdiv,1), :);
    P2 = P(t(trianglesToSubdiv,2), :);
    P3 = P(t(trianglesToSubdiv,3), :);
    for j = 1:IndexS
        currentIndices = ([1:sum(trianglesToSubdiv)] - 1)*IndexS + j; 
        Center_subdiv(currentIndices, :) = coeffS(1, j)*P1 + coeffS(2, j)*P2 + coeffS(3, j)*P3;
    end
    
    %% Calculate incident electric fields on subdivided triangles
    [E_subdiv, P_subdiv] = bemf3_inc_field_electric_plain(strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, Center_subdiv);
    
    %% Recover average electric field at whole triangles from subdivided triangles    
    %Every column contains the subdivided quantities for one full triangle
    P_subdiv_temp = reshape(P_subdiv, IndexS, []);
    Ex_temp = reshape(E_subdiv(:,1), IndexS, []);
    Ey_temp = reshape(E_subdiv(:,2), IndexS, []);
    Ez_temp = reshape(E_subdiv(:,3), IndexS, []);

    %% Write to output variables
    Ppri(trianglesToSubdiv) = transpose(weightsS*P_subdiv_temp);
    Epri(trianglesToSubdiv, 1) = transpose(weightsS*Ex_temp);
    Epri(trianglesToSubdiv, 2) = transpose(weightsS*Ey_temp);
    Epri(trianglesToSubdiv, 3) = transpose(weightsS*Ez_temp);
