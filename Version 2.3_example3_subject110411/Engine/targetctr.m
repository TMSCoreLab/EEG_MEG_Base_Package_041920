function [pointsline] = targetctr(name, Target)
%   Outputs centerline for Target vs skin in mm

    %%   Load shell
    load(name);
    center = meshtricenter(P, t);

    %%  Find nearest skin faces and construct the electrode position
    Rnumber     = 8;                          %   number of neighbor triangles 
    ineighborS  = knnsearch(center, Target, 'k', Rnumber);
    tempc       = center(ineighborS(1), :);
    tempn       = mean(normals(ineighborS, :), 1);
    tempn       = tempn/(tempn(3));             %   normalization
    Nx          = tempn(1);
    Ny          = tempn(2);
    Nz          = tempn(3);
    dir         = tempn/norm(tempn);
    tempc       = tempc + 10*dir;            %   10 mm away from skin
    MoveX       = tempc(1);
    MoveY       = tempc(2);
    MoveZ       = tempc(3);
    M = 10000;        
    argline      = linspace(0, 100, M);                 %   distance along a 100 mm long line   
    dirline      = -[Nx Ny Nz]/norm([Nx Ny Nz]);        %   line direction (inside)   
    offline      = 0e-3;                                %   offset from the electrode
    pointsline(1:M, 1) = MoveX + dirline(1)*(argline + offline);
    pointsline(1:M, 2) = MoveY + dirline(2)*(argline + offline);
    pointsline(1:M, 3) = MoveZ + dirline(3)*(argline + offline);
end

