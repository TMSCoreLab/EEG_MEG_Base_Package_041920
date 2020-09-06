%   This script computes the induced surface charge density for an
%   inhomogeneous multi-tissue object given the primary electric field, with
%   accurate neighbor integration
%
%   Copyright SNM/WAW 2017-2020

%%  Parameters of the iterative solution
iter         = 25;              %    Maximum possible number of iterations in the solution 
relres       = 1e-12;           %    Minimum acceptable relative residual 
weight       = 1/2;             %    Weight of the charge conservation law to be added (empirically found)

%%  Right-hand side b of the matrix equation Zc = b
%   Surface charge density is normalized by eps0: real charge density is eps0*c
tic
%   Gaussian integration is used here
%[Einc, Pinc] = bemf3_inc_field_electric_gauss(strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, P, t);

gaussRadius = 6 * R;
[Einc, Pinc] = bemf3_inc_field_electric_gauss_selective(strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, P, t, Center, Ctr, gaussRadius);
disp([newline 'Incident field calculated in ' num2str(toc) ' s']);

b        = 2*(contrast.*sum(normals.*Einc, 2));                         %  Right-hand side of the matrix equation

%%  GMRES iterative solution (native MATLAB GMRES is used)
h           = waitbar(0.5, 'Please wait - Running MATLAB GMRES'); 
tic
%   MATVEC is the user-defined function of c equal to the left-hand side of the matrix equation LHS(c) = b
MATVEC = @(c) bemf4_surface_field_lhs(c, Center, Area, contrast, normals, weight, EC);     
[c, flag, rres, its, resvec] = gmres(MATVEC, b, [], relres, iter, [], [], b); 
close(h);

%%  Plot convergence history
figure; 
semilogy(resvec/resvec(1), '-o'); grid on;
title('Relative residual of the iterative solution');
xlabel('Iteration number');
ylabel('Relative residual');

%%  Check charge conservation law (optional)
conservation_law_error = sum(c.*Area)/sum(abs(c).*Area)

%%  Check the residual of the integral equation
solution_error = resvec(end)/resvec(1)

%%   Topological low-pass solution filtering (repeat if necessary)
% c = (c.*Area + sum(c(tneighbor).*Area(tneighbor), 2))./(Area + sum(Area(tneighbor), 2));

%%  Save solution data (surface charge density, principal value of surface field)
tic
save('output_charge_solution', 'c', 'resvec', 'conservation_law_error', 'solution_error');

%%   Find and save surface electric potential
Padd = bemf4_surface_field_potential_accurate(c, Center, Area, PC);
%Padd = bemf4_surface_field_potential_subdiv(c, P, t, Area, 'barycentric', 3);
Ptot = Pinc + Padd;     %   Continuous total electric potential at interfaces
save('output_efield_solution.mat', 'Ptot');
