function [Epri, Ppri] = bemf3_inc_field_electric_gauss(strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, P, t)
%   Computes potential and electric field from the dipole distribution via the FMM
%   on traingles using Gaussian quadrature(s)
%
%   Copyright SNM 2018-2020

    M = size(t, 1);
    %   Gaussian quadrature of 5th order with seven integration points 
    N  = 7;
    a1          =   0.797426985353087;
    b1          =   0.101286507323456;       
    a2          =   0.059715871789770;
    b2          =   0.470142064105115;	
    coeff(:, 1)  = [1/3 1/3 1/3]';
    coeff(:, 2)  = [a1  b1  b1]';
    coeff(:, 3)  = [b1  a1  b1]';
    coeff(:, 4)  = [b1  b1  a1]';       
    coeff(:, 5)  = [a2  b2  b2]';
    coeff(:, 6)  = [b2  a2  b2]';
    coeff(:, 7)  = [b2  b2  a2]';    
    weights(1)    =   0.2250000;
    weights(2)    =   0.1259392;
    weights(3)    =   0.1259392;
    weights(4)    =   0.1259392;    
    weights(5)    =   0.1323942;
    weights(6)    =   0.1323942;
    weights(7)    =   0.1323942;
    
    P1 = P(t(:, 1), :); %   first triangle node
    P2 = P(t(:, 2), :); %   second triangle node
    P3 = P(t(:, 3), :); %   third triangle node
    
    Set1 = coeff(1, 1)*P1 + coeff(2, 1)*P2 + coeff(3, 1)*P3;
    Set2 = coeff(1, 2)*P1 + coeff(2, 2)*P2 + coeff(3, 2)*P3;
    Set3 = coeff(1, 3)*P1 + coeff(2, 3)*P2 + coeff(3, 3)*P3;
    Set4 = coeff(1, 4)*P1 + coeff(2, 4)*P2 + coeff(3, 4)*P3;
    Set5 = coeff(1, 5)*P1 + coeff(2, 5)*P2 + coeff(3, 5)*P3;
    Set6 = coeff(1, 6)*P1 + coeff(2, 6)*P2 + coeff(3, 6)*P3;
    Set7 = coeff(1, 7)*P1 + coeff(2, 7)*P2 + coeff(3, 7)*P3;
    Set  = [Set1; Set2; Set3; Set4; Set5; Set6; Set7];
    
    [Eall, Pall] = bemf3_inc_field_electric_plain...
    (strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, Set);
    
    Ppri = zeros(M, 1);
    Epri = zeros(M, 3);
    for m = 1:N 
        index = 1+(m-1)*M:m*M;
        Ppri    = Ppri + weights(m)*Pall(index);
        Epri    = Epri + weights(m)*Eall(index, :);
    end
end