%==========================================================================
%
% legendre_recursion  Recursive evaluation of coefficients related to the 
% Legendre polynomials required for evaluating the gravitational 
% acceleration.
%
%   [V,W] = legendre_recursion(r_ecef,R,N)
%
% Author: Tamas Kis
% Last Update: 2022-02-15
%
% REFERENCES:
%   [1] Montenbruck and Gill, "Satellite Orbits" (pp. 66-68)
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   r_ecef  - (3×1 double) position resolved in ECEF frame [m]
%   R       - (1×1 double) Earth mean equatorial radius [m]
%   N       - (1×1 double) maximum order/degree of gravity model
%
% -------
% OUTPUT:
% -------
%   V       - (1×1 function_handle) function to be used as "V(n,m)" to
%             return Legendre coefficient of degree n and order m
%   W       - (1×1 function_handle) function to be used as "W(n,m)" to
%             return Legendre coefficient of degree n and order m
%
%==========================================================================
function [V,W] = legendre_recursion(r_ecef,R,N)
    
    % magnitude of r_ecef squared
    r_sqr = idot(r_ecef,r_ecef);

    % TODO
    r_frac = R*R/r_sqr;
    
    % extracts position components
    x = R*r_ecef(1)/r_sqr;
    y = R*r_ecef(2)/r_sqr;
    z = R*r_ecef(3)/r_sqr;

    % distance from center of the Earth [m]
    r = inorm(r_ecef);

    % initialize V and W
    V = zeros(N+2,N+2);
    W = zeros(N+2,N+2);

    % first zonal term (n = m = 0)
    V(1,1) = R/r;
    
    % second zonal terms (n = 1, m = 0)
    V(2,1) = z*V(1,1);

    % remaining zonal terms (2 ≤ n ≤ N+1, m = 0)
    for n = 2:(N+1)

        % degree (n) for indexing
        ni = n+1;
        
        % coefficients
        V(ni,1) = ((2*n-1)*z*V(ni-1,1)-(n-1)*r_frac*V(ni-2,1))/n;

    end

    % tesseral and sectorial terms
    for m = 1:(N+1)
        
        % order (m) for indexing
        mi = m+1;
        
        % sectorial terms (n = m)
        V(mi,mi) = (2*m-1)*(x*V(mi-1,mi-1)-y*W(mi-1,mi-1));
        W(mi,mi) = (2*m-1)*(x*W(mi-1,mi-1)+y*V(mi-1,mi-1));

        % tesseral terms (n > m)
        for n = (m+1):(N+1)

            % degree (n) for indexing
            ni = n+1;
            
            % coefficients
            V(ni,mi) = ((2*n-1)*z*V(ni-1,mi)-(ni+mi-1)*r_frac*V(ni-2,...
                mi))/(n-m);
            W(ni,mi) = ((2*n-1)*z*W(ni-1,mi)-(ni+mi-1)*r_frac*W(ni-2,...
                mi))/(n-m);

        end

    end

    % assign function handles to return coefficients of degree n + order m
    V = @(n,m) V(n+1,m+1);
    W = @(n,m) W(n+1,m+1);

end