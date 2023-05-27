function [F,U] = AltStatePotential(N,n1,n2,nhat,k,c,v)
% This is a function to make the potential for the polynomial growth
% function that has alternative stable states separated by an unstable
% equilibrium. 

% The program takes the following input
%   N = the vector of population densities over which to calculate the
%       potential
%   n1 and n2 = the lower two values of N that correspond to the 
%       alternative stable equilibria. 
%   nhat = the value of N that is the unstable equilibrium. Note that nhat
%       must be between n1 and n2.
%   k = an integer value (including 0) that scales the order of the
%       polynomial.
%   c = a positive real value that scales the speed of dynamics and the
%      overall steepness of the potential across all N
%   v = the value of the environmental forcing term at time t.

%% 
    %Need to put in error statements here
%%
 
tot_terms = 1+2*k+1;
Uterms = zeros(tot_terms,length(N));
Fterms = zeros(tot_terms,length(N));

for i = 0:(1+2*k)
    bincoef = nchoosek(1+2*k,i)*(-nhat)^(1+2*k-i);
    Fterms(i+1,:) = bincoef*N.^(3+i)...
                    -(n1+n2)*bincoef*N.^(2+i)...
                    +(n1*n2)*bincoef*N.^(1+i);
                
    Uterms(i+1,:) = bincoef*N.^(4+i)/(4+i)...
                    -(n1+n2)*bincoef*N.^(3+i)/(3+i)...
                    +n1*n2*bincoef*N.^(2+i)/(2+i);
end


F = -c*ones(1,tot_terms)*Fterms + v*N;
U = c*ones(1,tot_terms)*Uterms - v*N.^2/2;
end

