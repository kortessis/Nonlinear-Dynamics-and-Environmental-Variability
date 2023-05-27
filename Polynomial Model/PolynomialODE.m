function dNdt = PolynomialODE(t,N,stable_eq,c,k,delta,omega)
%ODE object for the polynomial growth model
%   This function requires the parameter values for the stable equilibria,
%  of which there are two. The unstable equilibrium is placed exactly
%  between these two. The growth rate parameter c (positive real number) 
%  and parameter k (integer value, 0 inclusive) are also inputs.

% c scales how fast the population grows across the entire function.
% Dynamics are faster with larger c. 

% k controls the steepness of the ridge separating the two stable
% equilibria. The ridge is steepest when k = 0 and becomes flatter as k
% increases.

unstable_eq = mean(stable_eq);
F = -c.*N.*(N-stable_eq(1)).*(N-stable_eq(2)).*(N-unstable_eq).^(1+2*k);

Forcing = delta*cos(omega*t);

dNdt = F + N.*Forcing;

end

