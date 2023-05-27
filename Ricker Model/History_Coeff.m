function coef = History_Coeff(abar,Omega)

% This function calculates the coefficient c found in eqn. (5) of the main
% text and derived in the supplementary material. 

% This function takes 
%   abar <- the average of the environmental driver
%   Omega <- the period of the cycle of the environmental driver.

% Exogenous model
gen = 50000;

t = 0:gen;
adiff = sqrt(2)*sin(2*pi*t/Omega);


% Autocorrelation structure

uvec = 0:1000; % Time lags
tsteps = gen/2;
rho = zeros(1,length(uvec)); %Empty vector to hold computed autocorrelation

% loop over lags and calculate the autcorrelation
for i = 1:length(uvec)
    u = uvec(i);
    t_indx = (uvec(end)+1):(tsteps+uvec(end)+1);
    tlag_indx = (uvec(end)+1- u):(uvec(end)+1 - u + tsteps);
    rho(i) = mean(adiff(tlag_indx).*adiff(t_indx));
end

%% Calculate terms in approximation
lead_coef = (1-abar).^uvec(2:end).*rho(2:end);
coef = 2*sum(lead_coef);

end

