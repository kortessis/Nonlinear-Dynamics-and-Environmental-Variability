clear
clc

%%% Code to evlauate the variance in the Ricker model with exogenous and 
%%% stochastic variation using both simulations and the approximation.

%%% This code looks at different values of abar and different levels of
%%% stochastic variability. 

%%% This code requires the function "History_Coeff", which calculates the
%%% parameter c in equation (5) of the main text. 
%% Global parameters

% Choose different values of abar
abar_vec = linspace(0.25,1.75,3);

% Plotting colors from the "viridis" colormap.
colors = viridis(length(abar_vec)+2);

% Time points in the simulation
gen = 100000;

% Values of stochastic variation to run simulations
sigma = sqrt(linspace(0,0.25,20));

% Generate random variability for reproducibility
rng(10)
X = randn([1,gen]);

% Burnin time after which we evaluate variance in simulations. 
tinit = 2000;

% Set the average population size and initial population size
Nbar = 1;
Ninit = 0.95*Nbar;

%% Exogenous Cycle

% Set the level of deviations in the exogenous cycle
sigma_a = sqrt(0.015);
% Set the period of the cycle
cycle_length = 20;

hold on
% Loop over each value of the average low-density growth rate
for i = 1:length(abar_vec)
    
    abar = abar_vec(i);
    % Define the time sequence of a_t values
    at = abar + sqrt(2*sigma_a^2)*sin(2*pi*[1:gen]/cycle_length);
    
    % Set b (Nbar = abar/b --> b = abar/Nbar)
    b = abar/Nbar;
    
    % Create empty vector to hold population sizes
    N = zeros(length(sigma),gen);
    % Initiate population vector
    N(:,1) = Ninit;
    
    % Simulate dynamics
    for t = 2:gen
        N(:,t) = N(:,t-1).*exp(at(t-1) - b*N(:,t-1) + sigma'*X(t-1));
    end
    
    % Calculate the coefficient representing history of the cycle and the
    % effect of autocorrelations
    coef = History_Coeff(abar, cycle_length);
    
    % Approximation of the variance
    varNpredict =  Nbar^2/(1 - (1-abar)^2)*(sigma_a^2*(1+coef) + sigma.^2);
    
    % Calculation of the variance from simulation
    varNactual = var(N',1);

    % Plot results
    p = plot(sigma.^2, varNpredict, ':');
    p.LineWidth = 3; p.Color = colors(i+1,:);
    p = plot(sigma.^2, varNactual);
    p.LineWidth = 3; p.Color = colors(i+1,:);
    txt = text(max(sigma.^2)+0.0025,varNpredict(end),[' $\bar{a} = $ ', ...
        num2str(abar), '; Approx'],'Interpreter', 'Latex',...
        'HorizontalAlignment', 'left', 'FontSize', 20);
    txt = text(max(sigma.^2)+0.0025,varNactual(end),[' $\bar{a} = $ ', ...
        num2str(abar), '; Sim'],'Interpreter', 'Latex',...
        'HorizontalAlignment', 'left', 'FontSize', 20); 
end
hold off

% Plot formatting
xlabel('\sigma^2');
ylabel('Var({\itN})');
title({'Variance Scaling with','Exogenous and Stochastic Variation'});
ax = gca; ax.FontSize = 25; ax.FontName = 'Times New Roman';
xlim([0,max(sigma.^2)])
