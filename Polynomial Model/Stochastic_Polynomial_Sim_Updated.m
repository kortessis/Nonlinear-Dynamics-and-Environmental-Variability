%% The sde version of the polynomial model

clear; clc;

%% Basic model parameters
% Parameters for the model
n1 = 10; n2 = 12; stable_eq = [n1,n2];
unstable_eq = mean(stable_eq);
deltan = max(n1,n2) - unstable_eq;
gamma = 0.001; delta = 0.005; % Magnitude of the seasonal forcing term.

% Parameters related to the seasonal forcing term
t_int = 0.005; % Step size
period = 500; % Length of the cycle in time
tstep_len = period/t_int; % Number of steps in a cycle
NumPeriod = 6; % Number of cycles to consider
period_steps = linspace(0,period,tstep_len+1);
timesteps = period_steps(2:end);

% Splice together the time steps to create a set of time points across the
% entire number of periods we want to run the simulation for.
for i = 1:NumPeriod
    newsteps = period_steps(2:end) + i*period;
    timesteps = [timesteps, newsteps];
end
clear newsteps

% Parameters related to variability
% Set the number of replicate simulations to run
reps = 20;
% Create the environment
rng(1) % Set seed
% Sample from standard normal (these are the time increments)
X = randn(reps,length(timesteps));
% Set the standard deviation of the process in one unit time.
sd = linspace(0,0.005,30);
% Calculate the appropriate sigma based on the time interval.
sigma = sd/sqrt(t_int);

%% Create deterministic cycles
% Create a vector of k values to consider
kvec = 0:1:2;
Nend = zeros(size(kvec));
['Solving Deterministic Model']
for i = 1:length(kvec)
    k = kvec(i);
    [T,N] = ode15s(@(t,N) PolynomialODE(t,N,stable_eq,gamma/deltan.^(2*k),k,delta,2*pi/period), timesteps, unstable_eq);
    Nend(i) = N(end);
end

['Deterministic Model Solved']

%% Loop over the different values of k and reps

% To parallelize the computation, create a matrix of all combinations of
% reps and the values of k.
[K,REP] = meshgrid(kvec, 1:reps);
% Now reshape these into vectors where the pair of vectors contains all
% values of k for all reps.
K = reshape(K, 1, length(kvec)*reps); REP = reshape(REP, 1, length(kvec)*reps);

% Create empty vectors to hold global simulation output
Nvar = NaN(length(sigma),length(K));
NFivePercent = NaN(size(Nvar));

% Create a vector of initial population sizes
Nend_vec = NaN(size(K));
for i = 1:length(kvec)
    Nend_vec(K==kvec(i)) = Nend(i);
end
['Simulating Stochastic Model']
% Run the parallel computing over all combinations of reps and k.
parfor j = 1:length(K)
    N_stoc = zeros(length(sigma),length(timesteps));
    N_stoc(:,1) = Nend_vec(j);

    k = K(j)
    % Time step loop
    for t = 2:length(timesteps)

        % Calculate the polynomial part of the model
        F = -gamma*deltan^(-2*k).*N_stoc(:,t-1).*(N_stoc(:,t-1)-stable_eq(1)).*(N_stoc(:,t-1)-stable_eq(2)).*(N_stoc(:,t-1)-unstable_eq).^(1+2*k);

        % Calculate the seasonal forcing part of the model
        Forcing = N_stoc(:,t-1).*delta*cos(2*pi*timesteps(t-1)/period);

        % Combine the deterministic and stochastic parts together to get
        % the change in population size in dt
        dN = (F + Forcing)*t_int + sigma'.*N_stoc(:,t-1)*sqrt(t_int)*X(REP(j),t-1);

        % Project the next time step
        N_stoc(:,t) = N_stoc(:,t-1) + dN;
        
    end

    
    % Now calculate the total variance
    Nvar(:,j) = var(N_stoc,1,2);
    
    [K(j),REP(j)]
end
['Stochastic Model Simulation Done']

%% Plotting results
% Generate empty vectors to hold the quantiles for the simulation results
LowerVarQuantile = NaN(length(sigma), length(kvec));
UpperVarQuantile = NaN(length(sigma), length(kvec));
MeanVarQuantile = NaN(length(sigma), length(kvec));

% Open figure plot
figure()
plot(0,0,'')
xlabel('Environmental Variability, \sigma^2');
ylabel('Population Variability, Var({\itN})');
hold on

% For each k value
for i = 1:length(kvec)
    indx = find(K == kvec(i));
    Nvar_set = Nvar(:,indx);
    % Calculate mean and 95% quantile from the replicate sims.
    LowerVarQuantile(:,i) = quantile(Nvar_set,0.025,2);
    UpperVarQuantile(:,i) = quantile(Nvar_set,0.975,2);
    MeanVarQuantile(:,i) = mean(Nvar_set,2);

    % Plot the mean and 95% outcome ranges.
    fill([sigma.^2, fliplr(sigma.^2)], [LowerVarQuantile(:,i)', fliplr(UpperVarQuantile(:,i)')], 0.7*ones(1,3));
    plot(sigma.^2, MeanVarQuantile(:,i)', 'Color', 'black', 'LineWidth', 3);
end
hold off
ax = gca; ax.FontSize = 20; ax.FontName = 'Times New Roman';

