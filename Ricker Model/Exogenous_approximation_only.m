clear
clc
% Deterministic and stochastic Ricker models in the three cases.

%% Global parameters
abar_vec = linspace(0.5,1.5,3);
colors = viridis(length(abar_vec)+2);

Nbar = 1;
gen = 100000;
Ninit = 0.95*Nbar;

%% Exogenous Cycle Effect of abar
sigma_a = sqrt(linspace(0,0.5,20));
cycle_length = 20;

figure
subplot(1,2,1)
hold on
for i = 1:length(abar_vec)
    abar = abar_vec(i);
    at = abar + sqrt(2*sigma_a'.^2)*sin(2*pi*[1:gen]/cycle_length);
    b = abar;
    N = zeros(length(sigma_a),gen);
    N(:,1) = Ninit;
    
    for t = 2:gen
        N(:,t) = N(:,t-1).*exp(at(:,t-1) - b*N(:,t-1));
    end
    
    coef = History_Coeff(abar, cycle_length);
    
    varNpredict =  Nbar^2/(1 - (1-abar)^2)*(sigma_a'.^2*(1+coef));
    varNactual = var(N',1);
    p = plot(sigma_a.^2, varNpredict, ':');
    p.LineWidth = 3; p.Color = colors(i+1,:);
    s = scatter(sigma_a.^2, varNactual, 'filled');
    s.SizeData = 100; s.MarkerEdgeColor = 0.5*colors(i+1,:);
    s.MarkerFaceColor = colors(i+1,:);
    s.HandleVisibility = 'off';
    
    lgd = legend(strsplit(num2str(abar_vec)));
    title(lgd, 'Average Growth Rate, $\bar{a}$', 'Interpreter', 'Latex');
    
    
end
hold off

xlabel('\sigma_a^2, Exogenous Driver Variance');
ylabel('Var({\itN}), Population Size Variance');
title({'Approximation effect', 'of average gorwth rate'});
ax = gca; ax.FontSize = 25; ax.FontName = 'Times New Roman';
xlim([0,max(sigma_a.^2)])

%% Do it again with an effect of cycle length
clear
clc
% Deterministic and stochastic Ricker models in the three cases.

%% Global parameters
abar = 1.5;
cycle_length_vec = [10,20,50];

colors = viridis(length(cycle_length_vec)+2);

Nbar = 1;
gen = 100000;
plotgen = 100;
tinit = 2000;
Ninit = 0.95*Nbar;

%% Exogenous Cycle Effect of abar
sigma_a = sqrt(linspace(0,0.5,20));

subplot(1,2,2)
hold on
for i = 1:length(cycle_length_vec)
    cycle_length = cycle_length_vec(i);
    at = abar + sqrt(2*sigma_a'.^2)*sin(2*pi*[1:gen]/cycle_length);
    b = abar;
    N = zeros(length(sigma_a),gen);
    N(:,1) = Ninit;
    
    for t = 2:gen
        N(:,t) = N(:,t-1).*exp(at(:,t-1) - b*N(:,t-1));
    end
    
    coef = History_Coeff(abar, cycle_length);
    
    varNpredict =  Nbar^2/(1 - (1-abar)^2)*(sigma_a'.^2*(1+coef));
    varNactual = var(N',1);
    p = plot(sigma_a.^2, varNpredict, ':');
    p.LineWidth = 3; p.Color = colors(i+1,:);
    s = scatter(sigma_a.^2, varNactual, 'filled');
    s.SizeData = 100; s.MarkerEdgeColor = 0.5*colors(i+1,:);
    s.MarkerFaceColor = colors(i+1,:);
    s.HandleVisibility = 'off';

    lgd = legend(strsplit(num2str(cycle_length_vec)));
    title(lgd, 'Cycle length, \Omega');
    
end
hold off

xlabel('\sigma_a^2, Exogenous Driver Variance');
ylabel('Var({\itN}), Population Size Variance');
title({'Approximation effect', 'of exogenous cycle length'});
ax = gca; ax.FontSize = 25; ax.FontName = 'Times New Roman';
xlim([0,max(sigma_a.^2)])
