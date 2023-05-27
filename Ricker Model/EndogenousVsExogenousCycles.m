clear
clc
%% Simple Density-dependence and Exogenous Cycle
period = 20;
gen = 100*period;
abar = 1; sigma_a = sqrt([0,0.25]);
reps = 1000;
K = 1;
b = abar./K;
sigma = sqrt(linspace(0,0.1,100));
X = randn([reps,gen]);

figure(1)
figure(2)

TitleString = {'Simple Density-Dependence', 'Exogenous Cycle'};
TextString = {'a','b'};
for i = 1:2
    a = abar + sqrt(2)*sigma_a(i)*sin([1:gen]*2*pi/period);
    varN = zeros(reps,length(sigma));
    
    Ndet = zeros(1,gen);
    Ndet(1) = K(1);
    
    for t = 2:gen
         Ndet(t) = Ndet(t-1).*exp(a(t-1) - b.*Ndet(t-1));
    end
    
    for r = 1:reps
        Nstoc = zeros(length(sigma),gen);
        Nstoc(:,1) = K(1);
        
        for t = 2:gen
            Nstoc(:,t) = Nstoc(:,t-1).*exp(a(t-1) - b.*Nstoc(:,t-1) + sigma'*X(r,t-1));
        end
        varN(r,:) = var(Nstoc',1);
        StocVar(:,r) = mean((Nstoc - ones(length(sigma),1)*Ndet).^2,2);
    end
    
    DetVar = mean((Ndet - K).^2);
    IntVar = varN - StocVar' - DetVar;
    
    lowquant = quantile(varN,0.025);
    uppquant = quantile(varN,0.975);
    figure(1)
    subplot(1,3,i)
    fill([sigma.^2,fliplr(sigma.^2)], [uppquant, fliplr(lowquant)],0.5*ones(1,3));
    hold on
    plot(sigma.^2, mean(varN), 'LineWidth', 2, 'Color', zeros(1,3));
    hold off
    xlabel('Environmental Variance, \sigma^2');
    ylabel({'Population Variance','Var({\itn})'});
    title(TitleString(i))
    ax = gca; ax.FontSize = 20; ax.FontName = 'Times New Roman';
    t = text(0,1,TextString(i), 'Units', 'Normalized');
    t.FontSize = 40; t.FontName = 'Helvetica';
    t.HorizontalAlignment = 'right'; t.VerticalAlignment = 'bottom';
    
    figure(2)
    hold on
    plot(sigma.^2, mean(IntVar))
    hold off
    xlabel('Environmental Variance, \sigma^2');
    ylabel('Environment x Det Interaction');
        
end

%%
a = log(3)*2; b = a./K;
varN = zeros(reps,length(sigma));

Ndet = zeros(1,gen);
Ndet(1) = 1.5*K(1);

for t = 2:gen
    Ndet(t) = Ndet(t-1).*exp(a - b.*Ndet(t-1));
end

DetVar = var(Ndet,1);
StocVar = zeros(reps,length(sigma));
for r = 1:reps
    Nstoc = zeros(length(sigma),gen);
    Nstoc(:,1) = K(1)*1.1;
    
    for t = 2:gen
        Nstoc(:,t) = Nstoc(:,t-1).*exp(a - b.*Nstoc(:,t-1) + sigma'*X(r,t-1));
    end
    varN(r,:) = var(Nstoc',1);
    StocVar(r,:) = mean((Nstoc' - Ndet'*ones(1,length(sigma))).^2);
end

IntVar = varN - StocVar - DetVar;

figure(2)
hold on
plot(sigma.^2, mean(IntVar));
hold off

figure(1)
lowquant = quantile(varN,0.025);
uppquant = quantile(varN,0.975);
subplot(1,3,3)
fill([sigma.^2,fliplr(sigma.^2)], [uppquant, fliplr(lowquant)],0.5*ones(1,3));
hold on
plot(sigma.^2, mean(varN), 'LineWidth', 2, 'Color', zeros(1,3));
hold off
xlabel('Environmental Variance, \sigma^2');
ylabel({'Population Variance','Var({\itn})'});
title('Endogenous Cycles')
ax = gca; ax.FontSize = 20; ax.FontName = 'Times New Roman';
t = text(0,1,'c', 'Units', 'Normalized');
t.FontSize = 40; t.FontName = 'Helvetica';
t.HorizontalAlignment = 'right'; t.VerticalAlignment = 'bottom';
