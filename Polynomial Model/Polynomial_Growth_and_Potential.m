% This 

n1 = 8; n2 = 12; nhat = mean([n1,n2]); c = 0.01; v = 0;
deltan = n1 - nhat;
N = linspace(n1-0.5,n2+0.5,1000);

colors = viridis(5);

figure(1)
subplot(1,2,1)
scatter([n1,n2,nhat],zeros(1,3),'HandleVisibility', 'off');
yline(0,'--', 'HandleVisibility', 'off');
xlabel('$n$', 'Interpreter', 'Latex'); 
ylabel('$\frac{1}{n}\frac{dn}{dt}$', 'Interpreter', 'Latex');
ax = gca; ax.FontSize = 30;
ax.FontName = 'Times New Roman';
ax.XTick = [0,n1,nhat,n2];
ax.XTickLabels = {'0','$\hat{n}_1$','$\tilde{n}$','$\hat{n}_2$'};
title('Per-Capita Growth Rate')
ax.TickLabelInterpreter = 'Latex';
txt = text(0,1,'a', 'Units', 'normalized'); 
txt.FontSize = 50; txt.HorizontalAlignment = 'right'; txt.VerticalAlignment = 'bottom';
ylim(0.06*[-1,1])
xlim([n1-0.25,n2+0.25])

subplot(1,2,2)
plot([n1,n2,nhat],zeros(1,3),'Color','none', 'HandleVisibility', 'off');
xlabel('$n$', 'Interpreter', 'Latex'); 
ylabel('$U(N)$', 'Interpreter', 'Latex');
ax = gca; ax.FontSize = 30;
ax.FontName = 'Times New Roman';
ax.XTick = [0,n1,nhat,n2];
ax.XTickLabels = {'0','$\hat{n}_1$','$\tilde{n}$','$\hat{n}_2$'};
title('Dynamical Potential')
ax.TickLabelInterpreter = 'Latex';
txt = text(0,1,'b', 'Units', 'normalized'); 
txt.FontSize = 50; txt.HorizontalAlignment = 'right';
txt.VerticalAlignment = 'bottom';

ylim(0.5*[-0.025,1.1])
xlim([n1-0.5,n2+0.5])

for k = 0:1:2
    [F,U] = AltStatePotential(N,n1,n2,nhat,k,c*deltan^(-2*k),v);
    U = U-min(U);

    subplot(1,2,1)
    hold on; plot(N,F./N, 'LineWidth', 3, 'Color', colors(k+2,:)); hold off;
    
    subplot(1,2,2)
    hold on; plot(N,U, 'LineWidth', 3, 'Color', colors(k+2,:)); hold off;
end


subplot(1,2,1)
legend('{\itk} = 0', '{\itk} = 1', '{\itk} = 2')
% 
% 
% k = 0;
% 
% 
% N = linspace(n1-0.5,n2+0.5,1000);
% 
% figure()
% subplot(2,3,1:3)
% scatter([n1,n2,nhat],zeros(1,3),'HandleVisibility', 'off');
% hold on
% yline(0,'--');
% hold off
% xlabel('\itN'); ylabel({'{\itF}({\itN})','Population Growth Rate'});
% ax = gca; ax.FontSize = 20;
% title('Population Growth Rate')
% 
% 
% for i = 1:length(c)
%     [F,U] = AltStatePotential(N,n1,n2,nhat,k,c(i));
%     subplot(2,3,1:3)
%     hold on; plot(N,F); hold off;
%     
%     subplot(2,3,3+i)
%     plot(N,U);
%     xlabel('\itN'); ylabel('{\itU}({\itN})');
%     ax = gca; ax.FontSize = 20;
%     title(['Potential, {\itc} = ',num2str(c(i))])
% end
