% Plot figure as in Fig. 5 based on the data from DATA/*.mat
basePath = [fileparts(mfilename('fullpath')) filesep];
load([basePath 'DATA/VIA_Simulation.mat'])
load([basePath 'DATA/VIA_Analytical.mat'])

index=1:length(tau_N_inverse);
close all
figure( 'Position', [20 20 600 500])
a=subplot(1,2,1);
a.Position=[0.1530 0.1200 0.3046 0.7650];

Set1_Sr_Ana_cuv=fit(tau_N_inverse',MSE_S_ana','linearinterp');
Set1_Sr_Simu_cuv=fit(tau_N_inverse',MSE_S_simulation,'linearinterp');
Set1_Su_Ana_cuv=fit(tau_N_inverse',MSE_G_ana','linearinterp');
Set1_Su_Simu_cuv=fit(tau_N_inverse',MSE_G_simulation,'linearinterp');

p1=plot(tau_N_inverse(index),MSE_S_ana((index)),'-','Color',[0.5 0 0.5],'LineWidth',1.5,'MarkerSize',6);
hold on
p2=plot(tau_N_inverse(index),MSE_S_simulation((index)),'o--','Color',[0.5 0 0.5],'LineWidth',1.5,'MarkerSize',6);

grid on;
axis([tau_N_inverse(index(1)) tau_N_inverse(index(end)) -30 -10])
set(gca, 'YTick',-40:10:-10)
set(gca, 'XTick',-40:5:-20)
xlabel('$1/\tau_N$ (dB)','FontSize',22,'Interpreter','LaTex');
ylabel('MSE of $\bf S$ (dB)','FontSize',20,'Interpreter','LaTex');
set(gca,'FontSize',15);

lgd=legend([p1,p2],'K=40 - Analysis','K=40 - Simulation','FontSize',10,'Location','northeast');
set(gca,'FontSize',18);
ah=axes('position',get(gca,'position'),'visible','off');

lgd.Orientation ='horizontal';

b=subplot(1,2,2);
b.Position=[0.5934 0.1200 0.3116 0.7650];
plot(tau_N_inverse(index),MSE_G_ana((index)),'-','Color',[0.5 0 0.5],'LineWidth',1.5,'MarkerSize',6);
hold on
plot(tau_N_inverse(index),MSE_G_simulation((index)),'o--','Color',[0.5 0 0.5],'LineWidth',1.5,'MarkerSize',6);

grid on;
axis([tau_N_inverse(index(1)) tau_N_inverse(index(end)) -30 -10])
set(gca, 'YTick',-40:10:-10)
set(gca, 'XTick',-40:5:-20)
xlabel('$1/\tau_N$ (dB)','FontSize',22,'Interpreter','LaTex');
ylabel('MSE of $\bf G$ (dB)','FontSize',20,'Interpreter','LaTex');
set(gca,'FontSize',15);
