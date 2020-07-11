% Plot purple curves as in Fig. 8 based on the data from DATA/*.mat

warning('off')
basePath = [fileparts(mfilename('fullpath')) filesep];
load([basePath 'DATA/VIB_Simulation.mat'])

snr=70:5:110;
index=1:length(snr);


fitmethod='smoothingspline';
MP_RB=fit(snr',NMSE_H_RB',fitmethod);
MP_UR=fit(snr',NMSE_H_UR',fitmethod);

close all
figure( 'Position', [20 20 650 500])
a=subplot(1,2,1);
a.Position=[0.1230 0.1200 0.3346 0.7250];

% hold on
plot(snr(index),MP_RB(snr(index)),'o-','Color',[0.5 0 0.5],'LineWidth',1.5,'MarkerSize',7);
hold on
grid on;


set(gca, 'YTick',-30:10:0)
axis([70 110 -30 5])
set(gca, 'XTick',70:10:110)
xlabel({'$1/\tau_N$ (dB)'},'Interpreter','latex')




ylabel('NMSE of $\bf{H}_{RB}$ (dB)','Interpreter','latex','FontSize',18);
set(gca,'FontSize',15);

b=subplot(1,2,2);
b.Position=[0.5834 0.1200 0.3346 0.7250];
p1=plot(snr(index),MP_UR(snr(index)),'o-','Color',[0.5 0 0.5],'LineWidth',1.5,'MarkerSize',7);
hold on
grid on;
% set(gca, 'XTick',snr(index))
set(gca, 'YTick',-30:10:10)
axis([70 110 -30 5])
set(gca, 'YTick',-30:10:0)
set(gca, 'XTick',70:10:110)
xlabel({'$1/\tau_N$ (dB)'},'Interpreter','latex')
ylabel('Ave. NMSE of $\bf{h}_{UR,k}$ (dB)','Interpreter','latex','FontSize',18);
set(gca,'FontSize',15);
lgd=legend(p1,{'Algorithm 1'},'FontSize',18,'Location','northeast');
ah=axes('position',get(gca,'position'),'visible','off');
set(gca,'FontSize',15);
lgd.Orientation ='horizontal';
lgd.Position=[0.1120 0.8990 0.8 0.1];