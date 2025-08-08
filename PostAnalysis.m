% file1 = 'Varying_15_0.5.mat';
% file2 = 'Constant_15_0.5.mat';
% file1 = 'Varying_15_1.5.mat';
% file2 = 'Constant_15_1.5.mat';
% file1 = 'Varying_3.75_1.mat';
% file2 = 'Constant_3.75_1.mat';
file1 = 'Varying_3.75_2.mat';
file2 = 'Constant_3.75_2.mat';
% file1 = 'Varying_3.75_2_turbulence.mat';
% file2 = 'Constant_3.75_2_turbulence.mat';
% file2 = 'Varying_3.75_2.mat';
% file1 = 'Varying_3.75_2_turbulence.mat';

dirct = 'C:\Users\Jordy\Documents\RIT\Nebulizer research\My Documents\Papers\Bulletin of Mathematical Biology\Images\png files';

save_tog = 0;

load(file1)
dep1 = f;
sus1 = g;

load(file2)
dep2 = f;
sus2 = g;

x1 = linspace(.01,.1,10);
x2 = linspace(.2,1,9);
x3 = linspace(2,10,9);
x = [x1 x2 x3]; 

% Define particle size ranges
x_ultra = x(x<=0.1);
x_fine = x(x>0.1 & x<1);
x_med = x(x>=1 & x<5);
x_coarse = x(x>=5);

Out = Analyze(x,dep1,dep2);
% Out = Analyze(x,sus1,sus2);
% Out = Analyze(x,1-(dep1(5,:)+sus1(10,:)),1-(dep2(5,:)+sus2(10,:)));

set(groot,'defaultLineLineWidth',1, 'defaulttextInterpreter', 'latex', ...
     'defaultLegendInterpreter', 'latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultAxesFontSize', 16)
set(groot,'defaultFigureColor',[1 1 1])
set(groot,'defaultLegendFontSize',12,'DefaultLegendFontSizeMode','manual')

% DEPOSITION PLOTS=========================================================
figure(1)
% total deposition fraction after repeatability
hold on
p = plot(x,dep1(5,:),'-k>','MarkerSize',6,'DisplayName','varying compliance');
p.MarkerFaceColor='k';
p = plot(x,dep2(5,:),'-k>','MarkerSize',6,'DisplayName','constant compliance');
p.MarkerFaceColor='w';
set(gca,'XScale','log')
set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
set(gcf,'Color',[1 1 1])
ylabel('Deposition Fraction, $\eta_{total}$')
xlabel('Particle Size ($\mu$m)')
legend box off location northwest
grid on
ylim([0 1])
if save_tog==1
    saveas(gcf, fullfile(dirct, ['constantVSvarying_dep_tot_RR',num2str(res.RR),'_TV',num2str(vol.TV),'.png']));
end

figure(2)
% TB deposition after repeatability
hold on
p = plot(x,dep1(1,:)+dep1(2,:),'-ks','MarkerSize',6,'DisplayName','varying compliance');
p.MarkerFaceColor='k';
p = plot(x,dep2(1,:)+dep2(2,:),'-ks','MarkerSize',6,'DisplayName','constant compliance');
p.MarkerFaceColor='w';
set(gca,'XScale','log')
set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
set(gcf,'Color',[1 1 1])
ylabel('Tracheobronchial Deposition')
legend box off location northwest
grid on
ylim([0 1])
if save_tog==1
    saveas(gcf, fullfile(dirct, ['constantVSvarying_dep_TB_RR',num2str(res.RR),'_TV',num2str(vol.TV),'.png']));
end

figure(3)
% respiratory bronchioles deposition
hold on
p = plot(x,dep1(3,:),'-ko','MarkerSize',6,'DisplayName','varying compliance');
p.MarkerFaceColor='k';
p4 = plot(x,dep2(3,:),'-ko','MarkerSize',6,'DisplayName','constant compliance');
p4.MarkerFaceColor='w';
set(gca,'XScale','log')
set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
set(gcf,'Color',[1 1 1])
ylabel('Respiratory Bronchioles Deposition')
xlabel('Particle Size ($\mu$m)')
legend box off location northwest
grid on
ylim([0 1])
if save_tog==1
    saveas(gcf, fullfile(dirct, ['constantVSvarying_dep_rb_RR',num2str(res.RR),'_TV',num2str(vol.TV),'.png']));
end

figure(4)
% alveolar deposition
hold on
p = plot(x,dep1(4,:),'-kd','MarkerSize',6,'DisplayName','varying compliance');
p.MarkerFaceColor='k';
p = plot(x,dep2(4,:),'-kd','MarkerSize',6,'DisplayName','constant compliance');
p.MarkerFaceColor='w';
set(gca,'XScale','log')
set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
set(gcf,'Color',[1 1 1])
ylabel('Alveolar Deposition')
xlabel('Particle Size ($\mu$m)')
legend box off location northwest
grid on
ylim([0 1])
if save_tog==1
    saveas(gcf, fullfile(dirct, ['constantVSvarying_dep_alv_RR',num2str(res.RR),'_TV',num2str(vol.TV),'.png']));
end

figure(5)
% respiratory deposition
hold on
p = plot(x,dep1(3,:)+dep1(4,:),'-kd','MarkerSize',6,'DisplayName','varying compliance');
p.MarkerFaceColor='k';
p = plot(x,dep2(3,:)+dep2(4,:),'-kd','MarkerSize',6,'DisplayName','constant compliance');
p.MarkerFaceColor='w';
set(gca,'XScale','log')
set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
set(gcf,'Color',[1 1 1])
ylabel('Pulmonary Deposition')
xlabel('Particle Size ($\mu$m)')
legend box off location northwest
grid on
ylim([0 1])
if save_tog==1
    saveas(gcf, fullfile(dirct, ['constantVSvarying_dep_P_RR',num2str(res.RR),'_TV',num2str(vol.TV),'.png']));
end

% RETENTION/EXHALE PLOTS==========================================================
figure(6)
% total exhaled
hold on
p4 = plot(x,1-(dep1(5,:)+sus1(10,:)),'-k<','MarkerSize',6,'DisplayName','varying compliance');
p4.MarkerFaceColor='k';
p4 = plot(x,1-(dep2(5,:)+sus2(10,:)),'-k<','MarkerSize',6,'DisplayName','constant compliance');
p4.MarkerFaceColor='w';
set(gca,'XScale','log')
set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
set(gcf,'Color',[1 1 1])
ylabel('Total Exhaled Fraction')
xlabel('Particle Size ($\mu$m)')
legend box off location northwest
grid on
ylim([0 1])
if save_tog==1
    saveas(gcf, fullfile(dirct, ['constantVSvarying_exh_RR',num2str(res.RR),'_TV',num2str(vol.TV),'.png']));
end

% figure(7)
% % TB retention after repeatability
% hold on
% if comp_tog==0
%     p3 = plot(dp_list*10^4,f(1,:)+f(2,:)+g(6,:)+g(7,:),'-ks','MarkerSize',6,'DisplayName','constant compliance');
%     p3.MarkerFaceColor='w';
% elseif comp_tog==1
%     p3 = plot(dp_list*10^4,f(1,:)+f(2,:)+g(6,:)+g(7,:),'-ks','MarkerSize',6,'DisplayName','varying compliance');
%     p3.MarkerFaceColor='k';
% end
% set(gca,'XScale','log')
% set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
% set(gcf,'Color',[1 1 1])
% ylabel('Tracheobronchial Retention')
% legend box off location northwest
% grid on
% ylim([0 1])
% 
% figure(8)
% % respiratory bronchioles retention
% hold on
% if comp_tog==0
%     p4 = plot(dp_list*10^4,f(3,:)+g(8,:),'-ko','MarkerSize',6,'DisplayName','constant compliance');
%     p4.MarkerFaceColor='w';
% elseif comp_tog==1
%     p4 = plot(dp_list*10^4,f(3,:)+g(8,:),'-ko','MarkerSize',6,'DisplayName','varying compliance');
%     p4.MarkerFaceColor='k';
% end
% set(gca,'XScale','log')
% set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
% set(gcf,'Color',[1 1 1])
% ylabel('Respiratory Bronchioles Retention')
% xlabel('Particle Size ($\mu$m)')
% legend box off location northwest
% grid on
% ylim([0 1])
% 
% figure(9)
% % alveolar retention
% hold on
% if comp_tog==0
%     p4 = plot(dp_list*10^4,f(4,:)+g(9,:),'-kd','MarkerSize',6,'DisplayName','constant compliance');
%     p4.MarkerFaceColor='w';
% elseif comp_tog==1
%     p4 = plot(dp_list*10^4,f(4,:)+g(9,:),'-kd','MarkerSize',6,'DisplayName','varying compliance');
%     p4.MarkerFaceColor='k';
% end
% set(gca,'XScale','log')
% set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
% set(gcf,'Color',[1 1 1])
% ylabel('Alveolar Retention')
% xlabel('Particle Size ($\mu$m)')
% legend box off location northwest
% grid on
% ylim([0 1])
% 
% figure(10)
% % respiratory retention
% hold on
% if comp_tog==0
%     p4 = plot(dp_list*10^4,f(3,:)+f(4,:)+g(8,:)+g(9,:),'-kd','MarkerSize',6,'DisplayName','constant compliance');
%     p4.MarkerFaceColor='w';
% elseif comp_tog==1
%     p4 = plot(dp_list*10^4,f(3,:)+f(4,:)+g(8,:)+g(9,:),'-kd','MarkerSize',6,'DisplayName','varying compliance');
%     p4.MarkerFaceColor='k';
% end
% set(gca,'XScale','log')
% set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
% set(gcf,'Color',[1 1 1])
% ylabel('Pulmonary Retention')
% xlabel('Particle Size ($\mu$m)')
% legend box off location northwest
% grid on
% ylim([0 1])

% SUSPENSION PLOTS==========================================================
figure(11)
% total suspension
hold on
p = plot(x,sus1(10,:),'-k<','MarkerSize',6,'DisplayName','varying compliance');
p.MarkerFaceColor='k';
p = plot(x,sus2(10,:),'-k<','MarkerSize',6,'DisplayName','constant compliance');
p.MarkerFaceColor='w';
set(gca,'XScale','log')
set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
set(gcf,'Color',[1 1 1])
ylabel('Total Suspended Fraction')
xlabel('Particle Size ($\mu$m)')
legend box off location northwest
grid on
ylim([0 1])
if save_tog==1
    saveas(gcf, fullfile(dirct, ['constantVSvarying_sus_tot_RR',num2str(res.RR),'_TV',num2str(vol.TV),'.png']));
end

figure(12)
hold on
% TB suspension after repeatability
p = plot(x,sus1(6,:)+sus1(7,:),'-ks','MarkerSize',6,'DisplayName','varying compliance');
p.MarkerFaceColor='k';
p = plot(x,sus2(6,:)+sus2(7,:),'-ks','MarkerSize',6,'DisplayName','constant compliance');
p.MarkerFaceColor='w';
set(gca,'XScale','log')
set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
set(gcf,'Color',[1 1 1])
ylabel('Tracheobronchial Suspension')
legend box off location northwest
grid on
ylim([0 1])
if save_tog==1
    saveas(gcf, fullfile(dirct, ['constantVSvarying_sus_TB_RR',num2str(res.RR),'_TV',num2str(vol.TV),'.png']));
end

figure(13)
hold on
% respiratory bronchioles suspension
p = plot(x,sus1(8,:),'-ko','MarkerSize',6,'DisplayName','varying compliance');
p.MarkerFaceColor='k';
p = plot(x,sus2(8,:),'-ko','MarkerSize',6,'DisplayName','constant compliance');
p.MarkerFaceColor='w';
set(gca,'XScale','log')
set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
set(gcf,'Color',[1 1 1])
ylabel('Respiratory Bronchioles Suspension')
xlabel('Particle Size ($\mu$m)')
legend box off location northwest
grid on
ylim([0 1])
if save_tog==1
    saveas(gcf, fullfile(dirct, ['constantVSvarying_sus_rb_RR',num2str(res.RR),'_TV',num2str(vol.TV),'.png']));
end

figure(14)
hold on
% alveolar suspension
p = plot(x,sus1(9,:),'-kd','MarkerSize',6,'DisplayName','varying compliance');
p.MarkerFaceColor='k';
p = plot(x,sus2(9,:),'-kd','MarkerSize',6,'DisplayName','constant compliance');
p.MarkerFaceColor='w';
set(gca,'XScale','log')
set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
set(gcf,'Color',[1 1 1])
ylabel('Alveolar Suspension')
xlabel('Particle Size ($\mu$m)')
legend box off location northwest
grid on
ylim([0 1])
if save_tog==1
    saveas(gcf, fullfile(dirct, ['constantVSvarying_sus_alv_RR',num2str(res.RR),'_TV',num2str(vol.TV),'.png']));
end

figure(15)
hold on
% respiratory suspension
p = plot(x,sus1(8,:)+sus1(9,:),'-kd','MarkerSize',6,'DisplayName','varying compliance');
p.MarkerFaceColor='k';
p = plot(x,sus2(8,:)+sus2(9,:),'-kd','MarkerSize',6,'DisplayName','constant compliance');
p.MarkerFaceColor='w';
set(gca,'XScale','log')
set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
set(gcf,'Color',[1 1 1])
ylabel('Pulmonary Suspension')
xlabel('Particle Size ($\mu$m)')
legend box off location northwest
grid on
ylim([0 1])
if save_tog==1
    saveas(gcf, fullfile(dirct, ['constantVSvarying_sus_P_RR',num2str(res.RR),'_TV',num2str(vol.TV),'.png']));
end

figure(16)
hold on
% ybars = [ avg_diff_TB_ultra avg_diff_P_ultra;...
%     avg_diff_tot_fine avg_diff_TB_fine avg_diff_P_fine;...
%     avg_diff_tot_med avg_diff_TB_med avg_diff_P_med;...
%     avg_diff_tot_coarse avg_diff_TB_coarse avg_diff_P_coarse];

% b1 = barh(ybars,'BaseValue',0,'FaceAlpha',1,'LineStyle','-');
b1 = barh(Out,'BaseValue',0,'FaceAlpha',1,'LineStyle','-');
% set(get(b1,'Children'),'FaceAlpha',0.3)
set(gca, 'YTick',[1 2 3 4])
% set(gca, 'XTick',[-0.25 -0.15 -0.05 0 0.05 0.15 0.25],'XGrid','on')
% set(gca, 'XTick',-0.06:.02:0.06,'XGrid','on')

yticklabels({'ultrafine' 'fine' 'medium' 'coarse'})
ytickangle(90)
xtickangle(45)
xline(0.02,'--k','LineWidth',1.5)
xline(-0.02,'--k','LineWidth',1.5)
legend({'TB','rb','alv','total'})
% legend({'total','alv','rb','TB'})
legend box off location southwest
xlim([-0.06 0.06])
xlabel('$mean(\tilde{\eta}_i-\eta_i$)')
ylabel('Particle Size')

if save_tog==1
    saveas(gcf, fullfile(dirct, ['constantVSvarying_bar_RR',num2str(res.RR),'_TV',num2str(vol.TV),'.png']));
end
%============================================================
function Out = Analyze(x,array1,array2)

% The input array could be either deposition or suspension fractions
u1 = array1(1,:);
c1 = array1(2,:);
rb1 = array1(3,:);
alv1 = array1(4,:);
tot1 = array1(5,:);

u2 = array2(1,:);
c2 = array2(2,:);
rb2 = array2(3,:);
alv2 = array2(4,:);
tot2 = array2(5,:);

% for ultrafine particles
avg_diff_TB_ultra = MeanDifference(u1+c1,u2+c2,x<=0.1);
avg_diff_rb_ultra = MeanDifference(rb1,rb2,x<=0.1);
avg_diff_alv_ultra = MeanDifference(alv1,alv2,x<=0.1);
avg_diff_tot_ultra = MeanDifference(tot1,tot2,x<=0.1);

% for fine particles
avg_diff_TB_fine = MeanDifference(u1+c1,u2+c2,x>0.1 & x<1);
avg_diff_rb_fine = MeanDifference(rb1,rb2,x>0.1 & x<1);
avg_diff_alv_fine = MeanDifference(alv1,alv2,x>0.1 & x<1);
avg_diff_tot_fine = MeanDifference(tot1,tot2,x>0.1 & x<1);

% for medium particles
avg_diff_TB_med = MeanDifference(u1+c1,u2+c2,x>=1 & x<5);
avg_diff_rb_med = MeanDifference(rb1,rb2,x>=1 & x<5);
avg_diff_alv_med = MeanDifference(alv1,alv2,x>=1 & x<5);
avg_diff_tot_med = MeanDifference(tot1,tot2,x>=1 & x<5);

% for coarse particles
avg_diff_TB_coarse = MeanDifference(u1+c1,u2+c2,x>=5);
avg_diff_rb_coarse = MeanDifference(rb1,rb2,x>=5);
avg_diff_alv_coarse = MeanDifference(alv1,alv2,x>=5);
avg_diff_tot_coarse = MeanDifference(tot1,tot2,x>=5);

Out = ones(4,4);
% Out(1,:) = [avg_diff_TB_ultra avg_diff_rb_ultra avg_diff_alv_ultra avg_diff_tot_ultra];
% Out(2,:) = [avg_diff_TB_fine avg_diff_rb_fine avg_diff_alv_fine avg_diff_tot_fine];
% Out(3,:) = [avg_diff_TB_med avg_diff_rb_med avg_diff_alv_med avg_diff_tot_med];
% Out(4,:) = [avg_diff_TB_coarse avg_diff_rb_coarse avg_diff_alv_coarse avg_diff_tot_coarse];

Out(1,:) = [avg_diff_TB_ultra avg_diff_rb_ultra avg_diff_alv_ultra avg_diff_tot_ultra];
Out(2,:) = [avg_diff_TB_fine avg_diff_rb_fine avg_diff_alv_fine avg_diff_tot_fine];
Out(3,:) = [avg_diff_TB_med avg_diff_rb_med avg_diff_alv_med avg_diff_tot_med];
Out(4,:) = [avg_diff_TB_coarse avg_diff_rb_coarse avg_diff_alv_coarse avg_diff_tot_coarse];

function A = MeanDifference(v1,v2,size)
    A = mean(v1(size)-v2(size));
%     A = mean((v1(size)-v2(size))./v2(size));
end

end