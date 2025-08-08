function [f,g,t,s] = Solver(vol,res,Y,mech,comp_tog,turb_tog)
%==========================================================================
% This function returns the solutions and creates the deposition/suspension
% plots.
%==========================================================================

% FOR FIGURES--------------------------------------------------------------
set(groot,'defaultLineLineWidth',1, 'defaulttextInterpreter', 'latex', ...
     'defaultLegendInterpreter', 'latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultAxesFontSize', 18)
set(groot,'defaultFigureColor',[1 1 1])
set(groot,'defaultLegendFontSize',14,'DefaultLegendFontSizeMode','manual')
%--------------------------------------------------------------------------

% INITIAL CONDITIONS
% IC = [0 0.96*vol.FRC 0 vol.FRC]; % [Q Vra Pve V]
IC = [0 0.96*vol.FRC vol.FRC]; % [Q Vra V] Initial conditions

% STORE REMAINING PARAMETERS IN STRUCTURE
par = Parameters(vol,Y,turb_tog,IC);

% TIME SETTINGS
dt = 10^-3; % time step
tvals = 0:dt:res.T; % times at which to solve
% tvals = 0:dt:res.TI;

% GET CONSTANT FLOWS AND VOLUMES-------------------------------------------
% Before running the full simulation, we need constant flows and volumes to
% compute the probabilities. Run the breathing mechanics model (without
% deposition) for one breathing cycle and then compute the average flows
% and volumes.
options = odeset('RelTol',10^-10);
[tt,ss] = ode15s(@(tt,ss) LungMechanicsODEs(tt,ss,par,Y,res,vol,comp_tog,IC),tvals,IC,options); % solve breathing mechanics model without deposition
[Qin,Qex,Vin,Vex] = ConstantFlow(tt,ss,Y,vol,par,res,comp_tog,IC);
% 1/(tt(end)-tt(1))*trapz(tt,ss(:,1))
% stop
% figure(90)
% hold on
% plot(tt,ss(:,1))
% max(ss(:,1))
% % stop
% figure(91)
% hold on
% plot(tt,ss(:,3))
% 
% max(ss(:,3))-min(ss(:,3))
% plot(tt,del2(ss(:,2),tt))
% stop

% for i = 1:length(tt)-1
%     diff(i) = (ss(i,2)-ss(i+1,2))/tt(i);
% end
% figure(90)
% plot(tt,ss(:,2))
% stop



% % Rebuild PV curves
% Q = ss(:,1);
% VA = ss(:,2);
% V = ss(:,3);
% % [~,~,~,Vc,Vcmax,~,~,~,~,~,~,~,~,~,~,~,~,Vu] = Scaling(Y,vol,V,VA);
% Vrange = vol.RV:.1:vol.TLC;
% 
% % MORPHOMETRY
% % save each column as array
% g = Y{:,1}; % generation number
% n = Y{:,2}; % number of airways in each generation
% Ldata = Y{:,3}; % length of airways in each generation in cm
% rdata = Y{:,4}/2; % radius of airways in each generation in cm
% Ba = Y{:,5}.*pi/180; % branching angle converted FROM degrees TO radians
% Ga = Y{:,6}.*pi/180; % gravity angle converted FROM degrees TO radians
% Vdata = Y{:,8}/1e3; % volume in L
% 
% c1_ind = 1:2;
% c2_ind = 3:16;  
% c3_ind = 17:length(g)-1; 
% c4_ind = length(g);
% 
% % Scale to our desired lung volume to initialize
% sf = (vol.TLC/sum(Vdata))^(1/3); % scaled to TLC
% % sf = (vol.FRC/sum(Vdata))^(1/3); % scaled to FRC
% r0 = rdata*sf; % cm
% L0 = Ldata*sf; % cm
% 
% Vu = sum(pi*r0(c1_ind).^2.*L0(c1_ind).*n(c1_ind))/1e3; % L
% Vc = V-(VA+Vu);
% Vcmax = sum(pi*r0(c2_ind).^2.*L0(c2_ind).*n(c2_ind))/1e3; % L
% 
% Compliance(ss,par,vol,Vc,Vcmax)
% % stop
% % actual simulation
% % Vc0 = 0.04*vol.FRC-Vu;
% V0 = vol.FRC-par.Ccw*par.Ppl;
% VA0 = 0.96*vol.FRC;
% VuA = VA0+par.Cra*par.Ppl;
% if comp_tog==0
% %     Pcw = par.Ppl+(V-vol.FRC)/par.Ccw;
% %     Pcw_range = par.Ppl+(Vrange-vol.FRC)/par.Ccw;
%     Pcw = (1/par.Ccw)*(V-V0);
%     Pcw_range = (1/par.Ccw)*(Vrange-V0);
% %     Ptm = -par.Ppl+(Vc-Vc0)/par.Cc;
% %     Pel = -par.Ppl+(V-0.96*vol.FRC)/par.Cra;
% %     Pel_range = -par.Ppl+(Vrange-0.96*vol.FRC)/par.Cra;
%     Pel = (V-VuA)/par.Cra;
%     Pel_range = (Vrange-VuA)/par.Cra;
% elseif comp_tog==1
%     Pcw = par.Acw-par.Bcw*log((vol.TLC-vol.RV)./(V-vol.RV)-0.999); % cmh2o
%     Pcw_range = par.Acw-par.Bcw*log((vol.TLC-vol.RV)./(Vrange-vol.RV)-0.999); % cmh2o
% %     Ptm = par.Ac-par.Bc2*log(Vcmax./Vc-0.999);
%     Pel = par.Al*exp(par.Kl*V)+par.Bl;
%     Pel_range = par.Al*exp(par.Kl*Vrange)+par.Bl;
% end
% 
% 
% figure(80)
% hold on
% yline(max(V)/max(Vrange),'--k','linewidth',1.5,'HandleVisibility','off')
% yline(min(V)/max(Vrange),'--k','linewidth',1.5,'HandleVisibility','off')
% % pPcw = plot(Pcw,V,'-r','DisplayName','nonlinear compliance','linewidth',2);
% % plot(Pcw_range,Vrange/max(Vrange),'-k','HandleVisibility','off');
% % plot(Pcw,V/max(Vrange),'-r','HandleVisibility','off');
% % plot(Ptm,Vc/max(Vrange),'-r','HandleVisibility','off');
% % plot(Pel_range,Vrange/max(Vrange),'-k','HandleVisibility','off');
% % plot(Pel,V/max(Vrange),'-r','HandleVisibility','off');
% % p = plot(Pcw_range(1:4:end),Vrange(1:4:end),'k>','DisplayName','$P_{cw}$, nonlinear compliance');
% p.MarkerFaceColor = 'w';
% if comp_tog == 0
%     plot(Pcw_range+Pel_range,Vrange/max(Vrange),'--k','DisplayName','Constant compliance');
%     plot(Pcw+Pel,V/max(Vrange),'--r','HandleVisibility','off');
% elseif comp_tog == 1 
%     plot(Pcw_range+Pel_range,Vrange/max(Vrange),'-k','DisplayName','Varying compliance');
%     plot(Pcw+Pel,V/max(Vrange),'-r','HandleVisibility','off');
% end 
% legend box off location southeast
% ylabel('$V/V_{TLC}$')
% xlabel('$P_{cw}+P_{el}$')
% stop

% SET TIME STEP------------------------------------------------------------
dt = 10^-3;

% SOLVE ODE SYSTEM=========================================================
%==========================================================================
% The system is solved for each particle size in dp_list
dp_list1 = linspace(.01,.1,10);
dp_list2 = linspace(.2,1,9);
dp_list3 = linspace(2,10,9);
dp_list = [dp_list1 dp_list2 dp_list3]*10^-4; % particle sizes in centimeters

% Heyder1986 = readtable('data.xlsx','Sheet','Heyder1986','VariableNamingRule','preserve');
% x = Heyder1986.PD;
% x = x(find(x<=8));
% dp_list = x*10^-4;

% dp_list = .01*10^-4; % to cm
% dp_list = dp_list(1:2:end);

for i = 1:length(dp_list)
    dp = dp_list(i);
    tvals1 = 0:dt:res.T;
    IC = [0 0.96*vol.FRC vol.FRC 0 0 0 0 0 0 0 0];
    [t1,s1] = ode15s(@(t1,s1) DepositionODEs(t1,s1,par,Y,res,vol,dp,Qin,Qex,Vin,Vex,mech,comp_tog),tvals1,IC,options);
    V = s1(:,3); % L
    Nu = s1(:,4);
    Nc = s1(:,6); % mass in conducting
    Nrb = s1(:,8); % mass in respiratory
    Nalv = s1(:,10);
    Ntotal = Nu+Nc+Nrb+Nalv;
    breaths=1;
    Tind = find(abs((breaths-1)*res.T-t1)==min(abs((breaths-1)*res.T-t1)));

    breaths=2;
    tvals2 = t1(end):dt:breaths*res.T;
    IC = s1(end,:);
    diff = abs(Ntotal(end)/V(end)/par.Cin-Ntotal(Tind)/V(Tind)/par.Cin);
    
%     diff = 100;
    t = t1;
    s = s1;
    while diff>=.01
    [t2,s2] = ode15s(@(t2,s2) DepositionODEs(t2,s2,par,Y,res,vol,dp,Qin,Qex,Vin,Vex,mech,comp_tog),tvals2,IC,options);
    t = [t; t2(2:end)];
    s = [s; s2(2:end,:)];

    V = s(:,3); % L
    Nu = s(:,4);
    Nc = s(:,6); % mass in conducting
    Nrb = s(:,8); % mass in respiratory
    Nalv = s(:,10);

    Ntotal = Nu+Nc+Nrb+Nalv;
    
    Tind = find(abs((breaths-1)*res.T-t)==min(abs((breaths-1)*res.T-t)));
%     diff = abs(s(end,6)/s(end,2)/par.Cin-s(Tind,6)/s(Tind,2)/par.Cin);
%     diff = abs(s(end,7)/s(end,2)/par.Cin-s(Tind,7)/s(Tind,2)/par.Cin);
%     diff = abs(s(end,8)/s(end,3)/par.Cin-s(Tind,8)/s(Tind,2)/par.Cin);
    diff = abs(Ntotal(end)/V(end)/par.Cin-Ntotal(Tind)/V(Tind)/par.Cin);
    
    breaths = breaths+1;
    tvals2 = t2(end):dt:breaths*res.T;
    IC = s2(end,:);%[s2(end,1) s2(end,2) s2(end,3) s2(end,4) s2(end,5) s2(end,6) s2(end,7)];  

    end

    [f(:,i),g(:,i)] = CalculateDeposition(t,s,dp,par,res,Y,vol,breaths-1,Qin,Qex,Vin,Vex);
end

% DEPOSITION PLOTS=========================================================
figure(1)
% total deposition fraction after repeatability
hold on
if comp_tog==0
    p5 = plot(dp_list*10^4,f(5,:),'-k>','MarkerSize',6,'DisplayName','constant compliance');
    p5.MarkerFaceColor='w';
elseif comp_tog==1
    if turb_tog==0
        p5 = plot(dp_list*10^4,f(5,:),'-k>','MarkerSize',6,'DisplayName','turbulence off');
        p5.MarkerFaceColor='k';
    elseif turb_tog==1
        p5 = plot(dp_list*10^4,f(5,:),'-k>','MarkerSize',6,'DisplayName','turbulence on');
        p5.MarkerFaceColor='w';
    end
end
set(gca,'XScale','log')
set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
set(gcf,'Color',[1 1 1])
ylabel('Deposition Fraction, $\eta_{total}$')
xlabel('Particle Size ($\mu$m)')
legend box off location northwest
grid on
ylim([0 1])

figure(2)
% TB deposition after repeatability
hold on
if comp_tog==0
    p3 = plot(dp_list*10^4,f(1,:)+f(2,:),'-ks','MarkerSize',6,'DisplayName','constant compliance');
    p3.MarkerFaceColor='w';
elseif comp_tog==1
    if turb_tog==0
        p3 = plot(dp_list*10^4,f(1,:)+f(2,:),'-ks','MarkerSize',6,'DisplayName','turbulence off');
        p3.MarkerFaceColor='k';
    elseif turb_tog==1
        p3 = plot(dp_list*10^4,f(1,:)+f(2,:),'-ks','MarkerSize',6,'DisplayName','turbulence on');
        p3.MarkerFaceColor='w';
    end
end
set(gca,'XScale','log')
set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
set(gcf,'Color',[1 1 1])
ylabel('Tracheobronchial Deposition')
legend box off location northwest
grid on
ylim([0 1])

figure(3)
% respiratory bronchioles deposition
hold on
if comp_tog==0
    p4 = plot(dp_list*10^4,f(3,:),'-ko','MarkerSize',6,'DisplayName','constant compliance');
    p4.MarkerFaceColor='w';
elseif comp_tog==1
    if turb_tog==0
        p4 = plot(dp_list*10^4,f(3,:),'-ko','MarkerSize',6,'DisplayName','turbulence off');
        p4.MarkerFaceColor='k';
    elseif turb_tog==1
        p4 = plot(dp_list*10^4,f(3,:),'-ko','MarkerSize',6,'DisplayName','turbulence on');
        p4.MarkerFaceColor='w';
    end
end
set(gca,'XScale','log')
set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
set(gcf,'Color',[1 1 1])
ylabel('Respiratory Bronchioles Deposition')
xlabel('Particle Size ($\mu$m)')
legend box off location northwest
grid on
ylim([0 1])

figure(4)
% alveolar deposition
hold on
if comp_tog==0
    p4 = plot(dp_list*10^4,f(4,:),'-kd','MarkerSize',6,'DisplayName','constant compliance');
    p4.MarkerFaceColor='w';
elseif comp_tog==1
    if turb_tog==0
        p4 = plot(dp_list*10^4,f(4,:),'-kd','MarkerSize',6,'DisplayName','turbulence off');
        p4.MarkerFaceColor='k';
    elseif turb_tog==1
        p4 = plot(dp_list*10^4,f(4,:),'-kd','MarkerSize',6,'DisplayName','turbulence on');
        p4.MarkerFaceColor='w';
    end
end
set(gca,'XScale','log')
set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
set(gcf,'Color',[1 1 1])
ylabel('Alveolar Deposition')
xlabel('Particle Size ($\mu$m)')
legend box off location northwest
grid on
ylim([0 1])

figure(5)
% respiratory deposition
hold on
if comp_tog==0
    p4 = plot(dp_list*10^4,f(3,:)+f(4,:),'-kd','MarkerSize',6,'DisplayName','constant compliance');
    p4.MarkerFaceColor='w';
elseif comp_tog==1
    if turb_tog==0
        p4 = plot(dp_list*10^4,f(3,:)+f(4,:),'-kd','MarkerSize',6,'DisplayName','turbulence off');
        p4.MarkerFaceColor='k';
    elseif turb_tog==1
        p4 = plot(dp_list*10^4,f(3,:)+f(4,:),'-kd','MarkerSize',6,'DisplayName','turbulence on');
        p4.MarkerFaceColor='w';
    end
    
end
set(gca,'XScale','log')
set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
set(gcf,'Color',[1 1 1])
ylabel('Pulmonary Deposition')
xlabel('Particle Size ($\mu$m)')
legend box off location northwest
grid on
ylim([0 1])

% RETENTION/EXHALE PLOTS==========================================================
figure(6)
% total exhaled
hold on
if comp_tog==0
    p4 = plot(dp_list*10^4,1-(f(5,:)+g(10,:)),'-k<','MarkerSize',6,'DisplayName','constant compliance');
    p4.MarkerFaceColor='w';
elseif comp_tog==1
    p4 = plot(dp_list*10^4,1-(f(5,:)+g(10,:)),'-k<','MarkerSize',6,'DisplayName','varying compliance');
    p4.MarkerFaceColor='k';
end
set(gca,'XScale','log')
set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
set(gcf,'Color',[1 1 1])
ylabel('Total Exhaled Fraction')
xlabel('Particle Size ($\mu$m)')
legend box off location southwest
grid on
ylim([0 1])

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
if comp_tog==0
    p4 = plot(dp_list*10^4,g(10,:),'-k<','MarkerSize',6,'DisplayName','constant compliance');
    p4.MarkerFaceColor='w';
elseif comp_tog==1
    if turb_tog==0
        p4 = plot(dp_list*10^4,g(10,:),'-k<','MarkerSize',6,'DisplayName','turbulence off');
        p4.MarkerFaceColor='k';
    elseif turb_tog==1
        p4 = plot(dp_list*10^4,g(10,:),'-k<','MarkerSize',6,'DisplayName','turbulence on');
        p4.MarkerFaceColor='w';
    end
end
set(gca,'XScale','log')
set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
set(gcf,'Color',[1 1 1])
ylabel('Total Suspended Fraction')
xlabel('Particle Size ($\mu$m)')
legend box off location southwest
grid on
ylim([0 1])

figure(12)
% TB suspension after repeatability
hold on
if comp_tog==0
    p3 = plot(dp_list*10^4,g(6,:)+g(7,:),'-ks','MarkerSize',6,'DisplayName','constant compliance');
    p3.MarkerFaceColor='w';
elseif comp_tog==1
    if turb_tog==0
        p3 = plot(dp_list*10^4,g(6,:)+g(7,:),'-ks','MarkerSize',6,'DisplayName','turbulence off');
        p3.MarkerFaceColor='k';
    elseif turb_tog==1
        p3 = plot(dp_list*10^4,g(6,:)+g(7,:),'-ks','MarkerSize',6,'DisplayName','turbulence on');
        p3.MarkerFaceColor='w';
    end
end
set(gca,'XScale','log')
set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
set(gcf,'Color',[1 1 1])
ylabel('Tracheobronchial Suspension')
legend box off location northwest
grid on
ylim([0 1])

figure(13)
% respiratory bronchioles suspension
hold on
if comp_tog==0
    p4 = plot(dp_list*10^4,g(8,:),'-ko','MarkerSize',6,'DisplayName','constant compliance');
    p4.MarkerFaceColor='w';
elseif comp_tog==1
    if turb_tog==0
        p4 = plot(dp_list*10^4,g(8,:),'-ko','MarkerSize',6,'DisplayName','turbulence off');
        p4.MarkerFaceColor='k';
    elseif turb_tog==1
        p4 = plot(dp_list*10^4,g(8,:),'-ko','MarkerSize',6,'DisplayName','turbulence on');
        p4.MarkerFaceColor='w';
    end
end
set(gca,'XScale','log')
set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
set(gcf,'Color',[1 1 1])
ylabel('Respiratory Bronchioles Suspension')
xlabel('Particle Size ($\mu$m)')
legend box off location northwest
grid on
ylim([0 1])

figure(14)
% alveolar suspension
hold on
if comp_tog==0
    p4 = plot(dp_list*10^4,g(9,:),'-kd','MarkerSize',6,'DisplayName','constant compliance');
    p4.MarkerFaceColor='w';
elseif comp_tog==1
    if turb_tog==0
        p4 = plot(dp_list*10^4,g(9,:),'-kd','MarkerSize',6,'DisplayName','turbulence off');
        p4.MarkerFaceColor='k';
    elseif turb_tog==1
        p4 = plot(dp_list*10^4,g(9,:),'-kd','MarkerSize',6,'DisplayName','turbulence on');
        p4.MarkerFaceColor='w';
    end
end
set(gca,'XScale','log')
set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
set(gcf,'Color',[1 1 1])
ylabel('Alveolar Suspension')
xlabel('Particle Size ($\mu$m)')
legend box off location northwest
grid on
ylim([0 1])

figure(15)
% respiratory suspension
hold on
if comp_tog==0
    p4 = plot(dp_list*10^4,g(8,:)+g(9,:),'-kd','MarkerSize',6,'DisplayName','constant compliance');
    p4.MarkerFaceColor='w';
elseif comp_tog==1
    if turb_tog==0
        p4 = plot(dp_list*10^4,g(8,:)+g(9,:),'-kd','MarkerSize',6,'DisplayName','turbulence on');
        p4.MarkerFaceColor='k';
    elseif turb_tog==1
        p4 = plot(dp_list*10^4,g(8,:)+g(9,:),'-kd','MarkerSize',6,'DisplayName','turbulence off');
        p4.MarkerFaceColor='w';
    end
end
set(gca,'XScale','log')
set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
set(gcf,'Color',[1 1 1])
ylabel('Pulmonary Suspension')
xlabel('Particle Size ($\mu$m)')
legend box off location northwest
grid on
ylim([0 1])



% stop
% figure(7)
% % deposition at end inhale
% hold on
% p = plot(dp_list*10^4,f(6,:)+f(7,:),'-ks','MarkerSize',6,'DisplayName','u+c');
% p.MarkerFaceColor='w';
% p = plot(dp_list*10^4,f(8,:)+f(9,:),'-ko','MarkerSize',6,'DisplayName','rb+alv');
% p.MarkerFaceColor='w';
% p = plot(dp_list*10^4,f(10,:),'-k>','MarkerSize',6,'DisplayName','total');
% p.MarkerFaceColor='w';
% set(gca,'XScale','log')
% set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
% set(gcf,'Color',[1 1 1])
% ylabel('$\eta_{i,inhale}$')
% xlabel('Particle Size ($\mu$m)')
% legend box off location southwest
% grid on
% ylim([0 1])
% 
% figure(8)
% % deposition during exhale
% hold on
% p = plot(dp_list*10^4,f(11,:)+f(12,:),'-ks','MarkerSize',6,'DisplayName','u+c');
% p.MarkerFaceColor='w';
% p = plot(dp_list*10^4,f(13,:)+f(14,:),'-ko','MarkerSize',6,'DisplayName','rb+alv');
% p.MarkerFaceColor='w';
% p = plot(dp_list*10^4,f(15,:),'-k>','MarkerSize',6,'DisplayName','total');
% p.MarkerFaceColor='w';
% set(gca,'XScale','log')
% set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
% set(gcf,'Color',[1 1 1])
% ylabel('$\eta_{i,exhale}$')
% xlabel('Particle Size ($\mu$m)')
% legend box off location southwest
% grid on
% ylim([0 1])
% 
% figure(9)
% % suspension at end inhale
% hold on
% p = plot(dp_list*10^4,g(1,:)+g(2,:),'-ks','MarkerSize',6,'DisplayName','u+c');
% p.MarkerFaceColor='w';
% p = plot(dp_list*10^4,g(3,:)+g(4,:),'-ko','MarkerSize',6,'DisplayName','rb+alv');
% p.MarkerFaceColor='w';
% p = plot(dp_list*10^4,g(5,:),'-k>','MarkerSize',6,'DisplayName','total');
% p.MarkerFaceColor='w';
% set(gca,'XScale','log')
% set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
% set(gcf,'Color',[1 1 1])
% ylabel('$\alpha_{i,inhale}$')
% xlabel('Particle Size ($\mu$m)')
% legend box off location southwest
% grid on
% ylim([0 1])
% 
% figure(10)
% % suspension at end exhale
% hold on
% p = plot(dp_list*10^4,g(6,:)+g(7,:),'-ks','MarkerSize',6,'DisplayName','u+c');
% p.MarkerFaceColor='w';
% p = plot(dp_list*10^4,g(8,:)+g(9,:),'-ko','MarkerSize',6,'DisplayName','rb+alv');
% p.MarkerFaceColor='w';
% p = plot(dp_list*10^4,g(10,:),'-k>','MarkerSize',6,'DisplayName','total');
% p.MarkerFaceColor='w';
% set(gca,'XScale','log')
% set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
% set(gcf,'Color',[1 1 1])
% ylabel('$\alpha_{i,exhale}$')
% xlabel('Particle Size ($\mu$m)')
% legend box off location southwest
% grid on
% ylim([0 1])
% 
% figure(11)
% % efficienty of exhale in TB
% hold on
% p = plot(dp_list*10^4,f(6,:)+f(7,:),'-ks','MarkerSize',6,'DisplayName','$\eta_{u,inhale}+\eta_{c,inhale}$');
% p.MarkerFaceColor='w';
% p = plot(dp_list*10^4,f(11,:)+f(12,:),'--ks','MarkerSize',6,'DisplayName','$\eta_{u,exhale}+\eta_{c,exhale}$');
% p.MarkerFaceColor='w';
% p = plot(dp_list*10^4,f(16,:)+f(17,:),'-.ks','MarkerSize',6,'DisplayName','$\eta_{u,efficiency}+\eta_{c,efficiency}$');
% p.MarkerFaceColor='w';
% set(gca,'XScale','log')
% set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
% set(gcf,'Color',[1 1 1])
% ylabel('Deposition')
% xlabel('Particle Size ($\mu$m)')
% legend box off location southwest
% grid on
% ylim([0 1])
% 
% figure(12)
% % efficienty of exhale in resp bronch
% hold on
% p = plot(dp_list*10^4,f(8,:),'-ko','MarkerSize',6,'DisplayName','$\eta_{rb,inhale}$');
% p.MarkerFaceColor='w';
% p = plot(dp_list*10^4,f(13,:),'--ko','MarkerSize',6,'DisplayName','$\eta_{rb,exhale}$');
% p.MarkerFaceColor='w';
% p = plot(dp_list*10^4,f(18,:),'-.ko','MarkerSize',6,'DisplayName','$\eta_{rb,efficiency}$');
% p.MarkerFaceColor='w';
% set(gca,'XScale','log')
% set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
% set(gcf,'Color',[1 1 1])
% ylabel('Deposition')
% xlabel('Particle Size ($\mu$m)')
% legend box off location southwest
% grid on
% ylim([0 1])
% 
% figure(13)
% % efficienty of exhale in aveoli
% hold on
% p = plot(dp_list*10^4,f(9,:),'-ko','MarkerSize',6,'DisplayName','$\eta_{alv,inhale}$');
% p.MarkerFaceColor='w';
% p = plot(dp_list*10^4,f(14,:),'--ko','MarkerSize',6,'DisplayName','$\eta_{alv,exhale}$');
% p.MarkerFaceColor='w';
% p = plot(dp_list*10^4,f(19,:),'-.ko','MarkerSize',6,'DisplayName','$\eta_{alv,efficiency}$');
% p.MarkerFaceColor='w';
% set(gca,'XScale','log')
% set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
% set(gcf,'Color',[1 1 1])
% ylabel('Deposition')
% xlabel('Particle Size ($\mu$m)')
% legend box off location southwest
% grid on
% ylim([0 1])
% 
% figure(14)
% % efficienty of exhale total
% hold on
% p = plot(dp_list*10^4,f(10,:),'-k>','MarkerSize',6,'DisplayName','$\eta_{total,inhale}$');
% p.MarkerFaceColor='w';
% p = plot(dp_list*10^4,f(15,:),'--k>','MarkerSize',6,'DisplayName','$\eta_{total,exhale}$');
% p.MarkerFaceColor='w';
% p = plot(dp_list*10^4,f(20,:),'-.k>','MarkerSize',6,'DisplayName','$\eta_{total,efficiency}$');
% p.MarkerFaceColor='w';
% set(gca,'XScale','log')
% set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
% set(gcf,'Color',[1 1 1])
% ylabel('Deposition')
% xlabel('Particle Size ($\mu$m)')
% legend box off location southwest
% grid on
% ylim([0 1])
end

