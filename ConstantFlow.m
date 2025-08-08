function [Qin,Qex,Vin,Vex] = ConstantFlow(t,s,Y,vol,par,res,comp_tog,IC)
% NOTE: Need to be sure that the values for scaling and resistances are
% consistant with those used in Scaling.m and Resistance.m. The reason I
% didn't simply make calls to these functions is because they return single
% values. For the following computations, arrays are needed. I got too lazy
% to fix this annoyance.

%==========================================================================
% This function returns the constant flows and volumes in each compartment,
% which are needed to compute the deposition probabilities. 
%==========================================================================

% MORPHOMETRY
% save each column as array
g = Y{:,1}; % generation number
n = Y{:,2}; % number of airways in each generation
Ldata = Y{:,3}; % length of airways in each generation in cm
rdata = Y{:,4}/2; % radius of airways in each generation in cm
% Ba = Y{:,5}.*pi/180; % branching angle converted FROM degrees TO radians
% Ga = Y{:,6}.*pi/180; % gravity angle converted FROM degrees TO radians
Vdata = Y{:,8}/1e3; % volume in L

% Indices for each compartment (if groupings are changed, they also need to
% be changed in Scaling.m)
c1_ind = 1:2;
c2_ind = 3:16;  
c3_ind = 17:length(g)-1; 
c4_ind = length(g);

% Scale to our desired lung volume to initialize (if choice of scaling
% changes, it must also be changed in Scaling.m)
sf = (vol.TLC/sum(Vdata))^(1/3); % scaled to TLC
% sf = (vol.FRC/sum(Vdata))^(1/3); % scaled to FRC
% sc = (4.8/sum(Vdata))^(1/3); % scaled to TLC
r0 = rdata*sf; % cm
L0 = Ldata*sf; % cm

Q = s(:,1);
VA = s(:,2);
V = s(:,3);

% Q = s(1);
% Vra = s(2);
% Pve =s(3);
% V = s(4);

% r_TV = rdata*((vol.FRC+vol.TV)/sum(Vdata))^(1/3);
% L_TV = Ldata*((vol.FRC+vol.TV)/sum(Vdata))^(1/3);

% Constant volume of the upper airways 
Vu = sum(pi*r0(c1_ind).^2.*L0(c1_ind).*n(c1_ind))/1e3; % L

% Volume of the conducting airways at time t
Vc = V-(VA+Vu);

% Volume of the conducting airways at TLC
Vcmax = sum(pi*r0(c2_ind).^2.*L0(c2_ind).*n(c2_ind))/1e3; % L

% VDratio = sum(Vdata([c1_ind c2_ind]))/sum(Vdata);
% 
% % VDTLC = VDratio*vol.TLC; % L
% VRTLC = (1-VDratio)*vol.TLC; % L

% Volume of the deadspace (upper airways + conducting) at TLC
VDTLC = Vu+Vcmax;

% Volume of the respiratory airways at TLC
VRTLC = vol.TLC-VDTLC;

% cscale = ((Vu+Vc)/VDTLC).^(1/2); % updated
rscale = (VA/VRTLC).^(1/3); 

% cscale = ((Vu+Vc)/(0.04*vol.FRC)).^(1/2); % updated
% rscale = (VA/(0.96*vol.FRC)).^(1/3); % updated

% [~,~,n,~,Vcmax,c1_ind,c2_ind,c3_ind,c4_ind,~,~,VDTLC,VRTLC,r0,L0,Ba,Ga,Vu] = Scaling(Y,vol,V,VA);
% Vc = V-(VA+Vu);

% Initialize scaled geometry array to geometry at TLC
rs = r0;
Ls = L0;
rs = rs(:);
Ls = Ls(:);

for i=1:length(t)
rs(c3_ind) = rscale(i)*r0(c3_ind);
rs(c4_ind) = rscale(i)*r0(c4_ind);
Ls(c3_ind) = rscale(i)*L0(c3_ind);
Ls(c4_ind) = rscale(i)*L0(c4_ind);

% Respiratory resistance has to be recomputed since respiratory flow is not
% a state variable (if choice of resistances changes, must also be changed
% in Resistance.m)
R_pipe = (8/pi)*par.mu*(Ls([c3_ind c4_ind])./rs([c3_ind c4_ind]).^4)*(100/98.1);
R_gen = R_pipe./n([c3_ind c4_ind]);
Rra(i) = sum(R_gen);

Vrb(i) = sum(pi*rs(c3_ind).^2.*Ls(c3_ind).*n(c3_ind))/1e3;

end
% Rra = par.As*exp(par.Ks*(VA-vol.RV)./(vol.TLC-vol.RV))+par.Bs;
Rra = Rra(:);
Vrb = Vrb(:);

Valv = VA-Vrb;

% Switch function u=1 during inhale and u=0 during exhale
tt = mod(t,res.T);
for i = 1:length(tt)
    if tt(i)>=0 && tt(i)<=res.TI
        u(i) = 1;
    elseif tt(i)>res.TI && tt(i)<=res.T
        u(i) = 0;
    end
end
u = u(:);

[Pcw,Ptm,Pel] = PV(comp_tog,par,vol,V,Vc,VA,Vcmax,Vu,IC);

Qra = (Ptm-Pel)./Rra;
Qra = Qra(:);

% AVERAGED FLOWS

% inhale
avg_Qin = (1/res.TI)*trapz(t,u.*Q);
max_Qin = max(u.*Q);
Qcin = (avg_Qin+max_Qin)/2;

avg_Qrain = (1/res.TI)*trapz(t,u.*Qra);
max_Qrain = max(u.*Qra);
Qrain = (avg_Qrain+max_Qrain)/2;
% Qrain = avg_Qrain;

% exhale
avg_Qex = (1/res.TE)*trapz(t,(1-u).*abs(Q));
max_Qex = max((1-u).*abs(Q));
Qcex = (avg_Qex+max_Qex)/2;

avg_Qraex = (1/res.TE)*trapz(t,(1-u).*abs(Qra));
max_Qraex = max((1-u).*abs(Qra));
Qraex = (avg_Qraex+max_Qraex)/2;
% Qraex = avg_Qraex;

Qin = [Qcin,Qrain];
Qex = [Qcex,Qraex];

% VOLUMES

% inhale
% avg_Vcin = (1/res.TI)*trapz(t,u'.*(Vu+Vc));
% max_Vcin = max(u'.*(Vu+Vc));
% Vcin = (avg_Vcin+max_Vcin)/2;
% 
% avg_Vrain = (1/res.TI)*trapz(t,u'.*VA);
% max_Vrain = max(u'.*VA);
% Vrain = (avg_Vrain+max_Vrain)/2;
% 
% avg_Vin = (1/res.TI)*trapz(t,u'.*V);
% max_Vin = max(u'.*V);
% Vin = (avg_Vin+max_Vin)/2;

avg_Vuin = (1/res.TI)*trapz(t,u.*Vu);
max_Vuin = max(u.*Vu);
Vuin = (avg_Vuin+max_Vuin)/2;

avg_Vcin = (1/res.TI)*trapz(t,u.*Vc);
max_Vcin = max(u.*Vc);
Vcin = (avg_Vcin+max_Vcin)/2;

avg_Vrbin = (1/res.TI)*trapz(t,u.*Vrb);
max_Vrbin = max(u.*Vrb);
Vrbin = (avg_Vrbin+max_Vrbin)/2;

avg_Valvin = (1/res.TI)*trapz(t,u.*Valv);
max_Valvin = max(u.*Valv);
Valvin = (avg_Valvin+max_Valvin)/2;

avg_Vin = (1/res.TI)*trapz(t,u.*V);
max_Vin = max(u.*V);
Vtotin = (avg_Vin+max_Vin)/2;

% exhale
% avg_Vex = (1/res.TE)*trapz(t,(1-u)'.*V);
% max_Vex = max((1-u)'.*V);
% Vex = (avg_Vex+max_Vex)/2;
% 
% avg_Vcex = (1/res.TE)*trapz(t,(1-u)'.*(Vu+Vc));
% max_Vcex = max((1-u)'.*(Vu+Vc));
% Vcex = (avg_Vcex+max_Vcex)/2;
% 
% avg_Vraex = (1/res.TE)*trapz(t,(1-u)'.*VA);
% max_Vraex = max((1-u)'.*VA);
% Vraex = (avg_Vraex+max_Vraex)/2;

avg_Vuex = (1/res.TE)*trapz(t,(1-u).*Vu);
max_Vuex = max((1-u).*Vu);
Vuex = (avg_Vuex+max_Vuex)/2;

avg_Vcex = (1/res.TE)*trapz(t,(1-u).*Vc);
max_Vcex = max((1-u).*Vc);
Vcex = (avg_Vcex+max_Vcex)/2;

avg_Vrbex = (1/res.TE)*trapz(t,(1-u).*Vrb);
max_Vrbex = max((1-u).*Vrb);
Vrbex = (avg_Vrbex+max_Vrbex)/2;

avg_Valvex = (1/res.TE)*trapz(t,(1-u).*Valv);
max_Valvex = max((1-u).*Valv);
Valvex = (avg_Valvex+max_Valvex)/2;

avg_Vex = (1/res.TE)*trapz(t,(1-u).*V);
max_Vex = max((1-u).*V);
Vtotex = (avg_Vex+max_Vex)/2;

Vin = [Vuin,Vcin,Vrbin,Valvin,Vtotin];
Vex = [Vuex,Vcex, Vrbex,Valvex,Vtotex];

end