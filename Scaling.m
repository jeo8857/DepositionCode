function [rs,Ls,n,Vc,Vcmax,c1_ind,c2_ind,c3_ind,c4_ind,Vrb,Valv,VDTLC,VRTLC,r0,L0,Ba,Ga,Vu] = Scaling(Y,vol,V,VA)
% NOTE: Any changes to this function must also be duplicated in ConstantFlow.m
%==========================================================================
% This function outputs the scaled geometries and the volumes computed
% using the scaled geometries
%==========================================================================

% MORPHOMETRIC DATA
% save each column as array
g = Y{:,1}; % generation number
n = Y{:,2}; % number of airways in each generation
Ldata = Y{:,3}; % length of airways in each generation in cm
rdata = Y{:,4}/2; % radius of airways in each generation in cm
Ba = Y{:,5}.*pi/180; % branching angle converted FROM degrees TO radians
Ga = Y{:,6}.*pi/180; % gravity angle converted FROM degrees TO radians
Vdata = Y{:,8}/1e3; % volume in L

% SPECIFY HOW GENERATIONS ARE LUMPED
c1_ind = 1:2; % first compartment (upper)
c2_ind = 3:16; % second compartment (conducting)
c3_ind = 17:length(g)-1; % third compartment (respiratory bronchioles)
c4_ind = length(g); % fourth compartment (alveoli)

% Scale to our desired lung volume to initialize
sf = (vol.TLC/sum(Vdata))^(1/3); % scaled to total lung capacity
% sf = (vol.FRC/sum(Vdata))^(1/3); % scaled to functional residual capacity
r0 = rdata*sf; % initial radii in centimeters
L0 = Ldata*sf; % initial lengths in centimeters

% r_TV = rdata*((vol.FRC+vol.TV)/sum(Vdata))^(1/3);
% L_TV = Ldata*((vol.FRC+vol.TV)/sum(Vdata))^(1/3);

% Constant volume of the upper airways 
Vu = sum(pi*r0(c1_ind).^2.*L0(c1_ind).*n(c1_ind))/1e3; % L

% Volume of the conducting airways at time t
Vc = V-(VA+Vu);

% Volume of the conducting airways at TLC
Vcmax = sum(pi*r0(c2_ind).^2.*L0(c2_ind).*n(c2_ind))/1e3; % L

% Volume of the deadspace (upper airways + conducting) at TLC
VDTLC = Vu+Vcmax;

% Volume of the respiratory airways at TLC
VRTLC = vol.TLC-VDTLC;

% Scale factors
cscale = ((Vu+Vc)/VDTLC)^(1/2); % conducting airways scale factor
rscale = (VA/VRTLC)^(1/3); % respiratory airways scale factor

% Scale the airways by volume (dimensions in cm)
rs(c1_ind) = r0(c1_ind); % Doesn't get scaled because no compliance
rs(c2_ind) = cscale*r0(c2_ind);
rs(c3_ind) = rscale*r0(c3_ind);
rs(c4_ind) = rscale*r0(c4_ind);

Ls(c1_ind) = L0(c1_ind);
Ls(c2_ind) = L0(c2_ind);
Ls(c3_ind) = rscale*L0(c3_ind);
Ls(c4_ind) = rscale*L0(c4_ind);

rs = rs(:);
Ls = Ls(:);

% Compute scaled volume of the respiratory bronchioles
Vrb = sum(pi*rs(c3_ind).^2.*Ls(c3_ind).*n(c3_ind))/1e3;

% Compute volume of alveoli
Valv = VA-Vrb;
end