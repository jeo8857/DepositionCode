function par = Parameters(volumes,Y,turb_tog,IC)

vol = volumes;

% % MORPHOMETRY
% % save each column as array
% gnum = Y{:,1}; % generation number
% n = Y{:,2}; % number of airways in each generation
% Ldata = Y{:,3}; % length of airways in each generation in cm
% rdata = Y{:,4}/2; % radius of airways in each generation in cm
% % Ba = Y{:,5}.*pi/180; % branching angle converted FROM degrees TO radians
% % Ga = Y{:,6}.*pi/180; % gravity angle converted FROM degrees TO radians
% Vdata = Y{:,8}/1e3; % volume in L

% Indices for each compartment
% c1_ind = 1:5; 
% c1_ind = 1:2; 
% % c2_ind = 6:16;
% c2_ind = 3:16;
% c3_ind = 17:length(gnum);

% Scale to our desired lung volume to initialize
% sc = (vol.TLC/sum(Vdata))^(1/3); % scaled to TLC
% sc = (4.8/sum(Vdata))^(1/3); % scaled to TLC
% r0 = rdata*sc; % cm
% L0 = Ldata*sc; % cm

% Vu = sum(pi*r0(c1_ind).^2.*L0(c1_ind).*n(c1_ind))/1e3
% Vcmax = sum(pi*r0(c2_ind).^2.*L0(c2_ind).*n(c2_ind))/1e3 % L
% stop

[~,~,~,~,Vcmax,~,~,~,~,~,~,~,~,~,~,~,~,Vu] = Scaling(Y,vol,vol.TLC,IC(2));

% subjects are taken from Athanasiades(2000)
%--------------------------------------------------------------------------
par.Cin = 3*10^8;%.2; % concentration determined by device, mass/volume

par.Pao = 0;

% Lung compliance at rest
par.Ccw = .2445;
par.Cc = .0131;%.0131;
par.Cra = .208;

% Pleural pressure at rest
% par.Ppl = -8.8;
par.Ppl = -(vol.FRC-vol.RV)/(par.Cc+par.Cra);

% CHEST WALL
% par.Acw = 0.65; % I picked this value 
% par.Acw = 1.4;% subjects 1,2,3
% par.Acw = 4.4; % subject 4
% par.Bcw = 2.1; % (I adjusted this value to get desired compliance)
% par.Bcw = 3.5; % same for all
% par.Bcw = 2.5; % I picked
par.Bcw = (vol.FRC-vol.RV)*(vol.TLC-vol.FRC)/(par.Ccw*(vol.TLC-vol.RV));
par.Acw = par.Ppl+par.Bcw*log((vol.TLC-vol.RV)/(vol.FRC-vol.RV)-1);

% RIGID AIRWAYS
% par.Au = 0.34; % subject 1
par.Au = 0.31; % subjects 2,3,4

if turb_tog==0
    par.Ku = 0;
    par.Iu = .0000001;%.0001;
elseif turb_tog==1
    par.Ku = 0.32; % subject 3
%     par.Iu = 0.33;
    par.Iu = 0.05;
end
% par.Ku = 0.46; % subject 1
% par.Ku = 0.4; % subject 2
% par.Ku = 0.2; % subject 4


% par.Iu = 0.1;
 


% CONDUCTING AIRWAYS
% par.Kc = 0.21; % subject 1
% par.Kc = 0.49; % subject 2
par.Kc = 0.5; % subject 3
% par.Kc = 0.24; % subject 4

% par.Ac = 1.96;%1.1588;%1;%7.09;
% par.Ac = 0.7;%1;%7.09;
% par.Bc2 = 5.5;%3.73;

% par.Ac = 1;%10;
% par.Bc2 = 1;%3.73;

% Vcmax = 0.2;
% par.Bc2 = (0.04*vol.FRC)*(Vcmax-(0.04*vol.FRC))/(par.Cc*Vcmax);
Vc0 = (vol.FRC-IC(2))-Vu; % Volume of conducting airways at rest
par.Bc = Vc0*(Vcmax-Vc0)/(par.Cc*Vcmax);

% par.Ac = -par.Ppl+par.Bc2*log(Vcmax/(0.04*vol.FRC)-1);
par.Ac = -par.Ppl+par.Bc*log(Vcmax/Vc0-1);

% RESPIRATORY AIRWAYS
par.Kl = 1; % subjects 1,2,3
% par.Kl = 0.8; % subject 4

% par.Al = 0.2; % subject 1
% par.Al = 0.04; % subject 2
% par.Al = 0.1; % subject 3
% par.Al = 0.57; % subject 4
% par.Al = 0.19;%0.126265;%0.1235;%0.1;%0.2;
par.Al = exp(-par.Kl*IC(2))/(par.Cra*par.Kl);


% par.Bl = -0.5; % subject 1
% par.Bl = 1; % subject 2
% par.Bl = 1.5; % subject 3
% par.Bl = 0; % subject 4
% par.Bl = -2;
par.Bl = -par.Ppl-par.Al*exp(par.Kl*IC(2));

% par.As = 2.2; % subject 1
% par.As = 2.8; % subject 2
par.As = 2.47; % subject 3
% par.As = 5.47; % subject 4

% par.Ks = -10.9; % subject 1
% par.Ks = -9.9; % subject 2
par.Ks = -6.5; % subject 3
% par.Ks = -5.13; % subject 4

par.Bs = 0.02; % same for all subjects

par.Cve = 0.5;
par.Rve = 1;

% Air properties
% par.mu = 1.86*10^(-5); % viscosity of air at 25 degrees celcius in P*s (pascal = Newtons/m^2 = kg/(s^2*m))
par.mu = 1.84*10^(-4); % viscosity of air g cm^-1 s^-1 (Martonen, 1993)
par.g = 980; % gravitational acceleration cm/s^2
par.temp = 300; % normal adult body temperature in kelvins

% Particle properties
par.rho = 1; % particle density in g/cm^3
% par.rho = 1000; % particle density in kg/m^3
par.k = 1.381*1e-16;%*10^(-23); % Boltzmann constant (Oakes,Martonen)


end
