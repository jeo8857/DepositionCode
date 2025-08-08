function D = diffusion(n,r,L,par,dp,Q)
%==========================================================================
% This function generates the likelihood of Brownian diffusion in each
% generation.
%==========================================================================

r  = r(:);
L = L(:);
v = 1e3*abs(Q)./(n.*(pi*r.^2));
% v = 1e3*abs(Q)./((pi*r.^2)); % cm^3/s - velocity of air through each individual airway in each generation, m/s

t_res = L./v; % residence time

% cunningham
lam = 7*10^-6; %cm (martonen); % @ 37 deg Cel, 100% humidity and 76 cmhg atmospheric pressure (icrp 1994)
c = 1+(lam/dp)*(2.514+0.8*exp(-0.55*(dp/lam)));

Dc = c.*par.k*par.temp/(3*pi*par.mu*dp); % diffusion coefficient
Sd = 2*sqrt(Dc.*t_res); % diffusion length
h = Sd.^2./(2*(2*r).^2); % ratio of diffusion length to airway diameter

% D = ones(length(r),1);

% laminar
D = 1-0.819*exp(-7.315*h)-0.0976*exp(-44.61*h)-0.0325*exp(-114*h)-0.0509*exp(-79.31*h.^(2/3));% yeh 1980

% turbulent
% D(1:5) = 2.828*h(1:5).^(1/2).*(1-0.314*h(1:5).^(1/2));


end

