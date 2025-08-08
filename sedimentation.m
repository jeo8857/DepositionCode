function S = sedimentation(n,r,L,phi,par,dp,Q)
%==========================================================================
% This function generates the likelihood of gravitational sedimentation in each
% generation.
%==========================================================================

r = r(:);
L = L(:);
v = 1e3*abs(Q)./(n.*(pi*r.^2)); %cm/s
% v = 1e3*abs(Q)./((pi*r.^2)); % cm^3/s - velocity of air through each individual airway in each generation, m/s

tau = par.rho*dp^2/(18*par.mu); % relaxation time
vs = par.g*tau; % settling/terminal velocity
t_res = L./v; % residence time in an individual airway for each generation

% cunningham
lam = 7*10^-6; %cm (martonen); % @ 37 deg Cel, 100% humidity and 76 cmhg atmospheric pressure (icrp 1994)
c = 1+(lam/dp)*(2.514+0.8*exp(-0.55*(dp/lam)));
% Re = 0.001.*v.*L./par.mu;
% Re>2900
% stop
S = 1-exp(-(2/pi)*cos(phi).*c.*vs.*t_res./r);

end

