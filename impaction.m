function I = impaction(r,theta,n,par,dp,Q)
%==========================================================================
% This function generates the likelihood of inertial impaction in each
% generation.
%==========================================================================

v = 1e3*abs(Q)./(n.*(pi*r.^2)); % cm^3/s - velocity of air through each individual airway in each generation, m/s
% v = 1e3*abs(Q)./((pi*r.^2)); % cm^3/s - velocity of air through each individual airway in each generation, m/s

tau = par.rho*dp^2/(18*par.mu); % s- relaxation time
Sp = v.*tau; % cm - stopping distance of particle in any individual airway for each generation, m

lam = 7*10^-6; %cm (martonen); 0.0712*10^-6; % @ 37 deg Cel, 100% humidity and 76 cmhg atmospheric pressure (icrp 1994)
c = 1+(lam/dp)*(2.514+0.8*exp(-0.55*(dp/lam)));

stk = c.*Sp./(2*r); % stokes number for an individual airway in each generation

% for i = 1:length(n)
%     if stk(i)*theta(i) < 1
%         I(i) = 1-(2/pi)*acos(theta(i)*stk(i))+(1/pi)*sin(2*acos(theta(i)*stk(i)));
%     else
%         I(i) = 1;
%     end
% end

I = 1-exp(-4*theta.*stk);

end

