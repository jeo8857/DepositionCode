function sys = LungMechanicsODEs(t,s,par,Y,res,vol,comp_tog,IC)
%==========================================================================
% This function outputs the solution to the breathing mechanics model at
% time t
%==========================================================================

% FORCING FUNCTION
% Pmus = pressure induced by respiratory muscles
% This function drives the system
% Taken from Albanese et al (2016)
tt = mod(t,res.T);
if tt>=0 && tt<=res.TI
    Pmus = -res.Pmusmin*tt^2/(res.TI*res.TE)+res.Pmusmin*res.T*tt/(res.TI*res.TE);
elseif tt>res.TI && tt<=res.T
    Pmus = res.Pmusmin/(1-exp(-res.TE/res.tau))*(exp(-(tt-res.TI)./res.tau)-exp(-res.TE/res.tau));
end

Q = s(1);
VA = s(2);
V = s(3);

% Q = s(1);
% Vra = s(2);
% Pve =s(3);
% V = s(4);

% SCALE GEOMETRIC DATA TO DESIRED LUNG VOLUME (Currently scaled to total
% lung capacity)
[rs,Ls,n,Vc,Vcmax,c1_ind,c2_ind,c3_ind,c4_ind,~,~,~,~,~,~,~,~,Vu] = Scaling(Y,vol,V,VA);
% [rs,Ls,n,Vc,Vcmax,c1_ind,c2_ind,c3_ind,c4_ind,~,~,~,~,~,~,~] = Scaling(Y,vol,V,VA);

% GET PV CURVES
[Pcw,Ptm,Pel] = PV(comp_tog,par,vol,V,Vc,VA,Vcmax,Vu,IC);

% COMPUTE THE RESISTANCES
[Ru,Rc,Rra] = Resistance(par,vol,Ls,rs,n,c1_ind,c2_ind,c3_ind,c4_ind,Q,VA,Vc,Vcmax);

% STATE EQUATIONS
sys(1) = 1/par.Iu*(par.Pao-Pmus-(Ru+Rc)*Q-Ptm-Pcw); % Q
sys(2) = (Ptm-Pel)/Rra; % Va
sys(3) = Q; % V

sys = sys';

end


