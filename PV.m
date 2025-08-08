function [Pcw,Ptm,Pel] = PV(tog,par,vol,V,Vc,VA,Vcmax,Vu,IC)
%==========================================================================
% This function returns the pressure-volume (PV) curves. If tog==0, then it
% returns the linear PV curves (constant compliance). If tog==1, then it
% returns the nonlinear PV curves (varying compliance).
%==========================================================================

if tog==0
    V0 = vol.FRC-par.Ccw*par.Ppl;
    Pcw = (1/par.Ccw)*(V-V0);
   
    Vc0 = (vol.FRC-IC(2))-Vu; % conducting volume at rest (set by chosen initial conditions)
    
    Vuc = Vc0+par.Cc*par.Ppl;
    Ptm = (Vc-Vuc)/par.Cc;
    
    VA0 = IC(2); % respiratory volume at rest (set by the chosen initial conditions)
    VuA = VA0+par.Cra*par.Ppl;
    Pel = (VA-VuA)/par.Cra;
    
elseif tog==1
    Pcw = par.Acw-par.Bcw*log((vol.TLC-vol.RV)./(V-vol.RV)-0.999);
    Ptm = par.Ac-par.Bc*log(Vcmax./Vc-0.999);
    Pel = par.Al*exp(par.Kl*VA)+par.Bl;
    
end

end