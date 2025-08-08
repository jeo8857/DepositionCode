function [Ru,Rc,Rra] = Resistance(par,vol,Ls,rs,n,c1_ind,c2_ind,c3_ind,c4_ind,Q,VA,Vc,Vcmax)
% Expressions need to be changed manually in ConstantFlow.m
%==========================================================================
% This function returns the compartmental resistances.
%==========================================================================

% Poiseuille resistance
R_pipe = (8/pi)*par.mu*(Ls./rs.^4)*(100/98.1); 
R_gen = R_pipe./n;

% Ru = par.Au+par.Ku*abs(Q); % cmh2o*s/L
R_lam = sum(R_gen(c1_ind)); % laminar resistance
Ru = R_lam+par.Ku*abs(Q); % cmh2o*s/L laminar resistance + turbulant resistance

% Rc = par.Kc*(Vcmax/Vc)^2;
Rc = sum(R_gen(c2_ind));
Rra = sum(R_gen([c3_ind c4_ind]));
end