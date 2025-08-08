function sys = DepositionODEs(t,s,par,Y,res,vol,dp,Qin,Qex,Vin,Vex,mech,comp_tog)
%==========================================================================
% This functions solves the deposition equations at each time.
%==========================================================================

% STATE VARIABLES----------------------------------------------------------
Q = s(1); % L/s
% Vra = s(2); % L
VA = s(2); % L
V = s(3); % L
Nu = s(4);
Du = s(5);
Nc = s(6); % mass in conducting
Dc = s(7); % deposited in conducting
Nrb = s(8); % mass in respiratory
Drb = s(9); % deposited in respiratory
Nalv = s(10);
Dalv = s(11);

% SCALE GEOMETRIES TO DESIRED LUNG VOLUME
[rs,Ls,n,Vc,Vcmax,c1_ind,c2_ind,c3_ind,c4_ind,Vrb,Valv,VDTLC,VRTLC,r0,L0,Ba,Ga,Vu] = Scaling(Y,vol,V,VA);

% BREATHING MECHANICS------------------------------------------------------
% PRESSURE INDUCED BY MUSCLES
tt = mod(t,res.T);

if tt>=0 && tt<=res.TI
    Pmus = -res.Pmusmin*tt^2/(res.TI*res.TE)+res.Pmusmin*res.T*tt/(res.TI*res.TE);
    u = 1;
elseif tt>res.TI && tt<=res.T
    Pmus = res.Pmusmin/(1-exp(-res.TE/res.tau))*(exp(-(tt-res.TI)./res.tau)-exp(-res.TE/res.tau));
    u = 0;
end

% GET PV CURVES
[Pcw,Ptm,Pel] = PV(comp_tog,par,vol,V,Vc,VA,Vcmax,Vu,Pmus);

% COMPUTE RESISTANCES
[Ru,Rc,Rra] = Resistance(par,vol,Ls,rs,n,c1_ind,c2_ind,c3_ind,c4_ind,Q,VA,Vc,Vcmax);

% PROBABILITIES------------------------------------------------------------

% During inhale
if tt>=0 && tt<=res.TI
    Qvec(c1_ind) = Qin(1);
    Qvec(c2_ind) = Qin(1);
    Qvec(c3_ind) = Qin(2);
    Qvec(c4_ind) = Qin(2);
    Qvec = Qvec(:);
    
    cscale = ((Vin(1)+Vin(2))/VDTLC)^(1/2); % conducting scale factor
    rscale = ((Vin(3)+Vin(4))/VRTLC)^(1/3); % respiratory scale factor
    
elseif tt>res.TI && tt<=res.T
    Qvec(c1_ind) = Qex(1);
    Qvec(c2_ind) = Qex(1);
    Qvec(c3_ind) = Qex(2);
    Qvec(c4_ind) = Qex(2);
    Qvec = Qvec(:);

    cscale = ((Vex(1)+Vex(2))/VDTLC)^(1/2); % conducting scale factor
    rscale = ((Vex(3)+Vex(4))/VRTLC)^(1/3); % respiratory scale factor
    
end

% Geometries based on averaged Q and V used to compute probabilities
rc = ones(size(r0));
Lc = ones(size(L0));

% group generations into compartments
rc(c1_ind) = r0(c1_ind); % cm
Lc(c1_ind) = L0(c1_ind); % cm

rc(c2_ind) = r0(c2_ind)*cscale; % cm
Lc(c2_ind) = L0(c2_ind); % cm

rc(c3_ind) = r0(c3_ind)*rscale; % cm
Lc(c3_ind) = L0(c3_ind)*rscale; % cm

rc(c4_ind) = r0(c4_ind)*rscale;
Lc(c4_ind) = L0(c4_ind)*rscale;

pI = impaction(rc,Ba,n,par,dp,Qvec); % impaction in each compartment
pS = sedimentation(n,rc,Lc,Ga,par,dp,Qvec); % sedimentation in each compartment
pD = diffusion(n,rc,Lc,par,dp,Qvec); % diffusion in each compartment

if mech==1
    pS = 0; pD = 0;
elseif mech==2
    pI = 0; pD = 0;
elseif mech==3
    pI = 0; pS = 0;
end

p = pI+pS+pD-pI.*pS-pI.*pD-pS.*pD+pI.*pS.*pD;

pu = sum(cond(p(c1_ind)));
pc = sum(cond(p(c2_ind)));
prb = sum(cond(p(c3_ind)));
palv = sum(cond(p(c4_ind)));

%==========================================================================
 
% STATE EQUATIONS----------------------------------------------------------
sys(1) = 1/par.Iu*(par.Pao-Pmus-(Ru+Rc)*Q-Ptm-Pcw); % Q
sys(2) = (Ptm-Pel)/Rra; % Va
sys(3) = Q; % V
Qra = sys(2);

% lung generations 1-2
sys(4) = u*(Q*par.Cin-Q*(Nu/Vu)*pu-Q*(Nu/Vu)*(1-pu))... %inhale
    -(1-u)*(Q*(Nc/Vc)*(1-pc)-Q*(Nu/Vu)*pu-Q*(Nu/Vu)*(1-pu)); %exhale

sys(5) = u*Q*(Nu/Vu)*pu-(1-u)*Q*(Nu/Vu)*pu;

% lung generations 3-16
sys(6) = u*(Q*(Nu/Vu)*(1-pu)-Q*(Nc/Vc)*pc-Qra*(Nc/Vc)*(1-pc))... %inhale
    -(1-u)*(Qra*(Nrb/Vrb)*(1-prb)-Qra*(Nc/Vc)*pc-Q*(Nc/Vc)*(1-pc)); %exhale

sys(7) = u*Q*(Nc/Vc)*pc-(1-u)*Qra*(Nc/Vc)*pc;

% lung generations 17-24
sys(8) = u*(Qra*(Nc/Vc)*(1-pc)-Qra*(Nrb/Vrb)*prb-Qra*(Nrb/Vrb)*(1-prb))... %inhale
    -(1-u)*(Qra*(Nalv/Valv)*(1-palv)-Qra*(Nrb/Vrb)*prb-Qra*(Nrb/Vrb)*(1-prb)); % exhale % NOTE: Am I accounting for what deposited during inhale?

sys(9) = u*Qra*(Nrb/Vrb)*prb-(1-u)*Qra*(Nrb/Vrb)*prb;

% lung generation 25 (alveoli)
sys(10) = u*(Qra*(Nrb/Vrb)*(1-prb)-Qra*(Nalv/Valv)*palv)... %inhale
    -(1-u)*(-Qra*(Nalv/Valv)*palv-Qra*(Nalv/Valv)*(1-palv)); % exhale % NOTE: Am I accounting for what deposited during inhale?

sys(11) = u*Qra*(Nalv/Valv)*palv-(1-u)*Qra*(Nalv/Valv)*palv;

% Hysteresis
% sys(1) = 1/par.Iu*(par.Pao-Pmus-(Ru+Rc)*Q-Ptm-Pcw); % Q
% sys(2) = (Ptm-Pel-Pve)/Rra; % Vra
% Qra = sys(2);
% sys(3) = (Qra-Pve/par.Rve)/par.Cve; % Pve
% % sys(3) = sys(2)-(1/Rve)*(Vve/Cve); % Vve
% sys(4) = Q; % V
% 
% sys(5) = u*(Q*par.Cin-Q*(Nc/Vc)*pc-Qra*(Nc/Vc)*(1-pc))... %inhale
%     -(1-u)*(Qra*(Nra/Vra)*(1-pra)-Qra*(Nc/Vc)*pc-Q*(Nc/Vc)*(1-pc)); %exhale
% 
% sys(6) = u*Q*(Nc/Vc)*pc-(1-u)*Qra*(Nc/Vc)*pc;
% 
% sys(7) = u*(Qra*(Nc/Vc)*(1-pc)-Qra*(Nra/Vra)*pra)... %inhale
%     -(1-u)*(-Qra*(Nra/Vra)*pra-Qra*(Nra/Vra)*(1-pra)); % exhale % NOTE: Am I accounting for what deposited during inhale?
% 
% sys(8) = u*Qra*(Nra/Vra)*pra-(1-u)*Qra*(Nra/Vra)*pra;


sys = sys';

% FUNCTION TO COMPUTE CONDITIONAL PROBABILITIES
function c = cond(v)
% calculate conditional probability for each generation
    c(1) = v(1);
    for j = 2:length(v)
        pnot = prod(1-v(1:j-1));
        c(j) = pnot*v(j);
    end
end

end