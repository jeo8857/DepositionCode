function [f,g] = CalculateDeposition(t,s,dp,par,res,Y,vol,breaths,Qin,Qex,Vin,Vex)

tt = mod(t,res.T);
for i = 1:length(tt)
    if tt(i)>=0 && tt(i)<=res.TI
        u(i) = 1;
    elseif tt(i)>res.TI && tt(i)<=res.T
        u(i) = 0;
    end
end

Q = s(:,1); 

% Particles
% Nc = s(:,4);
% Dc = s(:,5);
% Nra = s(:,6);
% Dra = s(:,7);
% Nc = s(:,5);
% Dc = s(:,6);
% Nra = s(:,7);
% Dra = s(:,8);


% Nu = s(:,4);
% Du = s(:,5);
% Nc = s(:,6);
% Dc = s(:,7);
% Nra = s(:,8);
% Dra = s(:,9);

Nu = s(:,4);
Du = s(:,5);
Nc = s(:,6);
Dc = s(:,7);
Nrb = s(:,8);
Drb = s(:,9);
Nalv = s(:,10);
Dalv = s(:,11);

% Ntotal = Nc+Nra;
% Dtotal = Dc+Dra;
Ntotal = Nu+Nc+Nrb+Nalv;
Dtotal = Du+Dc+Drb+Dalv;



TIind = find(abs(res.TI-t)==min(abs(res.TI-t)));
In1 = trapz(t(1:TIind),u(1:TIind)'.*Q(1:TIind).*par.Cin); % total amount breathed in
    
Tind2 = find(abs((breaths-1)*res.T-t)==min(abs((breaths-1)*res.T-t))); % index at the beginning of the last breathing cycle  

In2 = In1+Ntotal(Tind2); % The amount available to deposit during last breathing cycle (particle in + particles suspended at beginning of cycle)
TIind2 = find(abs(((breaths-1)*res.T+res.TI)-t)==min(abs(((breaths-1)*res.T+res.TI)-t))); % index at end of steady state inhalation
 
% figure(20)
% hold on
% plot(t(Tind2:end),Ntotal(Tind2:end)./In2)
% 
% figure(21)
% hold on
% plot(t(Tind2:end),Dtotal(Tind2:end)./In2)
% stop
% stop

% plot(t,Dtotal)
% hold on
% plot(t(end),Dtotal(end),'ro')
% plot(t(Tind2),Dtotal(Tind2),'ko')
% plot(t(TIind2),Dtotal(TIind2),'bo')
% stop
% Deposited at end of last breath
f(1) = (Du(end)-Du(Tind2))/In2;
f(2) = (Dc(end)-Dc(Tind2))/In2;
f(3) = (Drb(end)-Drb(Tind2))/In2;
f(4) = (Dalv(end)-Dalv(Tind2))/In2;
f(5) = (Dtotal(end)-Dtotal(Tind2))/In2;

% Fraction of inhaled particles that are deposited during inhale
f(6) = (Du(TIind2)-Du(Tind2))/In2;
f(7) = (Dc(TIind2)-Dc(Tind2))/In2;
f(8) = (Drb(TIind2)-Drb(Tind2))/In2;
f(9) = (Dalv(TIind2)-Dalv(Tind2))/In2;
f(10) = (Dtotal(TIind2)-Dtotal(Tind2))/In2;

% Fraction of inhaled particles that are deposited during exhale after
% equilibrium
f(11) = (Du(end)-Du(TIind2))/In2;
f(12) = (Dc(end)-Dc(TIind2))/In2;
f(13) = (Drb(end)-Drb(TIind2))/In2;
f(14) = (Dalv(end)-Dalv(TIind2))/In2;
f(15) = (Dtotal(end)-Dtotal(TIind2))/In2;

% Efficiency of deposition deposited during exhale. Divide the total
% deposited during exhale by the particles that were suspended at the end
% of inhale, rather than by the total number of particles inhaled.
f(16) = (Du(end)-Du(TIind2))/Ntotal(TIind2);
f(17) = (Dc(end)-Dc(TIind2))/Ntotal(TIind2);
f(18) = (Drb(end)-Drb(TIind2))/Ntotal(TIind2);
f(19) = (Dalv(end)-Dalv(TIind2))/Ntotal(TIind2);
f(20) = (Dtotal(end)-Dtotal(TIind2))/Ntotal(TIind2);

% Fraction of inhaled particles that are suspended at end of inhale
g(1) = Nu(TIind2)/In2;
g(2) = Nc(TIind2)/In2;
g(3) = Nrb(TIind2)/In2;
g(4) = Nalv(TIind2)/In2;
g(5) = Ntotal(TIind2)/In2;

% Fraction of inhaled particles that are suspended at end of breath
g(6) = Nu(end)/In2;
g(7) = Nc(end)/In2;
g(8) = Nrb(end)/In2;
g(9) = Nalv(end)/In2;
g(10) = Ntotal(end)/In2;

end