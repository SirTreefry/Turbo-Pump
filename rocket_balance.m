%Andrew Trefry
%Rocket engine balance tool
import combustiontoolbox.database.NasaDatabase.*
import combustiontoolbox.core.*



go = 9.81;% gravity constant


F = 10; %kn
AcAt = 45;
Pc = 100; %chamber pressure
rc = 5.82;

%OX
%=================================
%input Fuel pump inputs (bar)
pi_o = 2.5;
effp_o = .76;
rho_o = 1141;

%Turbine ox supersonic inputs
Ptr_o = 13.6;
Teff_o = .27;

%Fuel
%=================================
%input pump inputs (bar)
pi_f = 3;
effp_f = .73;
rho_f = 70.85;


%Turbine Fuel supersonic inputs
Ptr_f = 17;
Teff_f = .59;


%Gas Generator
Tg = 871;
pg = 85;
rg = .9;


Cp = .19;
Ctp = 0; %is 0.5 if a booster pump is added
Cnozz = .1;
eps = 1;
mtp = Cp * Ctp;
mvalv = .02 * (F*Pc)^(.75);
minj = 0.25 * F^(0.85) ;
mcc = 0.75 * F^(.85);
mne = eps*F*(0.00225 * Cnozz + (0.225 - 0.075*Cnozz))/Pc;
meng = 1.34*(mtp + mvalv + minj + mcc + mne); %total mass flow balance for engine



%numerical section

%intial geuss
mt = .0604; %drive gass mdot

%engine mass flow
Cf = 1.8399;
Cstar = 2323.8;
meng = (F*1000)/(Cf*Cstar);
%total mass flow
mf = mt/(1+rg) + meng/(1+rc);
mo = (mt - mt/(1+rg)) + (meng - meng/(1+rc));
%fuel line pressure losses
pcool = 1*(0.15)*Pc; %cooling chamber for fuel lines
pinj_f = .2*Pc;
plines_f = -1*(0.05-0.1)*Pc;
pvalves_f = -1*(0.05-0.1)*Pc;

%ox line pressure losses
pinj_o = .2*Pc;
plines_o = -1*(0.05-0.1)*Pc;
pvalves_o = -1*(0.05-0.1)*Pc;


pd_o = Pc + pinj_o + plines_o + pvalves_o; %discharge pressure bar
pd_f = Pc + pcool + pinj_f + plines_f + pvalves_f; %discharge pressure


head_o = (pd_o*100000 - pi_o*100000)/(rho_o*go);
head_f = (pd_f*100000 - pi_f*100000)/(rho_f*go);


%turbine power (loop section now) cp values taken from CEA
cp = 8.0701*1000;
gamma = 1.3680;

rpo = pg/Ptr_o;
Pt_o = Teff_o*(mt)*cp*Tg*(1 - (1/rpo)^((gamma - 1)/gamma));


rpf = pg/Ptr_f;
Pt_f = Teff_f*(mt)*cp*Tg*(1 - (1/rpf)^((gamma - 1)/gamma));

%pump power flow
Pp_o = (go*head_o*mo)/effp_o;

Pp_f = (go*head_f*mf)/effp_f;


%================================
%for loop matlab section
Pt = (Pp_o +  Pp_f) - (Pt_f + Pt_o);




