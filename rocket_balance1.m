%Andrew Trefry
%Rocket engine balance tool
import combustiontoolbox.database.NasaDatabase.*
import combustiontoolbox.core.*

a = .25;
res_tol = 1*10^(-6);
conv_tol = 1*10^(-6);
go = 9.81;% gravity constant


F = 62.2; %kn
AcAt = 82.9;
Pc = 36; %chamber pressure
rc = 5.15;

%OX
%=================================
%input Ox pump inputs (bar)
pi_o = 2;
effp_o = .73;
rho_o = 1141;



%Fuel
%=================================
%input pump inputs (bar)
pi_f = 3;
effp_f = .6;
rho_f = 70.85;


%Turbine Fuel supersonic inputs
Ptr_f = 16.7;
Teff_f = .45;


%Gas Generator
Tg = 860;
pg = 70;
rg = .87;


Cp = .19;
Ctp = 0; %is 0.5 if a booster pump is added
Cnozz = .1;
eps = 1;
mtp = Cp * Ctp;
mvalv = .02 * (F*Pc)^(.75);
minj = 0.25 * F^(0.85);
mcc = 0.75 * F^(.85);
mne = eps*F*(0.00225 * Cnozz + (0.225 - 0.075*Cnozz))/Pc;
meng = 1.34*(mtp + mvalv + minj + mcc + mne); %total mass flow balance for engine



%numerical section

%intial geuss

mt(1) = 3; %drive gass mdot

i = 1;
max_iter = 1000;
while i<max_iter

   
%engine mass flow
Cf = 1.8435;
Cstar = 2314.5;
meng = (F*1000)/(Cf*Cstar);
%total mass flow
mf = mt(i)/(1+rg) + meng/(1+rc);
mo = (mt(i) - mt(i)/(1+rg)) + (meng - meng/(1+rc));

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
cp = 8.1791*1000;
gamma = 1.3692;



%Eqs setups per denominator
%======================================================================



rpf =pg/Ptr_f;
%==============================================================
Pt_f = Teff_f*(mt(i))*cp*Tg*(1 - (1/Ptr_f)^((gamma - 1)/gamma)); %derivable eq

%pump power flow
Pp_o = (go*head_o*mo)/effp_o;

Pp_f = (go*head_f*mf)/effp_f;



%Eqs setups for Deriavtives
%======================================================================



rpf = pg/Ptr_f;
%==============================================================
Pt_fd = Teff_f*cp*Tg*(1 - (1/rpf)^((gamma - 1)/gamma)); %derivable eq




mfd = 1/(1+rg) + meng/(1+rc);
mod = (1 - 1/(1+rg)) + (meng - meng/(1+rc));
%pump power flow
Pp_od = (go*head_o*mod)/effp_o;

Pp_fd = (go*head_f*mfd)/effp_f;




%Physical Newton Method
%================================
%for loop matlab section
Pt(i) = (Pp_o +  Pp_f) - (Pt_f); %derviative value
ptd(i) = (Pp_od +  Pp_fd) - (Pt_fd);
mt(i+1) = mt(i) - a*(Pt(i)/ptd(i));
conv(i+1) = abs(mt(i+1) - mt(i)); % calculate iterative convergence
res(i+1) = abs(Pt(i)); % calculate the residual
 
if abs(ptd(i)) < 1e-10
    break
end


    
i=i+1; % increment our iteration counter
 

end


fprintf("Increment %f\n",i)
fprintf("Conv %f\n",conv(i-1))
fprintf("residual %f\n",res(i-1))