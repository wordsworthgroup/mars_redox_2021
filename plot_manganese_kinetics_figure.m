%% plot_manganese_kinetics_figure.m
%
% Calculate manganese oxidation rates and compare with impactor data
%
% checked by RDW 13/2/20
%
% Data is from Hem (1981) Table 1 and Morgan (2005) Table 1

close all
clear all

Rstar     = 8.314462;        % ideal gas constant [J/K/mol]
time.d2s     = 86400;           % 1 Earth day [s]
time.y2sMars = 668.6*88775.244; % 1 Mars year [s]
% http://en.wikipedia.org/wiki/Timekeeping_on_Mars [s]

check_prm = 1; % check parametrization vs. existing results?

% Use Hem 1981 to get an activation energy for the reaction

% data for Reaction K17 from Hem (1981) Table 1
Ta = 37.4 + 273.15; % temperature [K]
Ka = 10^(-3.1);     % reaction coefficient [L/s/m2]

% data for Reaction K13 from Hem (1981) Table 1
Tb = 24.6 + 273.15; % temperature [K]
Kb = 10^(-4.42);    % reaction coefficient [L/s/m2]

% k = A*exp(-Ea/(R*T))
% k_a/k_b = exp(-(Ea/R)*(1/T_a - 1/T_b))
% log[k_b]-log[k_a] = (Ea/R)*(1/T_a - 1/T_b)
% Ea = R*(log[k_b]-log[k_a])/(1/T_a - 1/T_b) 

Ea = +Rstar*(log(Kb)-log(Ka))/(1/Ta-1/Tb); % activation energy of reaction [J/mol]

clear Ka Tb Kb

% Henry's Law coefficient for O2 [mol/kg/Pa]
% https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Units=SI&Mask=10#Solubility
KH0  = 1.3e-8; % mol/kg/Pa
T0   = 298.15; % K
KH   = @(T) KH0*exp(1700*(1./T - 1/T0));

% effective reaction rate at room temperature [1/M/s]
% From Morgan 2005 for a pH of 8 
% M is mol/L. 1 kg/L for H2O, so L & kg interchangeable here.
k_app0 = 1e-4;

% effective reaction rate vs temperature [L/mol/s]
k_app = @(T) k_app0*exp(-(Ea/Rstar)*(1./T-1/T0));

if(check_prm)
    
    % [L/mol/s]*[Pa]*[mol/kg/Pa]=[L/kg/s]=[1/s]

    T_test = 310; % test temperature [K]
    
    % baseline case (c.f. Morgan 2005 Table 1, last entry)
    pO2_earth  = 0.21*1e5; % pO2 [Pa]
    k_check    = k_app0*pO2_earth*KH0; % rate constant [1/s]
    tau_check  = 1/(k_check*time.d2s) % reaction time constant [days]
    
    % Henry's Law temperature variation only
    % this slows the reaction down b/c less dissolved oxygen
    k_check2   = k_app0*pO2_earth*KH(T_test); % rate constant [1/s]
    tau_check2 = 1/(k_check2*time.d2s) % reaction time constant [days]
    
    % all temperature variation
    k_check3   = k_app(T_test)*pO2_earth*KH(T_test); % rate constant [1/s]
    tau_check3 = 1/(k_check3*time.d2s) % reaction time constant [days]
    
end

% create anonymous function for tau_MnOx

tau_MnOx = @(T,pO2) 1./(k_app(T).*KH(T)*pO2); % manganese oxidation timescale [s]

% load Kathryn's 3D GCM impact modeling results

data = load('steakley_impact_results/temp_150mb_RAC_nomelt_bare.txt');
ta1  = data(:,1)*time.y2sMars; % time [s]
Ts1  = data(:,3); % surface temperature after 30 km impactor, 150 mbar atmosphere [K]

data = load('steakley_impact_results/temp_500mb_RAC_nomelt_bare.txt');
ta2  = data(:,1)*time.y2sMars; % time [s]
Ts2  = data(:,3); % surface temperature after 30 km impactor, 500 mbar atmosphere [K]

data = load('steakley_impact_results/temp_1bar_RAC_nomelt_bare.txt');
ta3  = data(:,1)*time.y2sMars; % time [s]
Ts3  = data(:,3); % surface temperature after 30 km impactor, 1 bar atmosphere [K]

% display results

% plot impactor temperature curves
semilogx(ta1/time.d2s,Ts1,'Color',[0.7 0.7 0.7]); hold on
plot(ta2/time.d2s,Ts2,'k')
plot(ta3/time.d2s,Ts3,'k','Color',[0.7 0.7 0.7])

% plot manganese reaction data
Ta = linspace(250,400,1e4); % temperature array [K]
pA = 3e2; % O2 partial pressure A [Pa]
pB = 3e3; % O2 partial pressure B [Pa]
pC = 3e4; % O2 partial pressure C [Pa]
plot(tau_MnOx(Ta,pA)/time.d2s,Ta,'Color',[1 0.7 0.7],'LineWidth',1.5)
plot(tau_MnOx(Ta,pB)/time.d2s,Ta,'Color',[1 0 0],'LineWidth',1.5)
plot(tau_MnOx(Ta,pC)/time.d2s,Ta,'Color',[1 0.7 0.7],'LineWidth',1.5)

% make labels etc.
rot_val = 352.5;
h = text(1.5,325,[num2str(pA/1e5) ' bar O_2'],'Color',[1 0.7 0.7]); set(h,'Rotation',rot_val)
h = text(1.5,340,[num2str(pB/1e5) ' bar O_2'],'Color',[1 0 0]); set(h,'Rotation',rot_val)
h = text(1.5,355,[num2str(pC/1e5) ' bar O_2'],'Color',[1 0.7 0.7]); set(h,'Rotation',rot_val)
h = text(8e2,210,'150 mbar CO_2','Color',[0.7 0.7 0.7]);
h = text(8e2,250,'500 mbar CO_2','Color','k');
h = text(9e2,270,'1 bar CO_2','Color',[0.7 0.7 0.7]);
xlabel('t [Earth days]','FontSize',12)
ylabel('T_{surf} [K]','FontSize',12)
axis tight
scatter(60.03,316.31,65,'bo','filled')
axis([1 3000 200 450])
