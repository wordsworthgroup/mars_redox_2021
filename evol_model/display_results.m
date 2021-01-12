%% display_results.m
%
% display results of redox evolution model
%
% checked by RDW 1/11/21

close all
clearvars -except run_type

% save the data 
if(run_type==1)
    load mars_redox_data_Fig1
elseif(run_type==2)
    load mars_redox_data_FigS2
elseif(run_type==3)
    load mars_redox_data_Fig2
end

if(run_type==1)
    for ip=1:length(pCO2_3p5Gya)
        display_results_A(time,params,mars,N_a(:,ip),Tsurf_a(:,ip),uCO2_a(:,ip),cv(ip),ip,ip==2,ip==2);
    end
elseif(run_type==2)
    for ik=1:length(beta)
        display_results_B(time,params,mars,N_a(:,ik),Tsurf_a(:,ik),ik)
    end
elseif(run_type==3)
    num_tot = 1; num_now = 1; tag = '';
    display_results_C(time,params,mars,avgN,sigN,beta,Notot,num_tot,num_now,tag)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = display_results_A(time,params,mars,N_a,Tsurf_a,uCO2_a,cv,ip,plot_Fig1A,plot_Fig1B)

% displays results for Fig. 1 in main text

% inputs
%   time: time structure
%   params: general parameter structure
%   mars: mars parameter structure
%   N_a: number of moles oxidizing power in atmosphere [mol]
%   Tsurf_a: surface temperature [K]
%   uCO2_a: mass column of CO2 in atmosphere [kg/m2]
%   plotting variables

N_a_oxid  = N_a.*(N_a>0);                               % moles of oxidizing power [mol]
N_a_redu  = -N_a.*(N_a<0);                              % moles of reducing power [mol]
p_O2_a    = (N_a_oxid/4)*params.muCO2*mars.g/mars.Area; % O2 partial pressure (assuming CO2 dominant) [Pa]
p_H2_a    = (N_a_redu/2)*params.muCO2*mars.g/mars.Area; % H2 partial pressure (assuming CO2 dominant) [Pa]
time.t_Gya   = (time.T - time.t_a)/1e3;                 % time from present [Gya]

% p_O2_a and p_H2_a are the partial pressures if CO2 is the dominant gas in
% the atmosphere. To see this, write u = (mu_i/mu_CO2)*(p_i/g) (e.g. PPC p.
% 97), and also note that u = (N/A)*mu_i, where u is mass column [kg/m2].
% Then (N/A)*mu_i = (mu_i/mu_CO2)*(p_i/g), or p_i = N*mu_CO2*g/A, where N
% is the # of moles of O2 or H2 (taking stoichiometry into account).

figure(1)

if(plot_Fig1A)
    % plot partial pressures of O2 and H2 in bars (Fig. 1B)
    h1 = subplot(3,1,1);
    semilogy(time.t_Gya,p_O2_a/1e5,'r',time.t_Gya,p_H2_a/1e5,'b'); hold on
    plot([4.1 4.1],[1e-4 1e3],'Color',[1 1 1]*0.5) % pre-N/Noachian boundary
    plot([3.5 3.5],[1e-4 1e3],'Color',[1 1 1]*0.5) % Noachian/Hesperian boundary
    plot([3.0 3.0],[1e-4 1e3],'Color',[1 1 1]*0.5) % Hesperian/Amazonian boundary
    xlabel('t [Gya]')
    ylabel('partial pressure [bar]')
    yticks([1e-4 1e-2 1e0])
    legend('oxidizing (O_2)','reducing (H_2)')
    axis([0 time.T/1e3 1e-4*0.99 1e0/0.99])
    set(h1, 'Xdir', 'reverse')
    hold off
end

if(plot_Fig1B)
    % plot surface temperature (Fig. 1B)
    h2 = subplot(3,1,2);
    plot(time.t_Gya,Tsurf_a,cv); hold on
    plot([4.1 4.1],[200 300],'Color',[1 1 1]*0.5)
    plot([3.5 3.5],[200 300],'Color',[1 1 1]*0.5)
    plot([3.0 3.0],[200 300],'Color',[1 1 1]*0.5)
    xlabel('t [Gya]')
    ylabel('T_{surf} [K]')
    axis([0 time.T/1e3 200 300])
    text(1,250+20,['mean T_{surf} = ' num2str(mean(Tsurf_a),3) ' K'])
    set(h2, 'Xdir', 'reverse')
end

% plot cumulative warm period (Fig. 1C)
h3 = subplot(3,1,3);
tmax = 5;
plot(time.t_Gya,cumtrapz(time.t_a,Tsurf_a>params.Tmelt),cv); hold on
text(1,ip,['p_{CO2} 3.5 Gya = ' num2str(uCO2_a(10009)*mars.g/mars.bar,2) ' bar'],'Color',cv)
if(ip==1)
    plot([4.1 4.1],[0 tmax],'Color',[1 1 1]*0.5)
    plot([3.5 3.5],[0 tmax],'Color',[1 1 1]*0.5)
    plot([3.0 3.0],[0 tmax],'Color',[1 1 1]*0.5)
end
xlabel('t [Gya]')
ylabel('t_{warm} [My]')
set(h3, 'Xdir', 'reverse')
axis([0 time.T/1e3 0 tmax])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = display_results_B(time,params,mars,N_a,Tsurf_a,tog)

% displays results for Figure S2

% inputs
%   time: time structure
%   params: general parameter structure
%   mars: mars parameter structure
%   N_a: number of moles oxidizing power in atmosphere [mol]
%   Tsurf_a: surface temperature [K]
%   tog: subplot location []

N_a_oxid  = N_a.*(N_a>0);                          % moles of oxidizing power [mol]
N_a_redu  = -N_a.*(N_a<0);                         % moles of reducing power [mol]
p_O2_a    = (N_a_oxid/4)*params.muCO2*mars.g/mars.Area; % O2 partial pressure (assuming CO2 dominant) [Pa]
p_H2_a    = (N_a_redu/2)*params.muCO2*mars.g/mars.Area; % H2 partial pressure (assuming CO2 dominant) [Pa]
time.t_Gya   = (time.T - time.t_a)/1e3;                     % time [Gya]

figure(1)

% plot partial pressures of O2 and H2
h1 = subplot(5,2,(tog-1)*2+1);
semilogy(time.t_Gya,p_O2_a/1e5,'r',time.t_Gya,p_H2_a/1e5,'b'); hold on
xlabel('t [Gya]')
ylabel('partial pressure [bar]')
legend('oxidizing (O_2)','reducing (H_2)')
axis([0 time.T/1e3 1e-5 1e0])
set(h1, 'Xdir', 'reverse')
yticks([1e-4 1e-2 1e0])

% plot surface temperature
h2 = subplot(5,2,tog*2);
plot(time.t_Gya,Tsurf_a,'k'); hold on
xlabel('t [Gya]')
ylabel('T_{surf} [K]')
axis([0 time.T/1e3 200 300])
text(1,250,['mean T_{surf} = ' num2str(mean(Tsurf_a),3) ' K'])
text(1,235, ['t_{warm} = ' num2str(trapz(time.t_a,Tsurf_a>params.Tmelt)) ' My'])
set(h2, 'Xdir', 'reverse')

set(gcf, 'Position', [100 100 800 800])

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% displays results for Fig. 2 in main text and Supplementary Fig. 4

function [] = display_results_C(time,params,mars,avgN,sigN,beta,Notot,num_tot,num_now,tag)

% inputs
%   time: time structure
%   params: general parameter structure
%   mars: mars parameter structure
%   avgN: average cumulative warm duration in the Noachian [My]
%   sigN: std. dev. of cumulative warm duration in the Noachian [My]
%   beta: array of max. possible input increase in 1 timestep values
%   Notot: array of total reducing power added values [mol]
%   plotting variables

c            = 'krgbc';                   % line color string
nn           = newline;                   % shorthand notation
il1          = 1e2;                       % range for black lines in plot
il2          = 4.5e9;                     %
GEL_2_mols   = mars.N_e/(time.T*time.My); % convert GEL H2O equiv. to mol/s of H
My2y         = 1e6;                       % convert My to y

avgN = squeeze(avgN);
sigN = squeeze(sigN);

validate = 0;

% set minimum warm time to 5 years (~ brief impact warming episodes total)
avgN(avgN<5e-6) = 5e-6; 
sigN(sigN<5e-6) = 5e-6; 

% draw the bold colored lines
if(num_tot>1)
subplot(num_tot/2,2,num_now)
end
for ik=1:length(beta)
    loglog(Notot/(time.T*time.My),avgN(:,ik)*My2y,c(ik),'LineWidth',1.5); hold on
end


% draw the translucent shaded areas
for ik=1:length(beta)
    x  = Notot/(time.T*time.My);          % mean input rate of H [mol/s]
    y1 = (avgN(:,ik) - sigN(:,ik))'*My2y; % lower bound on total warm time [y] 
    y2 = (avgN(:,ik) + sigN(:,ik))'*My2y; % upper bound on total warm time [y]
    y1(y1<0) = 1e-8;
    y2(y2<0) = 1e-8;
    idx = y1>0 & y2>0;
    h   = patch([x(idx) fliplr(x(idx))], [y1(idx) fliplr(y2(idx))],c(ik),'LineStyle','none');
    set(h,'facealpha',.2)
end

xlabel('average H input [mol/s]')
ylabel('integrated warm period [y]')
N_a1_dot = 1e1; N_a2_dot = 1e7;
axis([N_a1_dot N_a2_dot 1e2 4.5e9])

geo_t_low  = 1e4; % lower time limit from geology [y]
geo_t_high = 1e7; % upper time limit from geology [y]

% set up colored patches for warm regime / geological constraints
patch([N_a1_dot N_a2_dot N_a2_dot N_a1_dot],[geo_t_low geo_t_low geo_t_high geo_t_high],[0.1 0.9 0.1],'FaceAlpha',0.1,'LineStyle','None')
patch([N_a1_dot N_a2_dot N_a2_dot N_a1_dot],[0.7e9 0.7e9 4.5e9 4.5e9],[1 0.75 0.0],'FaceAlpha',0.3,'LineStyle','None')
loglog([N_a1_dot N_a2_dot],[1e4 1e4],'Color',[0.5 0.5 0.5],'LineWidth',0.5)
loglog([N_a1_dot N_a2_dot],[1e7 1e7],'Color',[0.5 0.5 0.5],'LineWidth',0.5)
loglog([N_a1_dot N_a2_dot],[0.7e9 0.7e9],'Color',[0.5 0.5 0.5],'LineWidth',0.5)

text(1.5e1,1.7e9,'continuously warm regime','FontSize',11)
text(1.5e1,3e5,['geological' newline 'constraints'],'FontSize',11)

% add labels for various input processes

% Lower limit: 17.5 m GEL of H2O from Grott et al. 2011, QFM buffer, scaled
% to 73 ppmw H2O mantle limit from McCubbin et al. 2012.
% Upper limit: 61 m GEL of H2O from Grott et al. 2011, IW buffer, scaled
% to 290 ppmw H2O mantle limit from McCubbin et al. 2012.

% 0.025, 0.65 factors are the outgassed H2:H2O ratio
% we multiply by GEL_2_mols to get mol/s of H.
mantle_H2O_low   = 73; % ppmw H2O in Mars mantle (McCubbin et al. 2012)
mantle_H2O_high  = 290; % ppmw H2O in Mars mantle (McCubbin et al. 2012)
mantle_H2O_Grott = 100; % ppmw H2O in Mars mantle (Grott et al. 2011)
volc_lower_limit = 17.5*0.025*(mantle_H2O_low/mantle_H2O_Grott)*GEL_2_mols
volc_upper_limit = 61*0.65*(mantle_H2O_high/mantle_H2O_Grott)*GEL_2_mols

% calculation for impactors, following Haberle et al. (2019)
rhoc             = 2000; % assumed impactor density [kg/m3]
qFe              = 0.27; % mass fraction of Fe [kg/kg]
muFe             = 55.8450; % molar mass of Fe [kg/mol]
Dt               = @(Dc) (7e3)^0.15*Dc.^0.85 % transient crater diameter [m], from Haberle et al. 2019 eqn. S1
r_i              = @(Dc) (0.00409)*Dt(Dc).^(1/0.79) % impactor radius [m], from Haberle et al. 2019 eqn. S3
% checked for Hellas in their supp info spreadsheet: Dc = 2070 km, Dt = 881.7 km, ri = 274.4/2 km
impact_mol_H2    = @(Dc) (4/3)*pi*r_i(Dc).^3*rhoc*qFe/(muFe/1e3);
% factor two in equation below for Fe+H2O-->FeO+H2 reaction stoichiometry
impact_mol_Hs    = @(Dc) 2*impact_mol_H2(Dc)/(time.Gy*time.My*4.1);

if(validate==1)
    % for check vs. Haberle et al. 2019
    mol_1bar_CO2     = mars.N_CO2(mars.bar/mars.g);
    mol_CO2          = mol_1bar_CO2;
    ptot             = @(Dc) (mol_CO2*params.muCO2 + impact_mol_H2(Dc)*params.muH2)*mars.g/(mars.Area);

    impact_pi_H2     = @(Dc) ptot(Dc).*(impact_mol_H2(Dc)./(impact_mol_H2(Dc)+mol_1bar_CO2)); % partial pressure [bar]
    impact_pi_H2([381 1315 1352 2070]*1e3)/mars.bar
    
end

load haberle_crater_data.mat
Dc            = Craterdata*1e3; % crater diameters [m]
all_impactors = sum(impact_mol_Hs(Dc)) % [mol H/s]

% Antoniadi example for Methods text
r_i(381*1e3)/1e3 % [km]
impact_mol_H2(381*1e3)*2 % [mol H2]
impact_mol_H2(381*1e3)*2/(time.dt*time.My) % [mol H/s]

% beta values for Methods text
beta_Antoniadi = impact_mol_H2(381*1e3)*2/(time.dt*time.My)/all_impactors
beta_Argyre    = impact_mol_H2(1315*1e3)*2/(time.dt*time.My)/all_impactors
beta_Hellas    = impact_mol_H2(2070*1e3)*2/(time.dt*time.My)/all_impactors

crustal_alt    = 143*GEL_2_mols

dx             = 1.2

plot([1 1]*crustal_alt,[il1 il2],'k','LineWidth',1.5)
text(dx*crustal_alt,1e9,['crustal' newline 'alteration']);

plot([1 1]*all_impactors,[il1 il2],'k:','LineWidth',1.5)
text(dx*all_impactors,1e9,['bolide impacts']);

plot([1 1]*volc_lower_limit,[il1 il2],'k--','LineWidth',1.5)
text(dx*volc_lower_limit,1e9,'volcanism');

plot([1 1]*volc_upper_limit,[il1 il2],'k--','LineWidth',1.5)

% make font a bit bigger
text(3e0,5e0,tag,'FontSize',15)

text(2e6,4e6/3.5^0,'\beta = 1',   'Color','k')
text(2e6,4e6/3.5^1,'\beta = 10',  'Color','r')
text(2e6,4e6/3.5^2,'\beta = 10^2','Color','g')
text(2e6,4e6/3.5^3,'\beta = 10^3','Color','b')
text(2e6,4e6/3.5^4,'\beta = 10^4','Color','c')

end

