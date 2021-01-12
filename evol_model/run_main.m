%% run_main.m
%
% run main script for redox evolution box model
%
% checked by RDW 7/1/21

tic

close all
clear all

run_type = 1; % chose type of run to perform
    
% check for saved data data 
if(run_type==1)
    fname = 'mars_redox_data_Fig1.mat';
elseif(run_type==2)
    fname = 'mars_redox_data_FigS2.mat';
elseif(run_type==3)
    fname = 'mars_redox_data_ensemble.mat';
    %parpool('local',str2num(getenv('SLURM_CPUS_PER_TASK')))
end

if(exist(fname))
    error('Saved data already exists for this run type, please move it before continuing.')
end

% set up everything that doesn't change
model_setup
g_exp = -1.95; % power law distribution exponent []
runtime = time.My*time.T; % total model runtime [s]
weather_rate = 1.58e3*time.My; % weathering rate [mol/My]

if run_type==1 % varying pCO2 parametrization
    
    nR          = 1;             % number of realizations for each parameter set
    Notot       = 1.5e4*runtime; % total reducing power added [mol]
    beta        = 1e3;           % max. possible input increase in 1 timestep
    pCO2_3p5Gya = [2 1.5 1];
    cv          = 'mkc';
    
elseif run_type==2 % varying beta
    
    nR          = 1;
    Notot       = 1.5e4*runtime;
    beta        = [1 10 100 1000 10000];
    pCO2_3p5Gya = 1.5;
    
elseif run_type==3 % big ensemble runs
    
    nR          = 256;
    Notot       = logspace(2,6,128)*params.exa;
    beta        = logspace(0,4,5);
    pCO2_3p5Gya = 1.5;
    
end

nP = length(pCO2_3p5Gya);
nN = length(Notot);
nB = length(beta);

% create output arrays
warm_timeT = zeros(nP,nN,nB,nR); % total time where Ts > 273 K [My]
warm_timeN = zeros(nP,nN,nB,nR); % Noachian time where Ts > 273 K [My]
inp_tot    = zeros(nP,nN,nB,nR); % total reducing gas input [mol]
esc_tot    = zeros(nP,nN,nB,nR); % total reducing gas escape [mol]
wea_tot    = zeros(nP,nN,nB,nR); % total oxidative weathering [mol]
Ts_avg     = zeros(nP,nN,nB,nR); % mean surface temperature [K]

reducing_frac_noach  = zeros(nP,nN,nB,nR);
oxidizing_frac_noach = zeros(nP,nN,nB,nR);
reducing_frac_rest   = zeros(nP,nN,nB,nR);
oxidizing_frac_rest  = zeros(nP,nN,nB,nR);

N_a     = zeros(time.nt,nP);
Tsurf_a = zeros(time.nt,nP);
uCO2_a  = zeros(time.nt,nP);

if(run_type==1 || run_type==2)
    % initialize a Solver object
    solver = Solver(time);
end

% loop over pCO2_3p5Gya
for ip=1:length(pCO2_3p5Gya)
    
    ip
    
    % initialize a Climate object
    climate = Climate(params,mars,time,pCO2_3p5Gya(ip));
    
    % loop over reducing gas total input
    for im=1:nN
        
        im
        
        % loop over peak input rates
        for ik=1:nB
            
            ik
            
            % create necessary objects
            supply = Supply(time,g_exp,beta(ik),Notot(im));
            escape = Escape(time,params,mars);
            weathering = Weathering(weather_rate);
            
            % loop over realization
            %parfor in = 1:nR % use this for parallel calcs on Cannon
            %%%% comment this out for parallel runs %%%%
            for in = 1:nR
            
                if(run_type==3)
                    % initialize a Solver object
                    solver = Solver(time);
                end
            
                % solve system
                solver = solver.solve_system(time,params,mars,supply,escape,weathering,climate);
                
                % get fH2 and Tsurf vs. time
                for it=1:time.nt
                    tt                 = it*time.dt; % time elapsed [My]
                    Fsol               = climate.get_Fsol(tt,mars); % solar flux [W/m2]
                    solver.uCO2_a(it)  = climate.get_uCO2(tt,mars); % CO2 mass column [kg/m2]
                    solver.f_H2_a(it)  = (solver.N_a(it)<0).*abs(solver.N_a(it)/2)/(abs(solver.N_a(it)/2) + mars.N_CO2(solver.uCO2_a(it))); % H2 molar concentration [mol/mol]
                    solver.Tsurf_a(it) = climate.get_Tsurf(Fsol,solver.uCO2_a(it),solver.f_H2_a(it)); % surface temperature [K]
                end

                % freeze supply vector if performing variable CO2 runs
                if(run_type==1 && ip==1)
                    %%%% comment this out for parallel runs %%%%
                    params.fresh_i = 0;
                end
                
                % calculate diagnostics
                inp_tot(ip,im,ik,in)    = sum(solver.dNdt_s_a)*time.dt;
                esc_tot(ip,im,ik,in)    = sum(solver.dNdt_e_a)*time.dt;
                wea_tot(ip,im,ik,in)    = sum(solver.dNdt_w_a)*time.dt;
                warm_timeT(ip,im,ik,in) = sum((solver.Tsurf_a>params.Tmelt))*time.dt;
                warm_timeN(ip,im,ik,in) = sum((solver.Tsurf_a>params.Tmelt)&(time.t_a>time.tNoach))*time.dt;
                
                % calculate for 4.1 to 3.5 Gya, and 3.5 Gya to present day
                % if resolution changes, so must ncut!
                ncutA = 4000;
                ncutB = 10000;
                reducing_frac_noach(ip,im,ik,in)  = sum(solver.N_a(ncutA:ncutB)<0)/(ncutB-ncutA)
                oxidizing_frac_noach(ip,im,ik,in) = sum(solver.N_a(ncutA:ncutB)>0)/(ncutB-ncutA);
                reducing_frac_rest(ip,im,ik,in)   = sum(solver.N_a(ncutB:end)<0)/(time.nt-ncutB)
                oxidizing_frac_rest(ip,im,ik,in)  = sum(solver.N_a(ncutB:end)>0)/(time.nt-ncutB);
                
                % save display data specific to run types 1-2
                %%%% comment this out for parallel runs %%%%
                if(run_type==1)
                    N_a(:,ip)     = solver.N_a;
                    Tsurf_a(:,ip) = solver.Tsurf_a;
                    uCO2_a(:,ip)  = solver.uCO2_a;
                elseif(run_type==2)
                    N_a(:,ik)     = solver.N_a;
                    Tsurf_a(:,ik) = solver.Tsurf_a;
                    uCO2_a(:,ik)  = solver.uCO2_a;
                end

            end
        end
    end
end

% save the data 
if(run_type==1)
    save mars_redox_data_Fig1
elseif(run_type==2)
    save mars_redox_data_FigS2
elseif(run_type==3)
    % calculate some final diagnostics
    avgT = mean(warm_timeT,4);
    avgN = mean(warm_timeN,4);
    sigN = std(warm_timeN,0,4);    
    save mars_redox_data_ensemble
end
toc

if(run_type<3)
    display_results
end

%export_fig Fig_1.pdf
%print(gcf, '-dpdf', 'Fig_1.pdf')

