%% plot_TSR_figure.m
%
% Calculate thermochemical sulfate reduction rates
%
% checked by RDW 14/2/20

close all
clear all

diagnostic = 0

% define initial variables
nt    = 2000;         % temporal resolution []
nz    = 1440;         % vertical resolution []
nd    = 10;           % number of hot layer depths studied []
nc    = 5;            % number of temperature cases studied []
na    = 3;            % number of atmosphere cases studied []
ca    = 1;            % atmosphere case chosen (Stef-Boltz FYS) []
t     = zeros(1,nt);  % time [s]
d     = zeros(1,nz);  % depth [m]
T     = zeros(nz,nt); % temperature [K]
dl_a  = zeros(1,nd);  % hot layer depth [m]
r_d_a = zeros(nc,nd); % TSR alteration depth [m]

% load layer thickness data [m]
fname  = '../../scripts/mark_data/layer_heat/trials.csv';
s      = importdata(fname);
trial  = s.data(:,1);
dlayer = s.data(:,2);
Atype  = s.data(:,3);
Tlayer = s.data(:,4);

% load 1D thermal diffusion temperature data and compute TSR alteration thickness
for ii=1:na*nc*nd
    
    ii
    
    % load Mark's data: time, depth and temperature
    fname = ['../../scripts/mark_data/layer_heat/' num2str(ii-1) '.csv'];
    sd     = importdata(fname);
    d     = sd.data(:,1);     % depth [m]
    T     = sd.data(:,2:end); % temperature [K]
    for it=1:nt
        t1_data = sd.textdata{it+1};
        t2_data = t1_data(29:end);
        t(it)   = str2double(t2_data(1:end-8)); % time [s]
    end
    dt = mean(diff(t)); % time interval [s] 
    dz = mean(diff(d)); % gridscale thickness [m]
    t  = (1:nt)*dt;     % time [s]

    % sort into the various cases    
    j  = floor((ii-1)/(nc*na))+1; % hot layer thickness index []
    if(Atype(ii)==ca) % only select cases with chosen surface BC

        % hot layer thickness [m]
        dl_a(j) = dlayer(ii); 
        
        % layer temperature index []
        ic      = (Tlayer(ii) - 600)/200; 
        
        % solve the TSR equation to get total fraction of S reduced [mol/mol]
        fSred   = solve_chemistry(t,T); 
        
        % remove hot layer portion
        fSred(d<=dlayer(ii)) = 0.0; 
        
        % integrate to get TSR depth (eqn. (22) in Methods)
        r_d_a(ic,j) = trapz(d,fSred);

        % show some diagnostic information
        if diagnostic && Tlayer(ii)==1000
        
            subplot(1,2,1)
            plot(T(:,1),d,'k:',T(:,end),d,'k--',max(T,[],2),d,'k'); axis ij; 
            xlabel('T [K]'); ylabel('d [m]'); legend('start','end','max')

            subplot(1,2,2)
            plot(fSred,d); axis ij;
            xlabel('f_{S,red}'); ylabel('d [m]')
            title(['d_{S,red} = ' num2str(r_d_a(ic,j))])

            pause

        end        
        
    end
    
end

% plot TSR depth results for all temperatures
loglog(dl_a(:),r_d_a(1,:),'b','Linewidth',2); hold on
plot(dl_a(:),r_d_a(2,:),'c','Linewidth',2);
plot(dl_a(:),r_d_a(3,:),'g','Linewidth',2);
plot(dl_a(:),r_d_a(4,:),'r','Linewidth',2);
plot(dl_a(:),r_d_a(5,:),'m','Linewidth',2);
xlabel('hot layer depth [m]','FontSize',12)
ylabel('sulfate layer removal depth [m]','FontSize',12)
% add impactor lines: data from Toon et al. 2010 Table 1 column 5
depth_Toon = [1.67e0 3.37e1 2.12e2 3.47e2];
plot(([1 1]'*depth_Toon)',[1e-1 1e3],'k--')
axis([1e0 5e2 1e-1 1e3]); 
% add impactor labels
text(depth_Toon(1)/1.35,300,'Schiaparelli','Rotation', 90,'FontSize',12)
text(depth_Toon(2)/1.35,300,'Argyre',      'Rotation', 90,'FontSize',12)
text(depth_Toon(3)/1.35,300,'Isidis',      'Rotation', 90,'FontSize',12)
text(depth_Toon(4)/1.35,300,'Hellas',      'Rotation', 90,'FontSize',12)
legend('800 K','1000 K','1200 K','1400 K','1600 K')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fSred = solve_chemistry(t,T_a)

% data from Truche et al. 2009: see Table 2 expt. Ti13 and main text
Rstar = 8.314462;    % ideal gas constant [J/K/mol]
k0    = 10.^(-7.15); % reaction rate at reference temperature (280 C) [1/s]
Ea    = 131e3;       % activation energy [J/mol]
Tref  = 273.15+280;  % reference temperature [K]

% anonymous function for k(T) [1/s]
k     = @(T) k0*exp(-(Ea/Rstar)*(1./T - 1/Tref)); 

% calculate fraction of sulfate reduced [mol/mol]
% c.f. eqn. (25) in Methods
fSred = 1 - exp(-trapz(t,k(T_a(:,:)),2));

return

end

