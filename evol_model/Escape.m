%% Escape.m
%
% defines properties of the atmospheric escape
%
% checked by RDW 7/1/21

classdef Escape
    
    properties (SetAccess = public)
        
        Thomo_equi;  % homopause temperature [K]
        H_CO2;       % CO2 scale height [m]
        A;           % binary diffusion coeff. [molec./m/s/K^s]
        s;           % binary diffusion coeff. []
        b_H2CO2;     % H2-CO2 binary diffusion coeff. [molec./m/s]
        PhiH2_diff;  % diff-lim H2 loss flux function [H2 molecules/m2/s]
        PhiO_photoc; % photochemical escape O loss flux function [O atom/m2/s]
        PhiO_nther;  % nonthermal O loss flux function [O atom/m2/s]
        conv_Phi_F;  % Phi to dNdt conversion factor [molec./m2/s --> mol/My]
        td_DH;       % H from H2O escape decline timescale [My]
        alpha_DH;    % H from H2O escape present-day parameter []
        
    end
    
    methods
        
        function self = Escape(time,params,mars)
            
            % initializes an Escape object
            
            % inputs
            %   t: Time structure []
            %   p: Parameter structure []
            %   mars: Mars structure []
            
            % homopause temperature is not well constrained but likely
            % a few hundred Kelvin. Escape rate is pretty insensitive to it
            % in any case b/c of combined H_CO2 and b_H2CO2 dependence.
            
            self.conv_Phi_F = mars.Area*time.My/params.NA; % Phi to dNdt conversion factor [molec/m2/s --> mol/My]
            self.Thomo_equi = 500.0; % homopause temperature [K]
            self.H_CO2      = params.kB*self.Thomo_equi/(params.m_u*params.muCO2*1e3*mars.g); % scale height [m]
            self.A          = 22.3;  % coefficient for H2-CO2 diffusion
            self.s          = 0.75;  % coefficient for H2-CO2 diffusion
            self.b_H2CO2    = 1e2*self.A*self.Thomo_equi.^self.s/1e-16;   % binary diffusion coeff. [molec./m/s]
            % from Chamberlain & Hunten p. 439 Table VII.1
            
            self.td_DH    = 1.1e3; % values from plot_D_H_figure.m
            self.alpha_DH = 0.48;
            
            % non-thermal O escape rate [atom/m2/s]
            % Fig. 8 in Jakosky et al. 2018: 10^27 atoms/s 3.5 Ga, around
            % 10^26 today.
            %e.PhiO_nther = @(N,t) -1e26/mars.Area; % lower limit
            %e.PhiO_nther = @(N,t) -1e27/mars.Area; % upper limit
            
            % 4.3x10^25 atom/s today (Lillis et al. 2017)
            gamma = 1.7; % value from Heard & Kite (2020)
            self.PhiO_photoc = @(N,t) -4.3e25*(t/4.5e3)^(-1.2*gamma)/mars.Area;
            self.PhiO_nther  = @(N,t) self.PhiO_photoc(N,t)*(t>=1e3) + self.PhiO_photoc(N,1e3)*(t<1e3);
            
            self.PhiH2_diff = @(fH2) fH2.*self.b_H2CO2./self.H_CO2; % H2 diffusion-limited flux [molec./m2/s]
            
        end
        
        function dNdt_e = get_escape(self,climate,mars,time,N)
            
            % calculate oxidation rate through escape
            
            % inputs
            %   climate: Climate object []
            %   mars: Mars structure []
            %   time: time structure []
            %   N: atmospheric oxidizing power [mol]
            % output
            %   dNdt_e: time derivative of N due to preferential H loss to space [mol/My]
            
            uCO2 = climate.get_uCO2(time.t,mars);           % CO2 mass column [kg/m2]
            
            % hydrogen escape at diffusion limit when it is present
            % non-thermal O escape (only in oxidizing atmosphere) from MAVEN
            % oxidation via water loss constrained by D/H
            
            % increase in N due to H2 escape when N<0 [mol/My]
            if(N<0)
                
                % for increased accuracy, we take a semi-implicit step for
                % diffusion-limited escape, assuming it is the dominant process
                % (good approximation just after addition of a large amount of
                % H to the system)
                tau_diff = self.H_CO2*mars.N_CO2(uCO2)/(self.conv_Phi_F*self.b_H2CO2);
                dNdt_e   = N*(exp(-time.dt/tau_diff) - 1)/time.dt;
                
                % we don't take more hydrogen than we have
                if(dNdt_e>abs(N/time.dt))
                    dNdt_e = abs(N/time.dt);
                end
                
            else
                
                dNdt_e = 0;
                
            end
                        
            % add decrease in N due to O escape [mol/My]
            % O escape proceeds at all times (presence of CO2 etc.)
            dNdt_e = dNdt_e + 2*self.conv_Phi_F*self.PhiO_nther(N,time.t);
            
            % add increase in N due to H2O photolysis + H escape [mol/My]
            % H2O escape proceeds at all times
            dNdt_e = dNdt_e + self.alpha_DH*mars.H_loss_now*exp(-(time.t-time.T)/self.td_DH);
            
        end
        
    end
    
end

