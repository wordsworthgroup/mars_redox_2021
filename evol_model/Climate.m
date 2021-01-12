%% Climate.m
%
% defines properties of the climate
%
% checked by RDW 7/1/21

classdef Climate
    
    properties (SetAccess = public)
        
        F0    = 1366;    % solar flux today at 1 AU [W/m2]
        a_CO2 = 1.5e-3;  % parameter for CO2 pressure evolution [1/My]
        t_CO2 = 0.8e3;   % parameter for CO2 pressure evolution [My]
        
        A;               % parameter for CO2 pressure evolution [bar]

        Fsol_a;          % solar flux array [W/m2]
        ps_a;            % surface pressure array [Pa]
        fH2_a;           % H2 molar concentration array [mol/mol]
        Ts_a;            % surface temperature array [K]
        
        muavg            % anon. function for mean atmospheric molar mass [kg/mol]
        ps               % anon. function for total surface pressure [Pa]
        
        time_dep = 1;    % type of time-dependence for solar flux, pCO2?
        
    end
    
    methods
        
        function self = Climate(params,mars,time,pCO2_3p5_Gya)
            
            % initializes a Climate object
            
            % inputs
            %   params: general parameter structure
            %   mars: mars parameter structure
            %   pCO2_3p5_Gya: pCO2 3.5 Gya [bar]

            % load surface temperature data from PCM_LBL
            load PCM_data
            self.Fsol_a = Fsol.data;
            self.ps_a   = ps.data;
            self.fH2_a  = fH2.data;
            self.Ts_a   = Ts.data;
            
            % CO2 partial pressure multiplier
            self.A = pCO2_3p5_Gya/(1-tanh(self.a_CO2*(time.Gy - self.t_CO2)));
            
            % define anon fns. to calculate mu_avg and ps
            
            self.muavg  = @(fi) fi*params.muCO2 + (1-fi)*params.muH2;
            % input is CO2 molar concentration [mol/mol]
            % output is mean atmospheric molar mass [kg/mol]
            
            self.ps     = @(u,fi) (u*mars.g)*(self.muavg(fi)/params.muCO2)/fi;
            % input is CO2 mass path [kg/m2], CO2 molar concentration [mol/mol]
            % output is total surface pressure [Pa]
            % uses u = (mu_i/mu_avg)*(p_i/g), p_i = fi*ps (e.g. PPC p. 97)
            % u = (mu_i/mu_avg)*(fi*ps/g)
            % ps = u*g*(mu_avg/mu_i)/fi
            
        end
        
        function uCO2 = get_uCO2(self,t,mars)
            
            % calculates CO2 mass column in atm. as a function of time
            
            % inputs
            %   t: time from formation [My]
            %   mars: Mars parameter structure
            
            % output
            %   uCO2: CO2 mass column [kg/m2]
            
            % note we output mass column not pressure to avoid partial
            % pressure conversion confusions.
                        
            if(self.time_dep)
                % c.f. eqn. (20) in Methods
                % small extra term added to avoid numerical problems caused
                % by zero present-day pressure
                uCO2 = mars.bar*self.A* ...
                    (1-tanh(self.a_CO2*(t - self.t_CO2)))/mars.g ...
                    + 0.00301*mars.bar/mars.g;
            else
                % just constant 1.25 bar of CO2
                uCO2 = 1.25*mars.bar/mars.g;
            end
            
        end
        
        function Fsol = get_Fsol(self,t,mars)

            % calculates solar flux at Mars as a function of time

            % inputs
            %   t: time from formation [My]
            %   mars: Mars parameter structure 

            % output
            %   Fsol: solar flux [W/m2]

            % solar flux at Mars today [W/m2]
            Fmars0 = self.F0/mars.d^2;
            
            % time now [My]
            t0 = 4.5e3; 
            
            % mass of Sun today [kg]
            MSun = 1.9890e+30;
            
            % solar flux at Mars at time t [W/m2]
            if(self.time_dep==1)
                % uses Eq. (1) from Gough (1981)
                Fsol   = Fmars0/(1 + (2/5)*(1 - t/t0));
            elseif(self.time_dep==2)
                % uses Eq. (5), (7) and (8) from Minton & Malhotra (2007)
                Mdot = -1.26e-11*(1 + (2/5)*(1 - t/t0)).^(-0.852); % Eq. (8)
                M    = 0.974*MSun*(1 + (2/5)*(1 - t/t0)).^0.15.*(t<=2.39e3) ...
                + (MSun + Mdot.*(t-t0)).*(t>2.39e3); % Eq. (9)
                Fsol = (Fmars0./(1 + (2/5)*(1-t/t0))).*(M/MSun).^6.75;  % Eq. (5)
            elseif(self.time_dep==3)
                % constant solar flux assumption
                Fsol   = Fmars0;
            end
            
        end
        
        function Tsurf = get_Tsurf(self,Fsol,uCO2,fH2)
            
            % calculates surface temperature 
            
            % inputs
            %   Fsol: solar flux [W/m2]
            %   uCO2: CO2 mass column [kg/m2]
            %   fH2:  H2 molar concentration [mol/mol]
            %   Ts:   surface temperature array [K]
            
            % output
            %   Tsurf: surface temperature [K]
            
            % don't allow fH2 values outside the limits of array.
            % This could slightly bias results to lower temperatures in the 
            % very strongest warming episodes. Which is fine.
            fH2_lim = fH2;
            if(fH2_lim>self.fH2_a(end))
                fH2_lim=self.fH2_a(end);
            end
            
            % CO2 molar concentration [mol/mol]
            % (assumes only H2 can be non-minor additional gas)
            fCO2 = 1 - fH2_lim;
            
            % total surface pressure [Pa]
            ps_tmp = self.ps(uCO2,fCO2);
            
            % surface temperature [K]
            % interpolate database derived from PCM_LBL runs
            Tsurf  = interpn(self.Fsol_a,log10(self.ps_a),self.fH2_a,self.Ts_a,Fsol,log10(ps_tmp),fH2_lim,'spline');
            
        end
        
    end
    
end
