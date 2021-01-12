%% Solver.m
%
% defines core solver and associated saved variables
%
% checked by RDW 7/1/21

classdef Solver
    
    properties (SetAccess = public)
        
        N_a;      % number of moles oxidizing power in atmosphere [mol]
        dNdt_w_a; % rate of change in N due to oxidative weathering [mol/My]
        dNdt_e_a; % rate of change in N due to escape [mol/My]
        dNdt_s_a; % rate of change in N due to reducing gas supply [mol/My]
        uCO2_a;   % CO2 mass column [kg/m2]
        f_H2_a;   % H2 molar concentration [mol/mol]
        Tsurf_a;  % surface temperature [K]

    end
        
    methods

        function self = Solver(time)

            % initializes a Solver object

            % inputs
            %   time: Time structure []

            self.N_a      = zeros(1,time.nt);
            self.dNdt_w_a = zeros(1,time.nt);
            self.dNdt_e_a = zeros(1,time.nt);
            self.dNdt_s_a = zeros(1,time.nt);
            self.uCO2_a   = zeros(1,time.nt);
            self.f_H2_a   = zeros(1,time.nt);
            self.Tsurf_a  = zeros(1,time.nt);

        end
        
        function self = solve_system(self,time,params,mars,supply,escape,weathering,climate)
            
            % solves the redox system by Euler timestepping
            
            % inputs
            %   time: time structure
            %   params: general parameter structure
            %   mars: mars parameter structure
            %   supply: Supply object
            %   escape: Escape object
            %   weathering: Weathering object
            %   climate: Climate object
            
            % clear all the output
            self.N_a      = zeros(1,time.nt);
            self.dNdt_w_a = zeros(1,time.nt);
            self.dNdt_e_a = zeros(1,time.nt);
            if(params.fresh_i==1)
                self.dNdt_s_a = zeros(1,time.nt);
            end
            
            % N is moles of net oxidizing power in atmosphere [mol]
            N = 0.0;
            for it=1:time.nt
                
                % update the time structure
                time.it   = it;
                time.t    = it*time.dt;
                
                % get change due to gas input
                % with option to recalculate or just use previously saved values
                if(params.fresh_i==1)
                    dNdt_s = supply.get_redu_input(time);
                else
                    dNdt_s = self.dNdt_s_a(it);
                end
                
                % get change due to escape
                dNdt_e = escape.get_escape(climate,mars,time,N);
                
                % get change due to surface oxidation
                dNdt_w = weathering.get_weathering(time,N);
                
                % do a forward Euler timestep
                N      = N + (dNdt_w + dNdt_e + dNdt_s)*time.dt;
                
                % update arrays
                self.N_a(:,time.it)    = N;
                self.dNdt_w_a(time.it) = dNdt_w;
                self.dNdt_e_a(time.it) = dNdt_e;
                self.dNdt_s_a(time.it) = dNdt_s;
                
            end
            
        end
        
    end
    
end
