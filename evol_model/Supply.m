%% Supply.m
%
% defines properties of the reducing gas supply function
%
% checked by RDW 11/2/20

classdef Supply
    
    properties (SetAccess = public)
        td;        % volcanism/impactor decline timescale [My]
        g;         % power law exponent []
        beta;      % factor by which to increase x1 over mu []
        Notot;     % total reducing gas input [mol]
        mu0;       % initial expectation value of distribution [mol/My]
        x0_0       % initial minimum value of distribution [mol/My]
        x1_0       % initial maximum value of distribution [mol/My]
        S_scale    % temporal scaling factor (t.t is current time in My) []
        xs_fn;     % anon. distribution function [mol/My]
        mu_fn;     % fn. for expectation value of distribution [mol/My]
    end
    
    methods
        
        function self = Supply(time,g,beta,Notot)
            
            % initializes a Supply object

            % inputs
            %   time:      time structure []
            %   g:         power law exponent []
            %   beta:      peak input multiplier []
            %   Notot:     total reducing gas input [mol]
            
            self.td    = 1.0e3;
            self.g     = g;
            self.beta  = beta;
            self.Notot = Notot;
            gg         = (g+1)/(g+2);
            time_dep   = 1;
            
            % initial expectation value of distribution [mol/My]
            % uses eqn. (11) in Methods and Notot = muav*t0
            if(time_dep)
                self.mu0 = Notot/((1 - exp(-time.T/self.td))*self.td);
            else
                self.mu0 = Notot/time.T;
            end
            log10_mu0 = log10(self.mu0);
            
            % initial maximum value of distribution
            self.x1_0 = self.beta*self.mu0;
            
            % Now we calculate x0 via Newton's method. We only need it
            % once, because afterward we can simply scale by exp(-t/tau).
            
            % anon. fn. for expectation value of distribution vs. x0
            % see eqn. (5) in Methods
            mu_fn = @(x0) gg*(self.x1_0^(g+2)-x0^(g+2))/(self.x1_0^(g+1)-x0^(g+1));
            
            % log10 form of mu_fn 
            log_mu_fn = @(log_x0) log10(mu_fn(10^log_x0));
            
            % use Newton's Method to get x0_0 from mu0
            self.x0_0 = 10.^(fzero(@(x) log_mu_fn(x) - log10_mu0,-4));
            
            % anon. fn. to sample probability distribution
            % see eqn. (7) in Methods
            self.xs_fn    = @(y,x0,x1) ((x1^(g+1) - x0^(g+1))*y + x0^(g+1)).^(1/(g+1));
            
            if(time_dep)
                % temporal scaling factor (t.t is current time in My)
                % see eqn. (10) in Methods
                self.S_scale = @(tt) exp(-tt/self.td);
            else
                self.S_scale = @(t) 1;
            end
            
        end
        
        function dNdt_s = get_redu_input(self,time)

            % calculates change of N due to reducing gas supply
            
            % inputs
            %   time: Time structure []
            % output
            %   dNdt_s: time derivative of N due to reducing gas supply [mol/My]
            
            % scale x0 and x1 to account for decay of mu with time
            % it can be demonstrated that this is equivalent to scaling mu
            x0 = self.x0_0*self.S_scale(time.t);
            x1 = self.x1_0*self.S_scale(time.t);
            
            % sample distribution to get decrease in N due to reducing gas 
            % input [mol/My]
            dNdt_s = -self.xs_fn(rand,x0,x1);
            
        end
        
    end
    
end

