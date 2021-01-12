%% Weathering.m
%
% defines properties of the surface oxidative weathering
%
% checked by RDW 7/1/21

classdef Weathering
    
    properties (SetAccess = public)
        
        weather_rate; % weathering rate [mol/My]
        
    end
    
    methods
        
        function self = Weathering(weather_rate)

            % initializes a Weathering object

            % inputs
            %   weather_rate: weathering rate [mol/My]
            
            self.weather_rate = weather_rate;
            
        end
        
        function dNdt_w = get_weathering(self,time,N)
            
            % calculates loss of atm. oxidizing power due to weathering
            
            % inputs
            %   time: time structure []
            %   N: atmospheric oxidizing power [mol]
            % output
            %   dNdt_w: time derivative of N due to oxidation of basalt [mol/My]
            
            if(N>0)
                
                dNdt_w = (N>0)*(-self.weather_rate);
                
                % we don't take more oxygen than we have
                if(abs(dNdt_w)>N/time.dt)
                    dNdt_w = -(N/time.dt);
                end
                
            else
                dNdt_w = 0;
            end
            
        end
        
    end
    
end









