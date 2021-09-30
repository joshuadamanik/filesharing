%/////////////////////////////////////////////////////////////////////////%
%                                                                         %
%   - Name: Environment.m                                                 %
%                                                                         %
%                           - Created by C. H. Lee, 2020. 02. 10.         %
%                                                                         %
%/////////////////////////////////////////////////////////////////////////%

function out = Environment( Alt )

%.. Parameters 

    Re            	=   6370987.308 ; 
    E_mass       	=   5.973e24 ; 
    P_sl         	=   101325 ; 
    R               =  	287.053 ;   
    G           	=   6.673e-11 ;
	gamma       	=   1.4 ;    
    
%.. Compute Environment Variables

    if  Alt < 11000                                                         % Troposphere Case 
         
        T           = 	288.15 - 0.0065 * Alt ;                             % Temperature             	[K]
        Pressure    =   P_sl * ( T / 288.15 )^5.2559 ;                      % Atmospheric Pressure  	[Pa]
        
    else                                                                 	% Stratosphere Case
        
        T           = 	216 ;                                               % Temperature               [K]
        Pressure    = 	22630 * exp( - 0.00015769 * ( Alt - 11000 ) ) ;     % Atmospheric Pressure   	[Pa]
        
    end
    
    Rho             = 	Pressure / ( R * T ) ;                              % Atmospheric Density   	[kg/m^3]
    Grav            = 	G * E_mass/( Re + Alt )^2 ;                         % Gravitational Acc.        [m/s^2]           
	V_sound         = 	sqrt( gamma * R * T ) ;   
    
%.. Exporting Data

    out     =   [ Rho ; Grav ; V_sound ] ; 
    
end