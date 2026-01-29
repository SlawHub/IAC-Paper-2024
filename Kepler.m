function theta = Kepler(e,M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                          
% Solves Kepler's Equation for E using the Newton Raphson method, 
% where E(0) = M and Tol = 1e-10. Kepler's Equation: M = E - e*sin(E). 
% Using Eccentric Anomaly, outputs True Anomaly.
%                    
% INPUTS:                                                
% - Eccentricity, (-)                              
% - Initial Mean Anomaly (radians)                           
% - Tolerance: the conditional threshold for the Newton-Raphson equation.                                              
%
% OUTPUTS:                                                              
% - True Anomaly (radians/deg)                             
%                                                        
% Author: Sho Wright (sw01745)                           
% Date created: 02/07/23                                         
% Date updated: 08/07/23                                                                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defining variables & initialising
tol = 1e-12; %Tolerance
E_0 =  M; %Making initial Eccentric Anomaly equal to initial Mean Anomaly
f = E_0 - e*sin(E_0)-M; %Kepler's equation set for initial conditions
%% Newton-Raphson Iteration loop
while abs(f)>=tol %Condition for exiting loop
    df = 1-e*cos(E_0) ; %Differentiated Kepler's equation
    E = E_0-(f/df); %Newton-Raphson equation 
    %Updating variables  
    f = E - e*sin(E)-M; %Kepler's equation with new Eccentric anomaly
    E_0=E; %New Eccentric Anomaly
end
E = E_0;
%% Calculating True Anomaly
theta = 2*atan2(sqrt(1+e)*tan(E/2), sqrt(1-e)); %True Anomaly using final Eccentric Anomaly %can put abs here 
% theta = rad2deg(theta_rad); %Converting into degrees
%Keeping TA in range of 0 - 360 degrees
% if theta < 0 
%     theta = 360 + theta;
% end
if theta < 0
    theta = 2*pi + theta;
end
%% Clearing 'ans' in Command Window
if nargout<1 
    clear theta
end 
end