function [dXdt] = EOM2BP(t,X,mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 4). 
% i). Rewrite EOM as a system of 1st order ODEs. Write a function that 
% outputs the RHS of the newly found system of 1st order ODEs. Using the 
% time vector from W3.1 and the s/c ECI-frame initial conditions in W3.3,
% integrate the EOM of the ECI 2BP using ode45. 
%                                                         
% INPUTS:
% - Time, t (s)
% - Position and velocity components in vector X, (km,km/s)				
% - Earth's gravitational parameter, mu (km^3/s^2)
%                                                       
% OUTPUTS:                                               
% - System of first order ODEs, dxdt
%                                                                                
% Author: Sho Wright (sw01745)                           
% Date created: 09/11/22                                         
% Date updated: 25/10/22                                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defining position and velocity coordinates 
posX = X(1,1); posY = X(2,1); posZ = X(3,1);
velx = X(4,1); vely = X(5,1); velz = X(6,1);
%% Assigning Inputs
r = [posX;posY;posZ]; %Position vector (km)
rnorm = norm(r); %Norm of position vectors(km)
rdot = [velx;vely;velz]; %Velocity vector (km/s)
%% Initialise output 
dXdt = zeros(6,1); %Creating empty 6x1 vector
% Differential equations
dXdt(1:3,:) = rdot; %Velocity 
dXdt(4:6,:) = -(mu/rnorm^3)*r; %Acceleration
end