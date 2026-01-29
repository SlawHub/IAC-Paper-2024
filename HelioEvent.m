function [value,isterminal,direction] = HelioEvent(~,Y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is called upon to stop integration once the s/c enters 
% heliocentric apohelion and/or perihelion.
%
% INPUTS:
% - Normalised position and velocity components in vector X at given time
%  step t (Unitless)				                                                     
% - Time step of integration, t (s)
% 
% OUTPUTS:                                               
% - Value hits zero at entry into heliocentric orbit via perihelion
% - isterminal terminates at entry into heliocentric orbit
% - Direction determines when to locate 'Value'
%                                                                                
% Author: Sho Wright (sw01745)                           
% Date created: 01/04/23                                         
% Date updated: 19/04/23                                                                                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mubar = 3.0032080443e-06; %Mass-ratio parameter
x = Y(1,:); y = Y(2,:); z = Y(3,:); %Slicing out position components
r1 = [(x + mubar) ; y ; z]; %Position vector with respect to the Sun (Earth)
vel = Y(4:6,:); %Intertial velocity vector
event = dot(r1,vel); %Calculating point at heliocentric orbit

value = [event ;event]; %ode45 stops when these values hit zero
isterminal = [1 ; 1];  %Tagging local minimum, terminating at local maximum
direction  = [-1; 1]; %Locating zero at decreasing event function, and 
% locating zero at increasing event function. 

% isterminal = [0 ; 0]; 
% direction  = []; %Use this to tag a full homoclinic orbit dt = 95.6
end

% isterminal = [0 ; 1];  %Tagging local minimum, terminating at local maximum
% direction  = [-1; 1]; %Locating zero at decreasing event function, and 
% locating zero at increasing event function. 
