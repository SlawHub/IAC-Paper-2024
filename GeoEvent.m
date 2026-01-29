function [value,isterminal,direction] = GeoEvent(~,Y)
mubar = 3.0032080443e-06; %Mass-ratio parameter  %% Remove 
x = Y(1,:); y = Y(2,:); z = Y(3,:); %Slicing out position components
r1 = [(x - (1 - mubar)) ; y ; z]; %Position vector with respect to the (Earth)
vel = Y(4:6,:); %Intertial velocity vector
event = dot(r1,vel); %Calculating point at heliocentric orbit

value = [event ;event]; %ode45 stops when these values hit zero
isterminal = [1 ; 1];  %Tagging local minimum, terminating at local maximum
direction  = [-1; 1]; %Locating zero at decreasing event function, and 
% locating zero at increasing event function. 

% isterminal = [0 ; 0]; 
% direction  = []; %Use this to tag a full homoclinic orbit dt = 139.8
end