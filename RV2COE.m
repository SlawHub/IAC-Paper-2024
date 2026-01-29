function [coe] = RV2COE(r,v,mu,units)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 6). 
% i). Explain how to derive the direction cosine matrix [IT] from the 
% tangential frame (NTW) of the satellite T to the ECI Frame I 
%                                                         
% INPUTS:
% - ECI position components, r
% - ECI velocity components, v
%                                                       
% OUTPUTS:                                               
% - ECI Classical Orbital Elements (coe):
% - coe = [a e i RAAN ArgP theta]
%                                                                                 
% Author: Sho Wright (sw01745)                           
% Date created: 09/11/22                                         
% Date updated: 13/10/23                                                                                              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defining variables
K_hat = [0;0;1]; %Parallel to Earth's rotation axis
h = cross(r,v); %Angular Momentum
%% Semi-Major Axis
E = 0.5*norm(v)^2 - mu/norm(r);
a = -mu/(2*E);
%% Eccentricity 
e_vector = (1/mu)*cross(v,h)-r/norm(r); 
e = norm(e_vector); %Eccentricity (Unitless)
%% Inclination
i = acos(h(3,1)/norm(h)); %Inclination (rad)
i_deg = rad2deg(i); 
%% RAAN
n_hat = cross(K_hat,h)/norm(cross(K_hat,h));
RAAN = atan2(n_hat(2,1),n_hat(1,1)); 
RAAN_deg = rad2deg(RAAN);
%% Argument of Perigee
n_perp = cross(h,n_hat)/norm(cross(h,n_hat));
ArgP = atan2(dot(e_vector,n_perp),dot(e_vector,n_hat)); 
ArgP_deg = rad2deg(ArgP);
%% True Anomaly 
ie = e_vector/e;
ip = cross(h, e_vector)/norm(cross(h, e_vector));
theta = atan2(dot(r,ip),dot(r,ie));
theta_deg = rad2deg(theta);
% Keeping TA in range of 0 - 360 degrees
if theta < 0 
    theta = 2*pi + theta;
end
if theta_deg < 0 
    theta_deg = 360 + theta_deg;
end
%% Choice of Radians or Degrees
switch units
    case 'rad'
        coe = [a, e, i, RAAN, ArgP, theta];
    case 'deg'
        coe = [a, e, i_deg, RAAN_deg, ArgP_deg, theta_deg];
end
%% Outputting COE 
% fprintf(['Classical Orbital Elements from position and velocity vectors:' ...
%     '\n[a e i RAAN ArgP theta] \n' ...
%     '[%.3e %.3f %.3f %.3f %.3f %.3f] \n'],coe);
end