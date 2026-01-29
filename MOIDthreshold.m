function threshold = MOIDthreshold(RV_Earth,Ye,mubar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                          
% Calculates the maximum difference between old and new spacecraft semi-
% major axes before and after fuel burn of 150m/s.
%                    
% INPUTS:                                                
% - Position and Velocity Components of Earth during Departure period
% - Spacecraft State Space at Aphelion, Ye
% - Mass-ratio Parameter
%
% OUTPUTS:                                                              
% - MOID Threshold                           
%                                                        
% Author: Sho Wright (sw01745)                           
% Date created: 26/07/23                                         
% Date updated: 27/07/23                                                                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setting up script
close all; clc; format short;
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesFontSize', 16);

AU = 1.4959787070e08; mu = 1.32712440018e11; n = sqrt((mu)/AU^3);

r_sci = zeros(length(RV_Earth),3);
v_sci = zeros(length(RV_Earth),3);
for i = 1:length(RV_Earth)
    [r_sci(i,:),v_sci(i,:)] = Rotation(Ye(end,:),RV_Earth(i,:),mubar); %Rotates starting position at Perihelion into inertial frame
end
r = r_sci(1,:)*AU; v = v_sci(1,:)*AU*n;
vnew = v + 150e-03;

% coe1 = RV2COE(r',v',mu);
% coe2 = RV2COE(r',vnew',mu);
% MOIDlimit = abs(coe2(1,1) - coe1(1,1));

a1 = -mu/2*1/((norm(v))^2/2 - mu/norm(r));
a2 = -mu/2*1/((norm(vnew))^2/2 - mu/norm(r));
threshold = abs(a2 - a1);