function [r_sci,v_sci] = Rotation(Ye,RV_Earth,mubar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using rotation matrices to transport coordinates of spacecraft in rotating
% frame into Heliocentric inertial reference frame. 
%
% INPUTS:
% - Final position and velocity components of satellite in rotating
% reference frame at entry into Heliocentric orbit at the Perihelion, Ye
% - Mass-ratio parameter, mubar, (Unitless)
% 
% OUTPUTS:                                               
% - Final position vector of spacecraft with respect to the inertial
% Heliocentric frame at the Perihelion, rf
%                                                                                
% Author: Sho Wright (sw01745)                           
% Date created: 06/04/23                                         
% Date updated: 23/04/23                                                                                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Slicing inputs & Translating
translate = zeros(3,1); translate(1,1) = mubar;
r_sc_syn = Ye(end,1:3)' + translate;
v_sc_syn = Ye(end,4:6)';
Pos_Earth = RV_Earth(:,1:3)';
Vel_Earth = RV_Earth(:,4:6)';

%% Building Rotation Matrix
h = cross(Pos_Earth,Vel_Earth); %Calculating Specific Angular Momentum (km^2/s)

e_r = Pos_Earth/norm(Pos_Earth); %Axis parallel to the velocity direction
e_h = h/norm(h); %Axis perpendicular to orbital plane
e_theta = cross(e_h,e_r); %Axis completing the orthogonal set

IO = [e_r, e_theta, e_h];

%% Rotating into inertial frame
%Defining rotation matrices 
n = [0 ; 0 ; 1]; %Mean motion at unity

%Rotating position and velocity vectors into inertial frame 
r_sci = IO*r_sc_syn;
v_sci = IO*(v_sc_syn + cross(n,r_sc_syn));

%% Clearing 'ans' in Command Window
if nargout<1 
    clear r_sci v_sci
end