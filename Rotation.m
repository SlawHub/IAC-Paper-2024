function [r,v] = Rotation(Ye,RV_Earth,mubar,FrameChange)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using rotation matrices to transport coordinates of spacecraft current
% reference frame into another frame of choice
%
% INPUTS:
% - Final position and velocity components of satellite in rotating
% reference frame at entry into Heliocentric orbit at the Perihelion, Ye
% - Mass-ratio parameter, mubar, (Unitless)
% - Current frame-Frame of choice, FrameChange
% 
% OUTPUTS:                                               
% - Final position vector of spacecraft with respect to the inertial
% Heliocentric frame at the Perihelion, rf
%                                                                                
% Author: Sho Wright (sw01745)                           
% Date created: 06/04/23                                         
% Date updated: 14/10/23                                                                                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defining Common Variables
Pos_Earth = RV_Earth(:,1:3)';
Vel_Earth = RV_Earth(:,4:6)';
if strcmpi(FrameChange,'CR3BP-SCI') == 1
%% Slicing inputs & Translating
translate = zeros(3,1); translate(1,1) = mubar;
r_sc_syn = Ye(end,1:3)' + translate;
v_sc_syn = Ye(end,4:6)';
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
r = IO*r_sc_syn;
v = IO*(v_sc_syn + cross(n,r_sc_syn));
elseif strcmpi(FrameChange,'CR3BP-ECI') == 1
    %% Slicing inputs & Translating
    translate = zeros(3,1); translate(1,1) = -1+mubar;
    r_sc_syn = Ye(end,1:3)' + translate;
    v_sc_syn = Ye(end,4:6)';
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
    r = IO*r_sc_syn;
    v = IO*(v_sc_syn + cross(n,r_sc_syn));
elseif strcmpi(FrameChange,'ECI-CR3BP') == 1
    %% Slicing inputs & Translating
    translate = zeros(3,1); translate(1,1) = 1-mubar;
    r_eci = Ye(end,1:3)';
    v_eci = Ye(end,4:6)';
    %% Building Rotation Matrix
    h = cross(Pos_Earth,Vel_Earth); %Calculating Specific Angular Momentum (km^2/s)
    e_r = Pos_Earth/norm(Pos_Earth); %Axis parallel to the velocity direction
    e_h = h/norm(h); %Axis perpendicular to orbital plane
    e_theta = cross(e_h,e_r); %Axis completing the orthogonal set
    IO = [e_r, e_theta, e_h]; %Simply take transpose of this 
    OI = IO'; %Use this DCM to rotate back from ECI to CR3BP
    %% Rotating into inertial frame
    %Defining rotation matrices 
    n = [0 ; 0 ; 1]; %Mean motion at unity
    %Rotating position and velocity vectors into inertial frame 
    
    r = OI*r_eci;
    v = OI*v_eci - cross(n,r);
    r = r + translate;
end
%% Clearing 'ans' in Command Window
if nargout<1 
    clear r v
end
end