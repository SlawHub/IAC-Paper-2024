function [dXdt] = EOM(mubar,X,~)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations of Motion (EOM) for the Circular Restricted Three Body 
% Problem (CRTBP). EOM are non-dimensionalised using normalising units and 
% simplified by introducing the Jacobi integral.
% 
% To test this function copy and paste:
%   - EOM(0.1,[0.5;0.1;0;0;0;0])
%
% INPUTS:
% - Mass-ratio parameter, mubar, (Unitless)
% - Normalised position and velocity components in vector X, (Unitless)			
%                                                       
% OUTPUTS:                                               
% - System of first order ODEs, dXdt
%                                                                                
% Author: Sho Wright (sw01745)                           
% Date created: 20/02/23                                         
% Date updated: 19/04/23                                                                                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defining position and velocity coordinates 
x = X(1,1); y = X(2,1); z = X(3,1);
u = X(4,1); v = X(5,1); w = X(6,1);
%% Assigning Inputs
rdot = [u;v;w]; %Velocity vector (km/s)
r1 = sqrt((x+mubar)^2+y^2+z^2);
r2 = sqrt((x+mubar-1)^2+y^2+z^2);
%% Normalised Equations of Motion
dXdt = zeros(6,1); %Creating empty 6x1 vector
% Differential equations
dXdt(1:3,:) = rdot; %Velocity X-Y-Z components 
dXdt(4,:) = x - (1-mubar)/r1^3*(x+mubar)- mubar/r2^3*(x+(mubar-1))+ 2*v; %Acceleration X-component
dXdt(5,:) = y - (1-mubar)/r1^3*y- mubar/r2^3*y-2*u; %Acceleration Y-component
dXdt(6,:) = -(1-mubar)/r1^3*z-mubar/r2^3*z; %Acceleration Z-component
end

