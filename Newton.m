function [Xstar,i] = Newton(mubar,X0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Since the Jacobi integral is the only integral of motion in CRTBP, the 
% system cannot be solved analytically and has to be solved using iterative 
% methods instead. In this script, Newton's Method is applied.
% 
% To test this function copy and paste:
%   - Newton(0.1,[-1;0;0;0;0;0])
%
% INPUTS:
% - Mass-ratio parameter, mubar, (Unitless)
%   - mubar = m2/m1+m2
%   - where m1 and m2 are the masses of the primaries (kg)
%   - E.g., for a Sun-Earth system:
%        Mass of Sun:   m1 = 1.989E30 (kg)
%        Mass of Earth: m2 = 5.972E24 (kg)
% - Initial guess, X0, (1x6 vector)                                                     
% 
% OUTPUTS:                                               
% - Final normalised position and velocity components in vector X
%   after convergence (Unitless)
% - Final iteration count to convergence (Unitless)
%                                                                                
% Author: Sho Wright (sw01745)                           
% Date created: 02/04/23                                         
% Date updated: 02/04/23                                                                                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setting up script
% close all; clc; format short;
% set(0, 'DefaultTextInterpreter', 'latex');
% set(0, 'DefaultAxesFontSize', 16);

%% Initialising
tol = 1e-10;
i = 0; %Beginning tick for each initial guess

%Calculating first guess to set Newton's condition
J0 = Jacobi(mubar,X0); %Initital Jacobian matrix
J0inv = J0\eye(6,6); %Initial inverse Jacobian matrix
deltaX0 = -J0inv*EOM(mubar,X0);
X = X0 + deltaX0; %Initial guess via Newton's Method   

while norm(EOM(mubar,X))>=tol
    i=i+1; %Iteration count
    X = X0 - J0inv*EOM(mubar,X0); %Newton's method equation    
    %Updating variables 
    X0 = X; %Setting new X vector 
    J = Jacobi(mubar,X); %Calculating new Jacobian matrix for new X
    Jinv = J\eye(6,6); %New inverse Jacobian matrix for new X
    J0inv = Jinv; %Setting new Jinv matrix
end

Xstar = X; %Redefining output

%% Clearing 'ans' in Command Window
if nargout<1 
    clear Xstar i
end