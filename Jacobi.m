function [J] = Jacobi(mubar,X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Building the Jacobian matrix by linearisation of EOM. This script is 
% required to provide stability analysis of the Lagrange points, and is 
% called upon in 'LagrangeSolver.m'.
% 
% To test this function copy and paste:
%   - Jacobi(3e-06,[-1;0.0;0;0;0;0])
%
% INPUTS:
% - Mass-ratio parameter, mubar, (Unitless)
%   - mubar = m2/m1+m2
%   - where m1 and m2 are the masses of the primaries (kg)
%   - E.g., for a Sun-Earth system:
%        Mass of Earth: m1 = 5.972E24 (kg)
%        Mass of Sun:   m2 = 1.989E30 (kg)
% - Normalised position and velocity components in vector X, (Unitless)				                                                   
% 
% OUTPUTS:                                               
% - 6x6 Jacobian matrix, J (Unitless)
%                                                                                
% Author: Sho Wright (sw01745)                           
% Date created: 08/03/23                                         
% Date updated: 08/03/23                                                                                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defining submatrices
Zeros = zeros(3,3); %3x3 matrix of zeros
I = eye(3,3); %3x3 Identity matrix

x = X(1,1); %Indexing first row and first column of vector X
r1 = norm(x + mubar); %Position vector magnitude for 1st primary
r2 = norm(x-1+mubar); %Position vector magnitude for 2nd primary
E = (1-mubar)/r1^3+mubar/r2^3; %Constant 
F =   [2*E 0  0 ; %3x3 matrix of constants for [u;v;w] components
        0 -E  0 ; 
        0  0 -E];

skew = [0 -1 0 ; %3x3 unit skew-symmetric matrix
        1  0 0 ; 
        0  0 0];
%% Building Jacobian Matrix
J = [Zeros        I   ;
     F-skew^2 -2*skew];
%% Clearing 'ans' in Command Window
if nargout<1 
    clear J
end