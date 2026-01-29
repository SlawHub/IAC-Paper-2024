function [v1_vecs, v2_vecs] = lambert(r1_vec, r2_vec, tof, mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [v1_vecs, v2_vecs] = lambert_izzo(r1_vec, r2_vec, tof, mu)
%
% Solve the Lambert's problem via Izzo's method
%
% INPUTS:
% .r1_vec   departure position vector (3x1)
% .r2_vec   arrival positi on vector (3x1)
% .tof      time of flight
% .mu       gravitational parameter of attracting mass
%
% OUTPUTS:
% .v1_vecs  departure velocity vectors (3xN)
% .v2_vecs  arrival velocity vectors (3xN)
%
% Notes: N is the number of Lambert solver's solutions
% 
% Author: N. Baresi
% History:
% .05/07/2024 - File created
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(length(r1_vec)==3,'Stop here traveller! r1_vec must be three-dimensional');
assert(length(r2_vec)==3,'Stop here traveller! r2_vec must be three-dimensional');
assert(tof > 0,'Stop here traveller! tof must be positive');
assert(mu > 0,'Stop here traveller! mu must be positive');

[v1_vecs, v2_vecs] = lambert_solver(r1_vec, r2_vec, tof, mu);

v1_vecs = reshape(v1_vecs, 3, []);
v2_vecs = reshape(v2_vecs, 3, []);

v1_vecs(:, vecnorm(v1_vecs)==0) = NaN;
v2_vecs(:, vecnorm(v2_vecs)==0) = NaN;
