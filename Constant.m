function C = Constant(mubar,Y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the effective potential and hence the Jacobi Constant as a
% sanity check for the implemented EOM. 
%
% INPUTS:
% - Mass-ratio parameter, mubar, (Unitless)			                                                     
% - Integrated EOM vector, Y
% 
% OUTPUTS:                                               
% - Jacobi Integral, C
%                                                                                
% Author: Sho Wright (sw01745)                           
% Date created: 21/03/23                                         
% Date updated: 20/04/23                                                                                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Separating components
x = Y(:,1); y = Y(:,2); z = Y(:,3);
u = Y(:,4); v = Y(:,5); w = Y(:,6); 

%% Calculating Jacobi Constant
r1 = sqrt((x +mubar).^2+y.^2+z.^2); %Relative position vector from 1st primary
r2 = sqrt((x+mubar-1).^2+y.^2+z.^2); %Relative position vector from 2nd primary
ohm = (x.^2 + y.^2)/2 + (1-mubar)./r1 + mubar./r2; %Effective potential term
C = ohm - 0.5*(u.^2 + v.^2 + w.^2); %Jacobi Constant
% [minC,maxC] = bounds(C,'all');
% disp(maxC - minC) %Unhighlight to compute total delta-C 
%% Clearing 'ans' in Command Window
if nargout<1 
    clear C
end


