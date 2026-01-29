function [c,ceq] = constraintsinputD1(X,scenario,Y2,RV_Ast,r1,v1,mubar,choice2,TU,LU)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defines and outputs nonlinear equality and inequality constraints for
% the use of fmincon for leg 1 and 2 of design 1. 
%
% INPUTS:
% - Optimisation vector X (14 x 1), including:
%   - State space at SEL2, X0 (6 x 1)
%   - Time of SEL2 escape, TOF1 (1 x 1)
%   - Departure velocity, DV1 (3 x 1)
%   - Time of transfer, TOF2 (1 x 1)
%   - Relative velocity, DV2 (3 x 1)
% OUTPUTS:                                               
% - Nonlinear equality constraints, ceq
% - Nonlinear inequality constraints, c
%                                                                                
% Author: Sho Wright (sw01745)                           
% Date created: 08/08/23                                         
% Date updated: 15/10/23                                                                                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fundamental constraint - Asteroid Reached (Case 0)
ceq(1:3) = Y2(1:3) - RV_Ast(1:3); %(0)
c = [];
r1 = [r1(1) + mubar ; r1(2:3)]; %Position vector with respect to the Sun   

if scenario == 123 
    if choice2 == 1
        ceq(4) = dot(r1,v1); %Entry into Heliocentric orbit %(choice2)
        ceq(5:7) = Y2(4:6) + X(12:14) - RV_Ast(4:6); %Relative velocity evens out %(1)
    else
        ceq(4:6) = Y2(4:6) + X(12:14) - RV_Ast(4:6); %Relative velocity evens out %(1)
    end
    c(1) = norm(X(8:10)) - 0.15*TU/LU; %Searching for departure velocity compliance %(2)
    c(2) = norm(X(12:14)) - 3.0*TU/LU; %Searching for arrival (relative) velocity compliance %(3)
elseif scenario == 12
    if choice2 == 1
        ceq(4) = dot(r1,v1); %Entry into Heliocentric orbit %(choice2)
        ceq(5:7) = Y2(4:6) + X(12:14) - RV_Ast(4:6); %Relative velocity evens out %(1)
    else
        ceq(4:6) = Y2(4:6) + X(12:14) - RV_Ast(4:6); %Relative velocity evens out %(1)
    end
    c(1) = norm(X(8:10)) - 0.15*TU/LU; %(2)
elseif scenario == 13
    if choice2 == 1
        ceq(4) = dot(r1,v1); %Entry into Heliocentric orbit %(choice2)
        ceq(5:7) = Y2(4:6) + X(12:14) - RV_Ast(4:6); %Relative velocity evens out %(1)
    else
        ceq(4:6) = Y2(4:6) + X(12:14) - RV_Ast(4:6); %Relative velocity evens out %(1)
    end
    c(1) = norm(X(12:14)) - 3.0*TU/LU; %Searching for arrival (relative) velocity compliance %(3)
elseif scenario == 23
    if choice2 == 1
        ceq(4) = dot(r1,v1); %Entry into Heliocentric orbit %(choice2)
    end
    c(1) = norm(X(8:10)) - 0.15*TU/LU; %Searching for departure velocity compliance %(2)
    c(2) = norm(X(12:14)) - 3.0*TU/LU; %Searching for arrival (relative) velocity compliance %(4)
elseif scenario == 1
    if choice2 == 1
        ceq(4) = dot(r1,v1); %Entry into Heliocentric orbit %(choice2)
        ceq(5:7) = Y2(4:6) + X(12:14) - RV_Ast(4:6); %Relative velocity evens out %(1)
    else
        ceq(4:6) = Y2(4:6) + X(12:14) - RV_Ast(4:6); %Relative velocity evens out %(1)
    end
elseif scenario == 2
    if choice2 == 1
        ceq(4) = dot(r1,v1); %Entry into Heliocentric orbit %(choice2)
    end
    c(1) = norm(X(8:10)) - 0.15*TU/LU; %(2)
elseif scenario == 3
    if choice2 == 1
        ceq(4) = dot(r1,v1); %Entry into Heliocentric orbit %(choice2)
    end
    c(1) = norm(X(12:14)) - 3.0*TU/LU; %(3)
end
end