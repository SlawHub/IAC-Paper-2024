function [v0,vf,psi] = Lambert(r0,rf,DT,DM,psi,psi_d,psi_u,mu)
% =========================================================================
% =========================================================================
%
% [v0,vf,psi] = Lambert(r0,rf,DT,DM,psi,psi_d,psi_u,mu)
%
% Compute the initial and final velocity along the transfer trajectory
% found by solving the Lambert's problem
%
% INPUT:
% - r0      Spacecraft Initial Position Vector (3x1) [km]
% - rf      Spacecraft Final Position Vector (3x1) [km]
% - DT      Time of Flight [s]
% - DM      Direction of Motion (either +1 or -1)
% - psi     Newton's method Initial Guess
% - psi_d   Solution range lower bound
% - psi_u   Solution range upper bound
% - mu      Gravitational Parameter [km^3/s^2]
%
% OUPUT:
% - v0      Spacecraft Initial Velocity Vector (3x1) [km/s]
% - vf      Spacecraft Final Velocity Vector (3x1) [km/s]
% - psi     psi angle [rad]
% 
% Author: N.Baresi
% Date: Feb 2014
% 
% =========================================================================
% =========================================================================

% Check
assert(psi > psi_d && psi < psi_u)

% Magnitude of Position Vectors
R0 = norm(r0);
RF = norm(rf);
cosDV = dot(r0,rf)/(R0*RF);
A = DM*sqrt(R0*RF*(1 + cosDV));

if(acos(cosDV) == 0 && A == 0)
    fprintf('Sorry! Trajectory cannot be computed!\n');
    return;
end


% Defining c2 and c3
if(psi > 1e-6)
    c2 = (1-cos(sqrt(psi)))/psi;
    c3 = (sqrt(psi) - sin(sqrt(psi)))/sqrt(psi^3);
elseif(psi < -1e-6)
    c2 = (1-cosh(sqrt(-psi)))/psi;
    c3 = (sinh(sqrt(-psi)) - sqrt(-psi))/sqrt(-psi^3);
else
    c2 = 1/2;
    c3 = 1/6;
end

Tol  = 1e-12;

% Loop on the TOF
for iter = 1:100
    
    y = R0 + RF + A*(psi*c3 - 1)/sqrt(c2);
    
    if( A > 0.0 && y < 0.0)
        
        while(y < 0)
            psi = psi + 0.1;
            y = R0 + RF + A*(psi*c3 - 1)/sqrt(c2);
        end
        
    end
    
    chi = sqrt(y/c2);

    if(psi > 1e-6)
        dc2 = 1/(2*psi)*(1 - c3*psi - 2*c2);
        dc3 = 1/(2*psi)*(c2 - 3*c3);
    elseif(psi < -1e-6)
%         fprintf('Add dc2 and dc3 for hyperbolic case\n');
        dc2 = 1/(2*psi)*(1 - c3*psi - 2*c2);    % Not sure about this
        dc3 = 1/(2*psi)*(c2 - 3*c3);            % Not sure about this
    else
        dc2 = -1/24;
        dc3 = -1/120;
    end

    % Solving for Eccentric Anomaly
    K  = (chi^3*c3 + A*sqrt(y))/sqrt(mu) - DT;
    dK = ((chi^3*(dc3 - 3*c3*dc2/(2*c2)) + A/8*(3*c3*sqrt(y)/c2 + A/chi)))/sqrt(mu);
    if(abs(K/dK)<Tol)
        break;
    else
        psi = psi - K/dK;
    end
    
    if(psi < psi_d || psi > psi_u)
        iter = 100;
        break;
    end
    
    
%     fprintf('%f\t%f\t%f\n',K,dK,psi);
    
    if(psi > 1e-6)
        c2 = (1-cos(sqrt(psi)))/psi;
        c3 = (sqrt(psi) - sin(sqrt(psi)))/sqrt(psi^3);
    elseif(psi < -1e-6)
        c2 = (1-cosh(sqrt(-psi)))/psi;
        c3 = (sinh(sqrt(-psi)) - sqrt(-psi))/sqrt(-psi^3);
    else
        c2 = 1/2;
        c3 = 1/6;
    end
    
end

if(iter < 100)
    % Compute f and g functions
    f = 1 - y/R0;
    g = A*sqrt(y/mu);

    g_dot = 1 - y/RF;

    % Compute initial and final velocities
    v0 = (rf - f*r0)/g;
    vf = (g_dot*rf - r0)/g;
else
    v0  = NaN*zeros(3,1);
    vf  = NaN*zeros(3,1);
    psi = NaN;
end

end