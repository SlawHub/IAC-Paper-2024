function [MOID, RV1, RV2] = calculate_moid(coe1, coe2, GM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [MOID, RV1, RV2] = calculate_moid(coe1, coe2, GM)
%
% Calculate Minimum Orbital Intersection Distance given two
% eccentric Keplerian Orbits.
%
% Input:
% . coe1    First set of Keplerian orbit elements
% . coe2    Second set of Keplerian orbit elements
% . GMs     Gravitational parameter of attracting mass 
%
% Output:
% . MOID    Minimum Orbital Intersection Distance
% . RV1     State vector at closest approach along first orbit
% . RV2     State vector (pos+vel) at closest approach along second orbit 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define function handles for eccentric and mean anomalies
EA = @(e, TA) 2*atan2d(sqrt(1-e)*tand(TA/2), sqrt(1+e)); %put an abs here - probably not allowed
MA = @(e, EA) EA - e*sind(EA);

% convert from COE to RV
RV1 = COE2RV(coe1, GM);
r1 = RV1(1:3);
v1 = RV1(4:6);

% Calculate angular momentum and eccentricity vectors
h1 = cross(r1, v1);
ecc1 = cross(v1, h1)/GM - r1/norm(r1);

% convert from COE to RV
RV2 = COE2RV(coe2, GM);
r2 = RV2(1:3);
v2 = RV2(4:6);

% Calculate angular momentum and eccentricity vectors of second orbit
h2 = cross(r2, v2);
ecc2 = cross(v2, h2)/GM - r2/norm(r2);


% normalize unit vectors
eh1 = h1/norm(h1);
ee1 = ecc1/norm(ecc1);
ep1 = cross(eh1, ee1);

% normalize unit vectors
eh2 = h2/norm(h2);
if coe2(2,:) ~= 0
    ee2 = ecc2/norm(ecc2);
else
    ee2 = [0;0;0];
end
ep2 = cross(eh2, ee2);

% compute mutual nodal line
K = cross(eh1, eh2)/norm(cross(eh1, eh2));

% find true anomalies along first orbit
tht1a = atan2d(dot(K, ep1), dot(K, ee1));
M1a = MA(coe1(2), EA(coe1(2), tht1a));
RV1a = COE2RV([coe1(1:5); M1a], GM);

% tht1b = tht1a + pi;
tht1b = tht1a + 180;
M1b = MA(coe1(2), EA(coe1(2), tht1b));
RV1b = COE2RV([coe1(1:5); M1b], GM);

% find true anomalies along second orbit
tht2a = atan2d(dot(K, ep2), dot(K, ee2));
M2a = MA(coe2(2), EA(coe2(2), tht2a));
RV2a = COE2RV([coe2(1:5); M2a], GM);

% tht2b = tht2a + pi;
tht2b = tht2a + 180;
M2b = MA(coe2(2), EA(coe2(2), tht2b));
RV2b = COE2RV([coe2(1:5); M2b], GM);


% calculate distances across possible pairs
D = zeros(4,1);
D(1) = norm(RV1a(1:3) - RV2a(1:3));
D(2) = norm(RV1a(1:3) - RV2b(1:3));
D(3) = norm(RV1b(1:3) - RV2a(1:3));
D(4) = norm(RV1b(1:3) - RV2b(1:3));

% Store moid and retrieve position vectors
MOID = min(D);

I = find(D == MOID);

checksc = isscalar(I);
checkch = ischar(I);
if checksc == 0 && checkch == 0
    % assert(~isnan(I) && isscalar(I), 'Invalid input');
%     warning('Error caught at index %d', iter)
    I = 0;
end

switch I
    case 1
        RV1 = RV1a;
        RV2 = RV2a;
    case 2
        RV1 = RV1a;
        RV2 = RV2b;
    case 3
        RV1 = RV1b;
        RV2 = RV2a;
    case 4
        RV1 = RV1b;
        RV2 = RV2b;
%     case 0 
%         RV1 = NaN*ones(6,1);
%         RV2 = NaN*ones(6,1);
end