function [vInf_dep,vInf_arr,vInf_tot,Days_Dep,TOF,X0,Xf,traj] = PorkChopBuilder(mubar,Ye,RV_Earth,RV_Ast,JD_Dep,JD_Arr,psiguess,DM,psi_u,psi_d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcuates arrival, departure, and total V_infinities for each instance of
% flight time across a given time span. 
%
% INPUTS:
% - Mass-ratio parameter, mubar, (Unitless)
% - CR3BP Outputs, Ye (km & km/s)
% - Position and Velocity vectors of:
%   - Departing location
%   - Destination
% - Julian Dates in 1 day steps for:
%   - Depature window
%   - Arrival window
% - Eccentric anomaly guess, psiguess
% - Upper and lower thresholds
%   - psi_u
%   - psi_d
%
% OUTPUTS:                                               
% - Depature velocity, vInf_dep (km/s)
% - Arrival velocity, vInf_arr (km/s)
% - Total velocity, vInf_tot (km/s)
% - List of departing days, Days_Dep (days)
% - Time of Flight, TOF (s)
%                                                                                
% Author: Sho Wright (sw01745)                           
% Date created: 19/04/23                                         
% Date updated: 21/06/23                                                                                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Defining Variables
sidereal = 86400; AU = 1.4959787070e08; mu = 1.32712440018e11;
n = sqrt(mu/AU^3); %Mean motion
LU = AU; TU = 1/n; %Normalising units for Lambert solver operation 
mu = mu*TU^2/LU^3; sidereal = sidereal/TU;

%% Applying Lambert Solver
%Initialising matrices outside of loop
TOF = zeros(length(JD_Arr),length(JD_Dep)); Days_Dep = zeros(length(JD_Arr),length(JD_Dep)); Days_Arr = zeros(length(JD_Arr),length(JD_Dep));
vInf_dep = zeros(length(JD_Arr),length(JD_Dep)); vInf_arr = zeros(length(JD_Arr),length(JD_Dep)); vInf_tot = zeros(length(JD_Arr),length(JD_Dep));
r0x = zeros(length(JD_Arr),length(JD_Dep)); r0y = zeros(length(JD_Arr),length(JD_Dep)); r0z = zeros(length(JD_Arr),length(JD_Dep));
v0x = zeros(length(JD_Arr),length(JD_Dep)); v0y = zeros(length(JD_Arr),length(JD_Dep)); v0z = zeros(length(JD_Arr),length(JD_Dep));
rfx = zeros(length(JD_Arr),length(JD_Dep)); rfy = zeros(length(JD_Arr),length(JD_Dep)); rfz = zeros(length(JD_Arr),length(JD_Dep));
vfx = zeros(length(JD_Arr),length(JD_Dep)); vfy = zeros(length(JD_Arr),length(JD_Dep)); vfz = zeros(length(JD_Arr),length(JD_Dep));

% store = zeros(length(JD_Dep),4);
for i = 1:length(JD_Dep)
    [r_Aph,v_Aph] = Rotation(Ye(end,:),RV_Earth(i,:),mubar); %Rotates starting position at Perihelion into inertial frame day by day
%     r_Aph = r_sci; %*AU v_Aph = v_sci; %*(AU*n); %R and V in SI units at Perihelion in Sun-centered inertial frame
%     r_aph = RV_Earth(i,1:3)'/LU; v_aph = RV_Earth(i,4:6)'*TU/LU;

    JDi = JD_Dep(i);
   
    for j = 1:length(JD_Arr)
        rf= RV_Ast(j,1:3)'/LU; %shifting these to plot transfer time on y-axis not arrival date
        v_Ast = RV_Ast(j,4:6)'*TU/LU;
        JDf = JD_Arr(j);
        Days_Dep(j,i) = JDi - JD_Dep(1); %Building meshes for contour plot
        Days_Arr(j,i) = JDf - JD_Arr(1);
        TOF(j,i) = sidereal*(JDf - JDi); %1 day step sizes

        if TOF(j,i) > 0 
            %rf has to change for each step in TOF to account for new position 
            [v0,vf,~] = Lambert(r_Aph,rf,TOF(j,i),DM,psiguess,psi_d,psi_u,mu);
            
        else
            v0 = NaN*ones(3,1);
            vf = NaN*ones(3,1);
        end

        %Store transfer states individually for later
        r0x(j,i) = r_Aph(1,:); r0y(j,i) = r_Aph(2,:); r0z(j,i) = r_Aph(3,:);
        v0x(j,i) = v0(1,:); v0y(j,i) = v0(2,:); v0z(j,i) = v0(3,:);
        rfx(j,i) = rf(1,:); rfy(j,i) = rf(2,:); rfz(j,i) = rf(3,:);
        vfx(j,i) = vf(1,:); vfy(j,i) = vf(2,:); vfz(j,i) = vf(3,:);

        dv_dep = norm(v0 - v_Aph); 
        dv_arr = norm(v_Ast - vf); %subtracting velocities
        dv_tot = dv_dep + dv_arr;

        % Compute porkchop plot values
        vInf_dep(j,i) = dv_dep*LU/TU;
        vInf_arr(j,i) = dv_arr*LU/TU;
        vInf_tot(j,i) = dv_tot*LU/TU; 
              
    end
end
TOF = TOF*TU;

% minv = mink(min(vInf_dep),10);
% 
% traj = zeros(6,length(minv));
% X0 = zeros(6,length(minv));
% Xf = zeros(6,length(minv));
% for i = 1:length(minv)
%     [r,c] = find(vInf_dep == minv(:,i));
%     traj(:,i) = [r; c; minv(:,i); JD_Arr(:,r); JD_Dep(:,c); TOF(r,c)];
%     X0(:,i) = [r0x(r,c); r0y(r,c); r0z(r,c); v0x(r,c); v0y(r,c); v0z(r,c)];
%     Xf(:,i) = [rfx(r,c); rfy(r,c); rfz(r,c); vfx(r,c); vfy(r,c); vfz(r,c)]; 
%     if isnan(minv(:,i)) %|| any(isnan(X0(:,i)))
%         traj(:,i) = NaN*ones(6,1);
%         X0(:,i) = NaN*ones(6,1);
%         Xf(:,i) = NaN*ones(6,1);  
%     end
% end

minv = min(vInf_dep,[],'all');
[r,c] = find(vInf_dep == minv);
traj = [r; c; minv; JD_Arr(:,r); JD_Dep(:,c); TOF(r,c)];
X0 = [r0x(r,c); r0y(r,c); r0z(r,c); v0x(r,c); v0y(r,c); v0z(r,c)];
Xf = [rfx(r,c); rfy(r,c); rfz(r,c); vfx(r,c); vfy(r,c); vfz(r,c)];

if isnan(minv) %|| any(isnan(X0))
    traj = NaN*ones(6,1);
    X0 = NaN*ones(6,1);
    Xf = NaN*ones(6,1);
end