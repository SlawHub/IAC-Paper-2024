function PorkChopPlotMain(mubar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates Pork Chop plots for departure, arrival, and total V_infinity.
%
% INPUTS:
% - Mass-ratio parameter, mubar, (Unitless)
% - Position and Velocity vectors of:
%   - Departing location
%   - Destination
% 
% OUTPUTS:                                               
% - Final Pork Chop Plots
%                                                                                
% Author: Sho Wright (sw01745)                           
% Date created: 19/04/23                                         
% Date updated: 21/06/23                                                                                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setting up script
close all; clc; format short;
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesFontSize', 16);

choice1 = input('Enter 1 S/C design 1, or 2 for S/C design 2:'); %User input configuration selection
choice2 = input('Enter 1 for v-inf departure, or 2 for all 3 Pork chop plots:');
%% Loading CR3BP Outputs
switch choice1 
    case 1
        load("CR3BPOutput1.mat",'Y','Ye','te')
    case 2
        load("CR3BPOutput2.mat",'Y','Ye','te')
end

%% Retrieving Earth State
fid = fopen('txt Files\Earth.txt');
A = textscan(fid,'%f %s %f %f %f %f %f %f', 'Delimiter', ',')';
fclose(fid);
JD_Dep = A{1}';
DATEdep = A{2}';
RV_Earth = [A{3}, A{4}, A{5} A{6}, A{7}, A{8}];

%% Retrieving Asteroid State
fid = fopen('txt Files\2023GJ2.txt');
A = textscan(fid,'%f %s %f %f %f %f %f %f', 'Delimiter', ',')';
fclose(fid);
JD_Arr = A{1}';
DATEarr = A{2}';
RV_Ast = [A{3}, A{4}, A{5} A{6}, A{7}, A{8}];

%Slicing out tick labels for each axis
DATEarr = cellfun(@(x)x(1:end-8),DATEarr,'un',0);
DATEarr = cellfun(@(x)x(6:end),DATEarr,'un',0);

DATEdep = cellfun(@(x)x(1:end-8),DATEdep,'un',0);
DATEdep = cellfun(@(x)x(6:end),DATEdep,'un',0);

Ydate = DATEarr(:,1:12:end);
Xdate = DATEdep(:,1:12:end);

Ydate(:,2:2:end) = cell(1,1);
Xdate(:,2:2:end) = cell(1,1);
%% Plotting Pork Chop Plots
sidereal = 86400; AU = 1.4959787070e08; 
mu = 1.32712440018e11; n = sqrt(mu/AU^3); 
DMvec = [-1 1]; %Variables

% load('Transfer2015XK55.mat','X0min','Xfmin','trajmin')
% OrbitPlot(Y,Ye,RV_Earth,RV_Ast,X0min,Xfmin,trajmin);
switch choice2 
    case 1
        %% Quick Plotting V-inf departure 
        tic
        i = 0;

        figure      %% NOTE: For overview (0:0.05:5) / For precision (0:0.001:0.15)
        sgtitle('Aphelion - 2023 GJ2 Departure Velocity','FontSize',12);
        hold on
        for k = 1:length(DMvec)
            DM = DMvec(1,k);
            for m = 1:0.25:16 %1:0.25:16
                i = i + 1;
                psiguess = m*(pi)^2;
                %psi_u = (m+4)*(pi)^2; psi_d = (m-4)*(pi)^2; 
                psi_u = (m+1)*(pi)^2; psi_d = (m-1)*(pi)^2;
                [vInf_dep_new,~,~,Days_Dep,TOF,X0,Xf,traj] = PorkChopBuilder(mubar,Ye,RV_Earth,RV_Ast,JD_Dep,JD_Arr,psiguess,DM,psi_u,psi_d);
                
%                 if m == 1 && DM == -1 
%                     X0min = X0;
%                     Xfmin = Xf;
%                     trajmin = traj;
%                 else
%                     for n = 1:width(traj)
%                         if traj(3,n) < trajmin(3,n) && trajmin(3,n) ~= 0
%                             X0min(:,n) = X0(:,n);
%                             Xfmin(:,n) = Xf(:,n);
%                             trajmin(:,n) = traj(:,n);
%                         end
%                     end
%                 end

                if m == 1 && DM == -1
                    X0min = X0;
                    Xfmin = Xf;
                    trajmin = traj;
                elseif traj(3,:) < trajmin(3,:) && trajmin(3,:) ~= 0 && all(isfinite(X0))
                    X0min = X0;
                    Xfmin = Xf;
                    trajmin = traj;
                end

                if DM == -1
                    vInf_dep_short = vInf_dep_new;                    
                    [~,~] = contour(Days_Dep,TOF/sidereal,vInf_dep_short,0:0.001:0.15,'LineWidth',1);
                else
                    vInf_dep_long = vInf_dep_new;
                    [~,~] = contour(Days_Dep,TOF/sidereal,vInf_dep_long,0:0.001:0.15,'LineWidth',1);
                end
            end
        end
        hold off

        ylabel('Transfer Time (Months)')
        ylim([0 730]);
        yticks(0:730/24:730)
        yticklabels({'0','','','3','','','6','','','9','','','12','','','15','','','18','','','21','','','24'})
        xlabel('Departure Date (Months)')
        xlim([0 730]);
        xticks(0:730/24:730)
        xticklabels({'01/25','','','04/25','','','07/25','','','10/25','','','01/26','','','04/26','','','07/26','','','10/26','','','01/27'})

%         ylabel('Arrival Date (Days)')   
%         ylim([TOF(1,1) TOF(end,1)]) %Limit set at highest TOF in days
%         yticks(TOF(1,1):1:TOF(end,1))
%         yticklabels(Ydate)
%         xlabel('Departure Date (Days)')
%         xlim([0 Days_Dep(1,end)]) %Limit set at furthest departure day
%         xticks(Days_Dep(1,1):1:Days_Dep(1,end))
%         xticklabels(Xdate)

        grid on
        colormap("jet")
        hcb = colorbar;
        title(hcb,'Vinf (km/s)','FontSize',9)
        ax = gca;
        ax.FontSize = 9;
        ax.FontSize = 9;
        toc
        
    case 2
        %% Subplotting all 3
        tic
        figure     
        sgtitle('Aphelion - 2023 GJ2 Pork Chop Plots','FontSize',12); %Overall title
        subplot(1,3,1)
        hold on
        for k = 1:length(DMvec)
            DM = DMvec(1,k);
            for m = 1:0.25:16
                psiguess = m*(pi)^2;
                psi_u = (m+1)*(pi)^2; psi_d = (m-1)*(pi)^2;
                [vInf_dep_new,~,~,Days_Dep,TOF,~,~] = PorkChopBuilder(mubar,Ye,RV_Earth,RV_Ast,JD_Dep,JD_Arr,psiguess,DM,psi_u,psi_d);
                if DM == -1
                    vInf_dep_short = vInf_dep_new;
                    contour(Days_Dep,TOF/sidereal,vInf_dep_short,0:0.15,'LineWidth',1);
                else
                    vInf_dep_long = vInf_dep_new;
                    contour(Days_Dep,TOF/sidereal,vInf_dep_long,0:0.15,'LineWidth',1);
                end
            end
        end
        hold off
        ax = gca;
        ax.FontSize = 9;
        ax.FontSize = 9;
        title('Departure Velocity');
        ylabel('Transfer Time (Months)')
        ylim([0 730]);
        yticks(0:730/24:730)
        yticklabels({'0','','','3','','','6','','','9','','','12','','','15','','','18','','','21','','','24'})
        xlabel('Departure Date (Months)')
        xlim([0 730]);
        xticks(0:730/24:730)
        xticklabels({'01/25','','','04/25','','','07/25','','','10/25','','','01/26','','','04/26','','','07/26','','','10/26','','','01/27'})
        grid on
        colormap("jet")
        hcb = colorbar;
        title(hcb,'Vinf (km/s)','FontSize',9)
        
        subplot(1,3,2)
        hold on
        for k = 1:length(DMvec)
            DM = DMvec(1,k);
            for m = 1:0.25:16
                psiguess = m*(pi)^2;
                psi_u = (m+1)*(pi)^2; psi_d = (m-1)*(pi)^2;
                [~,vInf_arr_new,~,Days_Dep,TOF,~,~] = PorkChopBuilder(mubar,Ye,RV_Earth,RV_Ast,JD_Dep,JD_Arr,psiguess,DM,psi_u,psi_d);
                if DM == -1
                    vInf_arr_short = vInf_arr_new;
                    contour(Days_Dep,TOF/sidereal,vInf_arr_short,0:0.15,'LineWidth',1);           
                else
                    vInf_arr_long = vInf_arr_new;
                    contour(Days_Dep,TOF/sidereal,vInf_arr_long,0:0.15,'LineWidth',1);
                end
            end 
        end
        hold off
        ax = gca;
        ax.FontSize = 9;
        ax.FontSize = 9;
        title('Arrival Velocity');
        ylabel('Transfer Time (Months)')
        ylim([0 730]);
        yticks(0:730/24:730)
        yticklabels({'0','','','3','','','6','','','9','','','12','','','15','','','18','','','21','','','24'})
        xlabel('Departure Date (Months)')
        xlim([0 730]);
        xticks(0:730/24:730)
        xticklabels({'01/25','','','04/25','','','07/25','','','10/25','','','01/26','','','04/26','','','07/26','','','10/26','','','01/27'})
        grid on
        colormap("jet")
        hcb = colorbar;
        title(hcb,'Vinf (km/s)','FontSize',9)
        
        subplot(1,3,3)
        hold on
        for k = 1:length(DMvec)
            DM = DMvec(1,k);
            for m = 1:0.25:16
                psiguess = m*(pi)^2;
                psi_u = (m+1)*(pi)^2; psi_d = (m-1)*(pi)^2;
                [~,~,vInf_tot_new,Days_Dep,TOF,~,~] = PorkChopBuilder(mubar,Ye,RV_Earth,RV_Ast,JD_Dep,JD_Arr,psiguess,DM,psi_u,psi_d);
                if DM == -1
                    vInf_arr_short = vInf_tot_new;
                    contour(Days_Dep,TOF/sidereal,vInf_arr_short,0:0.15,'LineWidth',1);
                else
                    vInf_arr_long = vInf_tot_new;
                    contour(Days_Dep,TOF/sidereal,vInf_arr_long,0:0.15,'LineWidth',1);
                    
                end
            end
        end
        hold off
        ax = gca;
        ax.FontSize = 9;
        ax.FontSize = 9;
        title('Total Velocity');
        ylabel('Transfer Time (Months)')
        ylim([0 730]);
        yticks(0:730/24:730)
        yticklabels({'0','','','3','','','6','','','9','','','12','','','15','','','18','','','21','','','24'})
        xlabel('Departure Date (Months)')
        xlim([0 730]);
        xticks(0:730/24:730)
        xticklabels({'01/25','','','04/25','','','07/25','','','10/25','','','01/26','','','04/26','','','07/26','','','10/26','','','01/27'})
        grid on
        colormap("jet")
        hcb = colorbar;
        title(hcb,'Vinf (km/s)','FontSize',9) 
        toc
end

%% Calling Orbit Plotter
BestDep = DATEdep(:,trajmin(2,:));
BestArr = DATEarr(:,trajmin(1,:));
days2ast = trajmin(6,:)/sidereal;
finalTOF = days2ast + (te(end)*n^-1)/sidereal;
fprintf(['For minimal fuel transfer - \n' ...
         'Departure Date: %s \n' ...
         'Arrival Date: %s \n' ...
         'Transfer Time: %f days (~ %f months) \n' ...
         'Total Mission Duration: %f days (~ %f months) \n'], ...
         string(BestDep),string(BestArr),days2ast,days2ast/30,finalTOF,finalTOF/30)

save('Transfer2023GJ2.mat','trajmin','X0min','Xfmin')
OrbitPlot(Y,Ye,RV_Earth,RV_Ast,X0min,Xfmin,trajmin);

%% Extra Bits and Bobs
% % Returning Orbital Elements for Initial States
% fprintf('Spacecraft at Perhelion entry \n')
% coe_peri = RV2COE(r_sci, v_sci,mu); %#ok<NASGU> 
% 
% fprintf('Earth COE \n')
% coe_Earth = RV2COE(RV_Earth(1,1:3)',RV_Earth(1,4:6)',mu); %#ok<NASGU> 
% 
% fprintf('Asteroid COE \n')
% coe_Ast = RV2COE(RV_Ast(1,1:3)',RV_Ast(1,4:6)',mu); %#ok<NASGU> 
% 
% fprintf('Spacecraft COE \n')
% coe_sc = RV2COE(r0,v0,mu); %#ok<NASGU>

% %% Trying the search and replace method again
% tic
% vInf_depi = zeros(length(JD_Arr),length(JD_Dep));
% for k = 1:length(DMvec)
%     DM = DMvec(1,k);
%     for m = 1:0.25:16
%         psiguess = m*(pi)^2;
%         psi_u = (m+1)*(pi)^2; psi_d = (m-1)*(pi)^2;
%         [vInf_dep,~,~,Days_Dep,TOF,~,~] = PorkChopBuilder(mubar,Ye,RV_Earth,RV_Ast,JD_Dep,JD_Arr,psiguess,DM,psi_u,psi_d);
% 
%         if m ~= 1
%             for i = 1:length(JD_Arr)
%                 for j = 1:length(JD_Dep)
%                     if isfinite(vInf_dep(i,j)) == 1 && isnan(vInf_depi(i,j)) == 1
%                         vInf_depi(i,j) = vInf_dep(i,j);
%                     end
%                 end
%             end
%         else
%             Days_Depi = Days_Dep;
%             vInf_depi = vInf_dep;
%         end
%     end
% end
% figure      %% NOTE: Tried 0:0.03:0.15, got no data (Original: 0:0.5:5)
% sgtitle('Aphelion - 2006 SE6 Departure Velocity','FontSize',12); %Overall title
% contour(Days_Depi,TOF/sidereal,vInf_depi,0:0.05:5,'LineWidth',1);
% range = [min(vInf_depi)];
% clabel(range,'Color','k')
% ax = gca;
% ax.FontSize = 9;
% ax.FontSize = 9;
% % title('Departure Velocity');
% ylabel('Transfer Time (Months)')
% ylim([0 730]);
% yticks(0:730/24:730)
% yticklabels({'0','','','3','','','6','','','9','','','12','','','15','','','18','','','21','','','24'})
% xlabel('Departure Date (Months)')
% xlim([0 730]);
% xticks(0:730/24:730)
% xticklabels({'01/25','','','04/25','','','07/25','','','10/25','','','01/26','','','04/26','','','07/26','','','10/26','','','01/27'})
% grid on
% colormap("jet")
% hcb = colorbar;
% title(hcb,'Vinf (km/s)','FontSize',9)
% toc