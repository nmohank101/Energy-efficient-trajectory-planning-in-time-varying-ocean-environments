% clc
% close all
% clear all


%% Defining parameters 

% general 
envi_type = 'real';               % specify the type of environment 'real'
conv_factor = 1;                   
t_segment = 1;                    % time of each segment in seconds 
c_d = 3;                          % drag co-efficient
v_max = 1;                        % maximum velocity of the vehicle in m/s
min_dist = 1 ;                     % minimum distance treshold
 
% for sampling based methods 
no_of_paths = 20;          % number of noisy trajectories sampled around the current best trajectory 
num_its = 50;              % maximum number of iterations 
decay_fact = 0.99;         % decay factor to reduce the amount of exploration over iterations 
mag_step = .000009; %15;   % amount of exploration allowed .0007(2); .0001(1) 
threshold = 10000;         % used to test convergence 
N = 50;                    % number of waypoints including the start and end location 


% parameters for OGD 
delay = 0.1;               % percentage of delay allowed by the user   
treshold_conv = 1000;
param_test = 1;


%% Preparing the environment to visualize

[lat,lon,q_x,q_y,q_x_m,q_y_m,u,v,mag] = ocean_environment(envi_type,conv_factor) ;

% defining start and goal locations 
start_point = [239.9, 33.14, 0];
end_point = [241.15 , 32.99, 0];

% start_point = [239.6, 32.7, 0];
% end_point = [241.1 , 33.1, 0];


% selecting the boundary : adjust the free constants below (if any error) 

x_min = round((start_point(1) - lon(1))*length(lon)/(lon(end)-lon(1))) -20;
x_max = round((end_point(1) - lon(1))*length(lon)/(lon(end)-lon(1)))+  12;

y_min = round((start_point(2) - lat(1))*length(lat)/(lat(end)-lat(1))) -9; % 30  
y_max = round((end_point(2) - lat(1))*length(lat)/(lat(end)-lat(1))) + 20;       


% visualizing the ocean environment
figure 
contourf(q_x,q_y,mag,'LineColor','none');
caxis([0,max(max(mag))]); colormap (jet); 
htb = colorbar ;
htb.Label.String = 'Magnitude of ocean currents';
hold on;
idx = ~isnan(u) & ~isnan(v);
quiver(q_x(idx),q_y(idx),u(idx),v(idx),'LineWidth',1,'Color','k');
hold off
title('Ocean map in geographic space')

figure
contourf(q_x_m,q_y_m,mag,'LineColor','none');
caxis([0,max(max(mag))]); colormap (jet); 
htb = colorbar ;
htb.Label.String = 'Magnitude of ocean currents';
hold on;
idx = ~isnan(u) & ~isnan(v);
quiver(q_x_m(idx),q_y_m(idx),u(idx),v(idx),'LineWidth',1,'Color','k');
hold off
title('Ocean map in euclidean space');

figure 
contourf(q_x(y_min:y_max,x_min:x_max),q_y(y_min:y_max,x_min:x_max),...
         mag(y_min:y_max,x_min:x_max),'LineColor','none');

caxis([0,max(max(mag))]); colormap (jet); 
htb = colorbar ;
htb.Label.String = 'Magnitude of ocean currents';
hold on;
idx = ~isnan(u) & ~isnan(v);
quiver(q_x(1:4:length(lat),1:4:length(lon)),q_y(1:4:length(lat),1:4:length(lon)),...
                                              u(1:4:length(lat),1:4:length(lon)),...
                                v(1:4:length(lat),1:4:211),'LineWidth',1,'Color','k');


%% generating the initial trajectory
                   
[STLN_path,STLN_path_energy_] = generate_stln_path(lat,lon,q_x_m,q_y_m,N,start_point,end_point,...
                                                 v_max,u,v,c_d,delay);

% finding the energy incurred travel in straight line path 
disp('--------------------------------------------')
disp('--------------------------------------------')
temp1 = ['Total dist btwn start and goal locations = ',num2str(lldistkm([lat(1),lon(1)],[lat(end),lon(end)])), ' km'];
temp2 = ['Energy incurred in straight line path = ', num2str(STLN_path_energy_/1000), ' Kilo Joules'];
temp3 = ['Total time elapsed to travel in straight line path = ', num2str(sum(STLN_path(:,3))), ' sec'];
disp(temp1)
disp(temp2)
disp(temp3)


% finding the start and end points in euclidean coordinates

[~ , index_lon] = min(abs(lon - start_point(1)));
[~ , index_lat] = min(abs(lat - start_point(2)));
start_point_m = [q_x_m(index_lat,index_lon), q_y_m(index_lat,index_lon),0];

[~ , index_lon] = min(abs(lon - end_point(1)));
[~ , index_lat] = min(abs(lat - end_point(2)));
end_point_m = [q_x_m(index_lat,index_lon), q_y_m(index_lat,index_lon),0];

% generating the straight line 
[~, N_stline, STLN_path_energy] = generate_straight_line(start_point_m,end_point_m,...
                                                     u,v,q_x_m,q_y_m,v_max,t_segment,c_d,min_dist);

% Note that: 'N_stline' - denotes the total number of waypoints including start and end locations
N_seg = N_stline-1;                      % total number of slots/segments
T_tot = round(N_seg*t_segment);          % total time in seconds elapsed to travel in st.line 

T_tot_del = T_tot*(1+delay);             % total time after incorporating delay 
N_seg_del = round(T_tot_del/t_segment);  % total number of time slots with delay 
                                              
% finding the energy incurred travel in straight line path 
disp('--------------------------------------------')
disp('--------------------------------------------')
temp1 = ['Tot dist btwn start and goal = ',num2str(lldistkm([lat(1),lon(1)],[lat(end),lon(end)])), ' km'];
temp2 = ['Energy incurred in straight line path = ', num2str(STLN_path_energy/1000), ' Kilo Joules'];
temp3 = ['Total time elapsed to travel in straight line path = ', num2str(T_tot), ' sec'];
disp(temp1)
disp(temp2)
disp(temp3)





%% Running various algorithms
% STOMP EESTO and OGD

% total number of waypoints that excludes goal location
N_stomp = N;
N_eesto = N;
N_ogd = N_seg_del;


% considering the uncertainty in the ocean environment
noise_var = [0 0.01 0.02 0.04 0.05 0.07 0.09 0.1];

u_true = u;
v_true = v;


STOMP_path_min = [];
EESTO_path_min = [];

total_STOMP_energy_min = [];
total_EESTO_energy_min = [];

for j = 1:length(noise_var)
j

u_noisy = u_true + (noise_var(j))*randn(size(u_true));
v_noisy = v_true + (noise_var(j))*randn(size(v_true));

    
    
STOMP_path_all = [];
EESTO_path_all = [];

total_STOMP_energy = [];
total_EESTO_energy = [];

STOMP_vrel_all = [];
STOMP_vabs_all = [];

EESTO_vrel_all = [];
EESTO_vabs_all = [];

num_its_loop = 10;
for i = 1:num_its_loop
    
i     

% STOMP    
[STOMP_path, STOMP_energy,STOMP_cost,STOMP_vrel,STOMP_vabs] = STOMP_dylan_mohan(v_max,no_of_paths,...
                                          num_its,decay_fact,...
                                          N_stomp,threshold,mag_step,q_x_m,q_y_m,u_noisy,v_noisy,mag,STLN_path,c_d,...
                                          x_min,x_max,y_min,y_max,u_true,v_true);
STOMP_path_all(:,:,i) = STOMP_path; 
total_STOMP_energy = [total_STOMP_energy STOMP_energy]; 
STOMP_vrel_all(:,:,i) = STOMP_vrel;
STOMP_vabs_all(:,:,i) = STOMP_vabs;

% EESTO                                      
[EESTO_path, EESTO_energy,EESTO_cost,EESTO_vrel,EESTO_vabs] = EESTO_dylan_mohan(v_max,no_of_paths,num_its,decay_fact,...
                                          N_eesto,threshold,mag_step,q_x_m,q_y_m,u_noisy,v_noisy,mag,STLN_path,c_d,...
                                          x_min,x_max,y_min,y_max,u_true,v_true);

EESTO_path_all(:,:,i) = EESTO_path;                                      
total_EESTO_energy = [total_EESTO_energy  EESTO_energy];
EESTO_vrel_all(:,:,i) = EESTO_vrel;
EESTO_vabs_all(:,:,i) = EESTO_vabs;
                                      
end

[~,ind_stomp] = min(total_STOMP_energy);
[~,ind_eesto] = min(total_EESTO_energy);

STOMP_path_min(:,:,j) = STOMP_path_all(:,:,ind_stomp); 
EESTO_path_min(:,:,j) = EESTO_path_all(:,:,ind_eesto); 

total_STOMP_energy_min = [total_STOMP_energy_min; total_STOMP_energy]; 
total_EESTO_energy_min = [total_EESTO_energy_min; total_EESTO_energy]; 

end



% running OGD algorithm for different noise statistics 
total_OGD_energy = [];
total_OGD_path = [];

for i = 1:length(noise_var)
    
    % note that i==1 represents no noise case
    
    % adding zero mean noise to the ocean current measurements 
    u_noisy = u_true + (noise_var(i))*randn(size(u_true));
    v_noisy = v_true + (noise_var(i))*randn(size(v_true));
    
   
    [OGD_path, OGD_energy,OGD_vrel,OGD_vabs] = OGD_dylan_mohan(v_max,N_ogd,q_x_m,q_y_m,u_noisy,v_noisy,delay,lat,lon,c_d,...
                         start_point,end_point,treshold_conv,x_min,x_max,y_min,y_max,u_true,v_true);                                                       

    
    
    total_OGD_energy = [total_OGD_energy OGD_energy];
%     total_OGD_path(:,:,i) = OGD_path(:,1:2);
    
end

% when the goal is mobile 
[OGD_path, OGD_energy,OGD_vrel,OGD_vabs] = OGD_dylan_mohan_moving_goal(v_max,200,q_x_m,q_y_m,u,v,mag,delay,lat,lon,c_d,...
                               start_point,end_point,treshold_conv,x_min,x_max,y_min,y_max);                                                       

% [OGD_path, OGD_energy,OGD_vrel,OGD_vabs] = OGD_dylan_mohan_temp(v_max,N_ogd,q_x_m,q_y_m,u,v,mag,delay,lat,lon,c_d,...
%                                start_point,end_point,treshold_conv,x_min,x_max,y_min,y_max);                                                       

%% visualization 

% plotting the trajectories
figure
hold on
contourf(q_x_m(y_min:y_max,x_min:x_max),q_y_m(y_min:y_max,x_min:x_max),mag(y_min:y_max,x_min:x_max),'LineColor','none');
caxis([0,max(max(mag))]); colormap (parula); 
htb = colorbar ;
htb.Label.String = 'Magnitude of ocean currents';
quiver(q_x_m(y_min:4:y_max,x_min:4:x_max),q_y_m(y_min:4:y_max,x_min:4:x_max),u(y_min:4:y_max,x_min:4:x_max),v(y_min:4:y_max,x_min:4:x_max),'LineWidth',1,'Color','k');
% plotting the minimum energy trajectory
h1 = plot(STLN_path(:,1),STLN_path(:,2),'--w','LineWidth',3);
[~,ind] = min(total_STOMP_energy_min(5,:));
STOMP_path_min_energy = STOMP_path_all(:,:,ind);
h2 = plot(STOMP_path_min_energy(:,1),STOMP_path_min_energy(:,2),'k','LineWidth',3);
[~,ind] = min(total_EESTO_energy);
EESTO_path_min_energy = EESTO_path_all(:,:,ind);
h3 = plot(EESTO_path_min_energy(:,1),EESTO_path_min_energy(:,2),'g','LineWidth',3);
h4 = plot(OGD_path(:,1),OGD_path(:,2),'m','LineWidth',3);
legend([h1,h2,h3,h4],{'Straight line','STOMP','EESTO','OGD'})
hold off
 
% boxplots
figure
h=boxplot([total_STOMP_energy'/1000,total_EESTO_energy'/1000,repmat(OGD_energy/1000,size(total_STOMP_energy'))],'Labels',{'STOMP','EESTO','OTOGD'},'Widths',.2,'Positions',[0.1,0.4,0.7]);
set(h,'LineWidth',2);
title('Energy Cost');


% plotting the trajectories
ind = 5;
figure
hold on
contourf(q_x_m(y_min:y_max,x_min:x_max),q_y_m(y_min:y_max,x_min:x_max),mag(y_min:y_max,x_min:x_max),'LineColor','none');
caxis([0,max(max(mag))]); colormap (parula); 
htb = colorbar ;
htb.Label.String = 'Magnitude of ocean currents';
quiver(q_x_m(y_min:1:y_max,x_min:1:x_max),q_y_m(y_min:1:y_max,x_min:1:x_max),u(y_min:1:y_max,x_min:1:x_max),v(y_min:1:y_max,x_min:1:x_max),'LineWidth',1,'Color','k');
% plotting the minimum energy trajectory
h1 = plot(STLN_path(:,1),STLN_path(:,2),'--w','LineWidth',3);

STOMP_path_min_energy = STOMP_path_min(:,:,ind);
h2 = plot(STOMP_path_min_energy(:,1),STOMP_path_min_energy(:,2),'k','LineWidth',3);
EESTO_path_min_energy = EESTO_path_min(:,:,ind);
h3 = plot(EESTO_path_min_energy(:,1),EESTO_path_min_energy(:,2),'g','LineWidth',3);
h4 = plot(OGD_path(:,1),OGD_path(:,2),'m','LineWidth',3);
legend([h1,h2,h3,h4],{'Straight line','STOMP','EESTO','OGA'})
h5 = plot(OGD_path(:,1),OGD_path(:,2),'r','LineWidth',3);
legend([h1,h2,h3,h4,h5],{'Straight line','STOMP','EESTO','IOGA','OGA'})
hold off



% preparing the data
data = [total_STOMP_energy_min' total_EESTO_energy_min' ... 
                     repmat(total_OGD_energy,size(total_EESTO_energy_min',1),1)]; 


figure;
var_per = repmat({'                    0%' '                    1%' '                    2%' '                    4%' '                    5%' '                    7%' '                    9%' '                    10%'}, 1,3);
simobs = [repmat({'S '},1,8),repmat({'E '},1,8),repmat({'I '},1,8)];
boxplot(data,{var_per,simobs},'colors',repmat('krm',1,length(noise_var)),'factorgap',[5,2],'labelverbosity','minor');
h11 = findobj(gca,'Tag','Median');
set(findobj(get(h11(1), 'parent'), 'type', 'text'), 'fontsize', 20,'FontWeight','bold','VerticalAlignment','top');
X = [h11.XData]; Y = [h11.YData];
X_avg = reshape(X,2,24)';
X_avg = sum(X_avg,2)/2;
Y = Y(1:2:48);
hold on;
h11 = plot(X_avg(1:3:end),Y(1:3:end),'mo-','LineWidth',2);
h22 = plot(X_avg(2:3:end),Y(2:3:end),'ro-','LineWidth',2);
h33 = plot(X_avg(3:3:end),Y(3:3:end),'ko-','LineWidth',2);
legend([h11,h22,h33],{'IOGA','STOMP','EESTO'})


% relative velocity profiles
figure
hold on
[~,ind] = min(total_STOMP_energy);
plot(STOMP_vrel_all(:,:,ind),'b','LineWidth',3)
[~,ind] = min(total_EESTO_energy);
plot(EESTO_vrel_all(:,:,ind),'m','LineWidth',3)
plot(OGD_vrel,'g','LineWidth',3)
legend('STOMP','EESTO','OGD');
title('relative velocity profile') 


% absolute velocity profiles
figure
hold on
[~,ind] = min(total_STOMP_energy);
plot(STOMP_vabs_all(:,:,ind),'b','LineWidth',3)
[~,ind] = min(total_EESTO_energy);
plot(EESTO_vabs_all(:,:,ind),'m','LineWidth',3)
plot(OGD_vabs,'g','LineWidth',3)
legend('STOMP','EESTO','OGD');
title('absolute velocity profile') 




%% figure
% plot(STOMP_cost,'r');
% start_point = [239.65, 34.12, 0];
% end_point = [240.3, 33.8, 0];
% start and end point in meters
% 
% 
% 
% 
% plot([239.45 239.45],[32.5,33.32],'k','LineWidth',2)
% 
% plot([239.45 241.34],[33.32,33.32],'k','LineWidth',2)
% 
% plot([241.34 241.34],[32.5,33.32],'k','LineWidth',4)
% 
% plot([239.45 241.34],[32.5,32.5],'k','LineWidth',2)
% 
% 
% plot([239.65 239.65],[33.05,33.32],'k','LineWidth',2)
% plot([239.65 240.07],[33.05,33.05],'k','LineWidth',2)
% plot([240.07 240.07],[33.05,33.32],'k','LineWidth',2)
% 
% 
% 
% plot([238.8 239.43],[33.92,33.92],'k','LineWidth',2)
%  
% plot([239.42 239.42],[33.92,34.34],'k','LineWidth',2)
% 
% plot([238.8 239.43],[34.34,34.34],'k','LineWidth',2)
% 
% plot([238.8 238.8],[33.92,34.34],'k','LineWidth',2)
% 
% 
% 
% contourf(q_x(y_min:y_max,x_min:x_max),q_y(y_min:y_max,x_min:x_max),...
%          mag(y_min:y_max,x_min:x_max),'LineColor','none');




























