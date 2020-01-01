% OGD 
function [path, OGD_energy,V_rel_i,V_abs_i] = OGD_dylan_mohan(v_max,N,q_x_m,q_y_m,u,v,dely,lat,lon,c_d,...
                                              start_point,end_point,treshold_conv,...
                                              x_min,x_max,y_min,y_max,u_true,v_true)



% finding the magnitude of ocean currents
[n,m] = size(v);
mag = zeros(n,m);
for i = 1:n
    for j = 1:m
        flag = false;

        if u(i,j) == -9999
            u(i,j) = NaN;
            flag = true;
        end

        if v(i,j) == -9999
            v(i,j) = NaN;
            flag = true;
        end

        if flag == false
            mag(i,j) = sqrt(u(i,j)^2 + v(i,j)^2);
        else
            mag(i,j) = NaN;
        end
    end
end                                          
                                          
                                          
                                          
                                          
                                          
%% Creating the initial path
 
[~ , index_lon] = min(abs(lon - start_point(1)));
[~ , index_lat] = min(abs(lat - start_point(2)));
start_point_m = [q_x_m(index_lat,index_lon), q_y_m(index_lat,index_lon),0];

[~ , index_lon] = min(abs(lon - end_point(1)));
[~ , index_lat] = min(abs(lat - end_point(2)));
end_point_m = [q_x_m(index_lat,index_lon), q_y_m(index_lat,index_lon),0];


% Initializing the path variable
path = zeros(N,3);
path(1,:) = start_point_m;
path(N,:) = end_point_m;

% Calculating the stright line path step size
% x_step = abs(path(1,1)-path(num_waypoints,1)) / (num_waypoints - 1);
% y_step = abs(path(1,2)-path(num_waypoints,2)) / (num_waypoints - 1);

x_step = (path(N,1)-path(1,1)) / (N - 1);
y_step = (path(N,2)-path(1,2)) / (N - 1);

% Intializing all the waypoints in the path
for m = 2:N-1
    path(m,1) = path(m-1,1) + x_step;
    path(m,2) = path(m-1,2) + y_step;
end

% calculate distance
t_init = sqrt(x_step^2 + y_step^2) / v_max;
path(:,3) = t_init;
path(N,3) = 0;
path_init = path;

% updating the total number of waypoints 
T = t_init*N;
T_tot = T + dely*T;
N = round(T_tot/t_init);

% alternate init


%% Running OGD 

start_loc = path(1,1:2)';    
goal_loc = path(end,1:2)';
curr_loc = start_loc;

OGD_path = zeros(3,N);    % Waypoints that forms the optimal trajectory  
OGD_path(3,:) = t_init;
OGD_path(1:2,1) = start_loc;  % Start location  
alpah = 0;
m = 2;

norm_param = max(max(mag));

temp_eta = [];
temp_ang = [];


while(norm(curr_loc-goal_loc) >= treshold_conv && m <= N) % 
    
    vel_ocean = find_ocean_vel(OGD_path(1,m-1),OGD_path(2,m-1),u,v, q_x_m(1,:),q_y_m(:,1)'); % finding vel at previous waypoint 

%     vel_ocean = vel_ocean + 0.3*rand(size(vel_ocean));
    %   various strategies for lamda 
    %     lamda_t =    (i/N_final); % control parameter
    % by default use this : (1-(((goal_loc-curr_loc)/norm(goal_loc-curr_loc))'*vel_ocean))/2;
    temp_cos =   (((goal_loc-curr_loc)/norm(goal_loc-curr_loc))'*vel_ocean);
    lam_temp =  ( 1+ (((goal_loc-curr_loc)/norm(goal_loc-curr_loc))'*vel_ocean) )/2 ; 
    lamda_t =   (1 -  (norm(vel_ocean)/norm_param)* (2*lam_temp-1))*2^-(norm(vel_ocean)/norm_param);
     
   % coding new lamda from the paper
%    temppp = goal_loc-curr_loc;
   eta_ocean = (norm(vel_ocean)/norm_param);
%   lamda_t = 1 - eta_ocean * (1 + ((temppp)'*vel_ocean)/(norm(temppp)*norm(vel_ocean)))/2;
    
%      lamda_t = m/N;
    
    %(a) 1 - (norm(vel_ocean)/norm_param)* lam_temp ;  0.33*(1+2*m/N)*
    %(b) (1 -  (norm(vel_ocean)/norm_param)* (2*lam_temp-1))*2^-(norm(vel_ocean)/norm_param) 
    
    % evaluting the gradient 
    X_grad = (lamda_t*(curr_loc-goal_loc)/norm(curr_loc-goal_loc)) + ((1-lamda_t)*(-(vel_ocean)));
    X_grad = X_grad./t_init;
     
    
    
    
    
    % evaluating the step size
    temp11 = norm(curr_loc-goal_loc)/((N-m)*t_init);
    temp12 = 0.5*(1 + ((goal_loc-curr_loc)'*vel_ocean)/(norm(goal_loc-curr_loc)*norm(vel_ocean)) );
    alpah = [alpah  exp(-dely*(1/temp11)*sqrt(temp12))];
    
    
    % exp(-1.2*dely*sqrt(temp12))     *(norm(vel_ocean)/norm_param)
    % exp(-(dely)*sqrt(temp12 *(norm(vel_ocean)/norm_param) ))
    % if (2*lam_temp - 1)>=0
%         alpah = [alpah  exp(-3*dely*(1/temp11)*sqrt(temp12))]; % exp(-dely*(1/temp11)*sqrt(temp12)) %exp(-1.2*dely*sqrt(temp12))
%     else
%     alpah(m) = max(0.5,alpah(m));
     

%     vel_max_modified =   v_max*exp(-0.5*(norm(vel_ocean)/norm_param + sqrt(dely))*sqrt(temp12)); 
    
     vel_max_modified =   v_max*exp(-0.5* (dely + eta_ocean* (sqrt((1-lamda_t )/eta_ocean)) ));    
           
    
    % v_max*alpah(m)    *(1/temp11) ;
    
%     vel_ocean =  (vel_ocean'*X_grad.*t_init)*X_grad/norm(X_grad);   
        
        %      vel_max_modified = 0.5*(norm(start_point_m(1:2)-end_point_m(1:2))/T_tot +v_max)*alpah(m) ;
    % evaluating the step size 
    
    temp7 = 4*(X_grad'*vel_ocean)^2 - 4*(norm(vel_ocean)^2-vel_max_modified^2)*(norm(X_grad)^2); 
    
    ss = ((-2*X_grad'*vel_ocean) - sqrt(temp7))/(2*(norm(vel_ocean)^2-vel_max_modified ^2));
    if ss < 0
        ss = ((-2*X_grad'*vel_ocean)+sqrt(temp7))/(2*(norm(vel_ocean)^2-vel_max_modified ^2));
    end
    
    ss = 1/ss;
    X_grad = X_grad.*t_init;
    
    % performing the update 
    curr_loc = curr_loc - ss * X_grad;
    OGD_path(1:2,m) = curr_loc; 
    
    if m == N && norm(curr_loc-goal_loc) >= treshold_conv
       x_step = (goal_loc(1)-OGD_path(1,m));
       y_step = (goal_loc(2)-OGD_path(2,m));
       [ocean_xm, ocean_ym] = vel_ocean_opt(curr_loc(1), curr_loc(2), u, v, q_x_m(1,:),q_y_m(:,1)');
       curr_loc(1) = curr_loc(1) + x_step;
       curr_loc(2) = curr_loc(2) + y_step;
%        computing the segment time 
       t_init = sqrt(x_step^2 + y_step^2) / (v_max-norm([ocean_xm, ocean_ym])) ;
       m = m+1;
       OGD_path(1:2,m) = curr_loc; 
       OGD_path(3,m) = 0; 
       OGD_path(3,m-1) = t_init;
    end
    
     m = m + 1 ;
end




figure 
hold on
contourf(q_x_m(y_min:y_max,x_min:x_max),q_y_m(y_min:y_max,x_min:x_max),mag(y_min:y_max,x_min:x_max),'LineColor','none');
caxis([0,max(max(mag))]); colormap (parula); 
htb = colorbar ;
htb.Label.String = 'Magnitude of ocean currents';
quiver(q_x_m(y_min:y_max,x_min:x_max),q_y_m(y_min:y_max,x_min:x_max),u(y_min:y_max,x_min:x_max),v(y_min:y_max,x_min:x_max),'LineWidth',1,'Color','k');
h1 = plot(path_init(:,1),path_init(:,2),'--k','LineWidth',3);
h2 = plot(OGD_path(1,1:m-1),OGD_path(2,1:m-1),'m','LineWidth',3);
ylim([5*10^4,9.5*10^4]);
title('IOGA trajectory for static goal scenario','FontSize',24,'FontWeight','bold')
xlabel('Distance (m)','FontSize',24,'FontWeight','bold');
ylabel('Distance (m)','FontSize',24,'FontWeight','bold');

legend([h1,h2],{'Straight line','OGD'})

ax = gca;
ax.FontWeight = 'bold';
ax.FontSize  = 24;



  
 
 

%% finding the energy incurred in the STOMP trajectory

path = OGD_path(:,1:m-1)';
[ocean_x, ocean_y] = vel_ocean_opt(path(:,1), path(:,2), u_true, v_true, q_x_m(1,:),q_y_m(:,1)');
ocean_init = [ocean_x' ; ocean_y'];
V_rel_i =  zeros(1,m-1);
V_abs_i =  zeros(1,m-1);
E_c = zeros(m-1,1);
for i = 2:m-1
        ocean_init(:,i-1) = ocean_init(:,i-1)'*(path(i,1:2)'-path(i-1,1:2)')/norm(path(i,1:2)'-path(i-1,1:2)')*...
                                  (path(i,1:2)'-path(i-1,1:2)')/norm(path(i,1:2)'-path(i-1,1:2)');
       
        V_abs_i(i) = norm((path(i,1:2)'-path(i-1,1:2)')./path(i-1,3));
        V_rel_i(i) = norm((path(i,1:2)'-path(i-1,1:2)')./path(i-1,3)-ocean_init(:,i-1));
        E_c(i) =   (c_d*V_rel_i(i)^3*path(i-1,3));

end
OGD_energy  = sum(E_c);
disp('--------------------------------------------')
disp('--------------------------------------------')
temp6 = ['Energy incurred in OGD path = ', num2str(OGD_energy/1000), ' Kilo Joules'];
disp(temp6)




% figure
% plot(V_rel_i)
% title('relative velocity profile')
% figure
% plot(V_abs_i)
% title('aboslute velocity profile')





end




 












