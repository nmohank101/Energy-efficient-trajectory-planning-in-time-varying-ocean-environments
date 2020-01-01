% OGD 
function [path, OGD_energy,V_rel_i,V_abs_i] = OGD_dylan_mohan_moving_goal(v_max,N,q_x_m,q_y_m,u,v,mag,dely,lat,lon,c_d,...
                                              start_point,end_point,treshold_conv,...
                                              x_min,x_max,y_min,y_max,u_noisy,v_noisy)



 
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
path(end,1:2) = path(end,1:2)+ [0 2*10^4];
start_loc = path(1,1:2)';    
goal_loc = path(end,1:2)';
curr_loc = start_loc;

OGD_path = zeros(3,N);    % Waypoints that forms the optimal trajectory  
OGD_path(3,:) = t_init;
OGD_path(1:2,1) = start_loc;  % Start location  
m = 2;
norm_param = max(max(mag));

th = 0:.02:4*pi;
r=8000;

fig_h = figure('units','normalized','outerposition',[0 0 1 1]);
hold on;
vid = VideoWriter('dynamic_goal_new.mp4','MPEG-4');
vid.Quality = 95; 
open(vid)

hold on
contourf(q_x_m(y_min:y_max,x_min:x_max),q_y_m(y_min:y_max,x_min:x_max),mag(y_min:y_max,x_min:x_max),'LineColor','none');
caxis([0,max(max(mag))]); colormap (winter); 
htb = colorbar ;
htb.Label.String = 'Magnitude of ocean currents';
htb.FontSize = 24;
htb.FontWeight = 'bold';
quiver(q_x_m(y_min:2:y_max,x_min:2:x_max),q_y_m(y_min:2:y_max,x_min:2:x_max),u(y_min:2:y_max,x_min:2:x_max),v(y_min:2:y_max,x_min:2:x_max),'LineWidth',1,'Color','k');
ylim([5*10^4,9.5*10^4]);
% plot([start_loc(1) goal_loc(1)],[start_loc(2) goal_loc(2)],'--w','LineWidth',3)


% ylim([0,9.11*10^4]);
% xlim([0.65*10^5,2.382*10^5])

%--------------------------------34
% xlim([9.19*10^4,2.382*10^5])
% ylim([2.89*10^4 8.89*10^4]);

h = animatedline('Color','m','LineWidth',3);
g = animatedline('Color','w','LineWidth',3);
h1 = animatedline('Color','m','MaximumNumPoints',1,'Marker','>','MarkerSize',20,'MarkerFaceColor','k');
g_stop = animatedline('Color','w','MaximumNumPoints',1,'Marker','p','MarkerSize',20,'MarkerFaceColor','r');
g_start = animatedline('Color','w','MaximumNumPoints',1,'Marker','p','MarkerSize',20,'MarkerFaceColor','g');
hold on;
title('IOGA trajectory for static goal scenario','FontSize',24,'FontWeight','bold')
xlabel('Distance (m)','FontSize',24,'FontWeight','bold');
ylabel('Distance (m)','FontSize',24,'FontWeight','bold');

ax = gca;
ax.FontWeight = 'bold';
ax.FontSize  = 24;
%--------------------------------34

% scatter(path(1,1), path(1,2),10,'k','o', 'filled');
  
goal_loc(1) = path(end,1) + r*cos(th(1));
goal_loc(2) = path(end,2)+ r*sin(th(1)); 



while( norm(curr_loc-goal_loc) >= treshold_conv ) 
    vel_ocean = find_ocean_vel(OGD_path(1,m-1),OGD_path(2,m-1),u,v, q_x_m(1,:),q_y_m(:,1)'); % finding vel at previous waypoint     
   
    lam_temp =  (1+(((goal_loc-curr_loc)/norm(goal_loc-curr_loc))'*vel_ocean))/2 ; 
    lamda_t =   (1 -  (norm(vel_ocean)/norm_param)* (2*lam_temp-1))*2^-(norm(vel_ocean)/norm_param);
     
    % evaluting the gradient 
    X_grad = (lamda_t*(curr_loc-goal_loc)/norm(curr_loc-goal_loc)) + ((1-lamda_t)*(-(vel_ocean)));
    X_grad = X_grad./t_init;
    
    
    % evaluating the step size
    temp12 = 0.5*(1 + ((goal_loc-curr_loc)'*vel_ocean)/(norm(goal_loc-curr_loc)*norm(vel_ocean)) );
    
    vel_max_modified =   v_max*exp(-0.5*(norm(vel_ocean)/norm_param + sqrt(dely))*sqrt(temp12)); 
    
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
    
%     if m == N && norm(curr_loc-goal_loc) >= treshold_conv
%        x_step = (goal_loc(1)-OGD_path(1,m));
%        y_step = (goal_loc(2)-OGD_path(2,m));
%        [ocean_xm, ocean_ym] = vel_ocean_opt(curr_loc(1), curr_loc(2), u, v, q_x_m(1,:),q_y_m(:,1)');
%        curr_loc(1) = curr_loc(1) + x_step;
%        curr_loc(2) = curr_loc(2) + y_step;
% %        computing the segment time 
%        t_init = sqrt(x_step^2 + y_step^2) / (v_max-norm([ocean_xm, ocean_ym])) ;
%        m = m+1;
%        OGD_path(1:2,m) = curr_loc; 
%        OGD_path(3,m) = 0; 
%        OGD_path(3,m-1) = t_init;
%     end
    
    %----------------------1
    addpoints(g, goal_loc(1), goal_loc(2));
    drawnow;
    addpoints(g_stop,goal_loc(1), goal_loc(2));
    drawnow;
    addpoints(g_start,start_loc(1), start_loc(2));
    drawnow;
    
    addpoints(h, OGD_path(1,m), OGD_path(2,m));
    drawnow;
    addpoints(h1,OGD_path(1,m), OGD_path(2,m));
    drawnow;
    %---------------------1
    
    
    
%     legend({'Ocean velocities','Moving Goal','Energy efficient Trajectory'},'Location','southeast','FontSize',20)
%     title('\fontsize{20}Trajectory Planning in ocean environments')

%       if mod(m,5)==0 || m==2
%     % Capture the plot as an image 
%         frame = getframe(fig_h); 
%         im = frame2im(frame); 
%         [imind,cm] = rgb2ind(im,256);
%       
%         % Write to the GIF File 
%       if m == 2 
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%       else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%       end   
%       end
     
    %------------------2
    m_sb = 0;
    while (m == 2)
        frame = getframe(fig_h); 
        writeVideo(vid,frame)      
        m_sb = m_sb+1;
        if m_sb ==60
         break;
        end 
     end
     
     frame = getframe(fig_h); 
     writeVideo(vid,frame)      
     
     m_sb = 0;
     while (norm(curr_loc-goal_loc) < treshold_conv ) 
           frame = getframe(fig_h); 
           writeVideo(vid,frame)      
           m_sb = m_sb+1;
           if m_sb ==60
               break;
           end
     end
     %----------------2
   m = m + 1 ;

   
% comment below for static goal scenario
   goal_loc(1) = path(end,1) + r*cos(th(m)) ;
   goal_loc(2) = path(end,2) + r*sin(th(m)) ;


end
close(vid);
% scatter(OGD_path(1,m-1), OGD_path(2,m-1),10,'k','o', 'filled');
% scatter(OGD_path(1,1), OGD_path(2,1),10,'k','o', 'filled');


% hold on
% plot(OGD_path(1,1:m-1),OGD_path(2,1:m-1),'k','LineWidth',3)


%% finding the energy incurred in the STOMP trajectory
OGD_path(3,m-1) = 0 ;
path = OGD_path(:,1:m-1)';
[ocean_x, ocean_y] = vel_ocean_opt(path(:,1), path(:,2), u, v, q_x_m(1,:),q_y_m(:,1)');
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
temp6 = ['Energy incurred in OGD path = ', num2str(OGD_energy), ' Joules'];
disp(temp6)




% figure
% plot(V_rel_i)
% title('relative velocity profile')
% figure
% plot(V_abs_i)
% title('aboslute velocity profile')





end




% h = animatedline('Color','m','MaximumNumPoints',1,'Marker','p','MarkerSize',20,'MarkerFaceColor','m');
% g = animatedline('Color','m','LineWidth',3);
% axis([0,4*pi,-1,1])
% x = linspace(0,4*pi,1000);
% y = sin(x);
% for k = 1:length(x)
% addpoints(h,x(k),y(k));
% drawfig_ = figure('units','normalized','outerposition',[0 0 1 1]);now
% addpoints(g,x(k),y(k));
% drawnow
% frame = getframe(fig_);
% writeVideo(v,frame)
% end


%% useful
% fig_ = figure('units','normalized','outerposition',[0 0 1 1]);
% hold on;
% v = VideoWriter('test.mp4','MPEG-4');
% open(v)
% h = animatedline('Color','m','MaximumNumPoints',1,'Marker','p','MarkerSize',20,'MarkerFaceColor','m');
% g = animatedline('Color','m','LineWidth',3);
% axis([0,4*pi,-1,1])
% x = linspace(0,4*pi,1000);
% y = sin(x);
% for k = 1:length(x)
% addpoints(h,x(k),y(k));
% drawnow
% addpoints(g,x(k),y(k));
% drawnow
% frame = getframe(fig_);
% writeVideo(v,frame)
% end
% close(v); 




 