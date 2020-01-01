%% to find a straight line path with out incorporating delay 

function [path_init,N_init,STLN_path_energy] =  generate_straight_line(start_loc,goal_loc,u,v,X_loc,...
                                                        Y_loc,vel_max,min_time,c_d,min_dist)
% start and end points
X_start = start_loc(1:2)';
X_end = goal_loc(1:2)';

% to store straight line trajectory 
opt_traj = X_start;


vel_max_modified =  vel_max;
m_ = 0;

% while loop till it reaches end location
while (norm(X_start-X_end) > min_dist)
    m_ = m_ + 1; 
    vel_ocean = find_ocean_vel(X_start(1),X_start(2),u,v,X_loc(1,:),Y_loc(:,1)'); 
    if isnan(vel_ocean) 
        vel_ocean = zeros(2,1);
    end
    
% moving in straight line direction
    X_grad = (X_start - X_end)/norm(X_start - X_end);
    X_grad = X_grad/min_time;
    
% calculating step size
    temp7 = 4*(X_grad'*vel_ocean)^2 - 4*(norm(vel_ocean)^2-vel_max_modified ^2)*norm(X_grad)^2; 
    ss = ((-2*X_grad'*vel_ocean)-sqrt(temp7))/(2*(norm(vel_ocean)^2-vel_max_modified ^2));
    if ss < 0
        ss = ((-2*X_grad'*vel_ocean)+sqrt(temp7))/(2*(norm(vel_ocean)^2-vel_max_modified ^2));
    end
    ss = 1/ss;
    X_grad = X_grad*min_time;
    
% performing online gradient descent update     
    X_start = X_start - ss*X_grad;
    
    opt_traj = [opt_traj X_start];
    
end

opt_traj = [opt_traj X_end];
opt_traj(3,:) = min_time;
opt_traj(3,end) = 0;


% total number of waypoints including start and end locations 
N_init = size(opt_traj,2);


% Now finding the energy incurred by travelling in straight line path
path_init = opt_traj';
[ocean_x, ocean_y] = vel_ocean_opt(path_init(:,1), path_init(:,2), u, v, X_loc(1,:),Y_loc(:,1)');
ocean_init = [ocean_x' ; ocean_y'];
V_rel_i =  zeros(1,N_init);
V_abs_i =  zeros(1,N_init);
E_c = zeros(N_init,1);
for i = 2:N_init
        V_abs_i(i) = norm((path_init(i,1:2)'-path_init(i-1,1:2)')./path_init(i-1,3));
        V_rel_i(i) = norm((path_init(i,1:2)'-path_init(i-1,1:2)')./path_init(i-1,3)-ocean_init(:,i-1));
        E_c(i) =   (c_d*V_rel_i(i)^3*path_init(i-1,3));
end
idx = ~isnan(E_c);
STLN_path_energy  = sum(E_c(idx));

figure
plot(V_rel_i);
title('relative velocity of the vehicle travelling in st.line path')
xlabel('segement');
ylabel('m/sec');
figure
plot(V_abs_i);
title('absolute velocity of the vehicle travelling in st.line path');
xlabel('segment');
ylabel('m/sec');




end