function [STLN_path,STLN_path_energy] = generate_stln_path(lat,lon,q_x_m,q_y_m,N,start_point,...
                                                          end_point,v_max,u,v,c_d,per_delay)


%% Creating the initial path
[~ , index_lon] = min(abs(lon - start_point(1)));
[~ , index_lat] = min(abs(lat - start_point(2)));
start_point_m = [q_x_m(index_lat,index_lon), q_y_m(index_lat,index_lon),0];

[~ , index_lon] = min(abs(lon - end_point(1)));
[~ , index_lat] = min(abs(lat - end_point(2)));
end_point_m = [q_x_m(index_lat,index_lon), q_y_m(index_lat,index_lon),0];

% initializing the trajectory
path = zeros(N,3);
path(1,:) = start_point_m;
path(N,:) = end_point_m;

% calculating the straight line path step size
x_step = (path(N,1)-path(1,1)) / (N - 1);
y_step = (path(N,2)-path(1,2)) / (N - 1);

% computing the straight line trajectory 
for i = 2:N-1
    path(i,1) = path(i-1,1) + x_step;
    path(i,2) = path(i-1,2) + y_step;
end

% computing the segment time 
t_init = sqrt(x_step^2 + y_step^2) / v_max ;

t_init = t_init + per_delay * t_init;


path(:,3) = t_init;
path(N,3) = 0;
path_init = path;
STLN_path = path;

% finding the energy incurred in the straight line path
[ocean_x, ocean_y] = vel_ocean_opt(path_init(:,1), path_init(:,2), u, v, q_x_m(1,:),q_y_m(:,1)');
ocean_init = [ocean_x' ; ocean_y'];
V_rel_i =  zeros(1,N);
V_abs_i =  zeros(1,N);
E_c = zeros(N,1);
for i = 2:N
        V_abs_i(i) = norm((path_init(i,1:2)'-path_init(i-1,1:2)')./path_init(i-1,3));
        V_rel_i(i) = norm((path_init(i,1:2)'-path_init(i-1,1:2)')./path_init(i-1,3)-ocean_init(:,i-1));
        E_c(i) =   (c_d*V_rel_i(i)^3*path_init(i-1,3));
end
STLN_path_energy  = sum(E_c);

figure
plot(V_rel_i);
title('relative velocity profile') 

figure
plot(V_abs_i);
title('absolute velocity profile') 

end

















