function [STOMP_path, STOMP_energy,cost_STOMP,V_rel_i,V_abs_i] = STOMP_dylan_mohan(v_max,K,num_its,decay_fact,...
                                          N,threshold,mag_step,q_x_m,q_y_m,u,v,mag,path,c_d,...
                                          x_min,x_max,y_min,y_max,u_true,v_true)


 
%% Initializing required parameters
h = 10;                   
weight = 0.001;           
decay_it = decay_fact;
change_factor = .0005;

%% Initial calculations for STOMP algorithm

D = finitediff1(N);
Dinv = inv(D);
DinvTran = inv(D.');

T = eye(N);
Tinv = T;

% Handling the two end cases
T(1,1) = 1 / (path(1,3) ^ 2);
Tinv(1,1) = path(1,3) ^ 2;
T(N,N) = 1 / (path(N-1,3) ^ 2);
Tinv(N,N) = path(N-1,3) ^ 2;

% Taking the square of the average time for the T matrix 
for i = 2:N-1
    T(i,i) = 1 / (((path(i-1,3) + path(i,3)) / 2) ^ 2);
    Tinv(i,i) = (((path(i-1,3) + path(i,3)) / 2) ^ 2);
end

%% Running STOMP 

% creating R matrix
R = D.'* T * T * D;
cost_STOMP = [];
m = 1;
flag = true;

% run the following loop untill convergence
while (flag && m <= num_its)
    
    % creating the covariance matrix
    cov_mat = Dinv * Tinv * Tinv * DinvTran;
    
    % this is needed to handle numerical errors 
    cov_mat = (cov_mat + cov_mat.') / 2;
    
    % creating the scaling matrix for position update
    scale_M = max(cov_mat) * N;
    M = zeros(N,N);
    for i = 1:N
        M(:,i) = cov_mat(:,i) / scale_M(i);
    end
    
    % Array for creating the pertubations of the path
    cov_array(:,:,1) = cov_mat;
    cov_array(:,:,2) = cov_mat; 
    
    % variable to hold all noisy paths generated
    [noisy_paths,eps_mat] = noisy_path_gen(N,K,cov_array,decay_it,mag_step,path);
    
    % visualizing all the noisy trajectories around the current best trajectory
    figure(4)
    clf
    contourf(q_x_m(y_min:y_max,x_min:x_max),q_y_m(y_min:y_max,x_min:x_max),mag(y_min:y_max,x_min:x_max),'LineColor','none');
    caxis([0,max(max(mag))]); 
    colormap (jet); 
    htb = colorbar ;
    htb.Label.String = 'Magnitude of ocean currents';
    hold on;
    quiver(q_x_m(y_min:y_max,x_min:x_max),q_y_m(y_min:y_max,x_min:x_max),u(y_min:y_max,x_min:x_max),v(y_min:y_max,x_min:x_max),'LineWidth',1,'Color','k');
    plot(path(:,1),path(:,2),'r-x')
    x = zeros(N,1);
    y = x;
    for i = 1:K
        for j = 1:N
            x(j) = noisy_paths(i,1,j);
            y(j) = noisy_paths(i,2,j);
        end
        plot(x,y,'g-x')
    end
    plot(path(:,1),path(:,2),'r-x')
    hold off

    % loop to calculate cost for all the noisy trajectories
    cost_mat = zeros(N,K);
    for i = 1:K
        for j = 2:N-1
            waypoint1 = zeros(1,3);
            waypoint2 = zeros(1,3);
            waypoint3 = zeros(1,3);

            waypoint1(1) = noisy_paths(i,1,j-1);
            waypoint1(2) = noisy_paths(i,2,j-1);
            waypoint1(3) = noisy_paths(i,3,j-1);
            
            waypoint2(1) = noisy_paths(i,1,j);
            waypoint2(2) = noisy_paths(i,2,j);
            waypoint2(3) = noisy_paths(i,3,j);
            
            waypoint3(1) = noisy_paths(i,1,j+1);
            waypoint3(2) = noisy_paths(i,2,j+1);
            waypoint3(3) = noisy_paths(i,3,j+1);
           
            cost_mat(j,i) = cost_with_currents_modified(waypoint1,waypoint2,waypoint3,...
                                                      u,v,v_max, [max(max(q_y_m)),max(max(q_x_m))],...
                                                      q_x_m(1,:),q_y_m(:,1)');
        end   
    end
    
    % Finding the minimum and maximum costs for each watypoint number
    max_costs = max(cost_mat.');
    min_costs = min(cost_mat.');
    
    % Finding the probability contribution of each waypoint
    prob_mat = zeros(N,K);
    e_mat = zeros(N,K);
    for i = 1:N-1
        e_sum = 0;
        if max_costs(i) ~= min_costs(i)
            for j = 1:K
                e_mat(i,j) = exp(-h*((cost_mat(i,j) - min_costs(i))/(max_costs(i) - min_costs(i))));
                e_sum = e_sum + e_mat(i,j);
            end
        
            prob_mat(i,:) = e_mat(i,:) / e_sum;
        else
            prob_mat(i,:) = 1 / K;
        end
    end
    
    % updating the optimal trajectory
    update_vector = zeros(N,3);
    % Loop through and add all contributions into update vector
    for i = 1:N-1
        for j = 1:K
            update_vector(i,:) = update_vector(i,:) + eps_mat(j,:,i) * prob_mat(i,j);
        end
    end
    % Scaling the update vector and ensuring the endpoints do not move
    update_vector(:,1) = M * update_vector(:,1);
    update_vector(:,2) = M * update_vector(:,2);
    update_vector(1,:) = [0,0,update_vector(1,3)];
    update_vector(N,:) = [0,0,0];
    
    % updating the trajectory
    new_path = path(:,:) + update_vector;
    
    % calculating the cost of the new optimal trajectory
    cost_new = calculate_total_cost(R,new_path,u,v,v_max,[max(max(q_y_m)),max(max(q_x_m))],...
                                                         q_x_m(1,:),q_y_m(:,1)',weight);
    cost_STOMP = [cost_STOMP cost_new];
    
    if m > 1
        if (~((cost_STOMP(m-1) - cost_STOMP(m)) > (change_factor * cost_STOMP(m))) && (cost_STOMP(m) < threshold))
            flag = false;
            display('stopping for cost reasons')
        end
    end
    
    % updating the parameters
    path(:,:) = new_path;
    m = m + 1;
    decay_it = decay_it * decay_fact;
    pause(.25)
end 


%% finding the energy incurred in the STOMP trajectory
[ocean_x, ocean_y] = vel_ocean_opt(path(:,1), path(:,2), u, v, q_x_m(1,:),q_y_m(:,1)');
ocean_init = [ocean_x' ; ocean_y'];
V_rel_i =  zeros(1,N);
V_abs_i =  zeros(1,N);
E_c = zeros(N,1);
for i = 2:N
     ocean_init(:,i-1) = ocean_init(:,i-1)'*(path(i,1:2)'-path(i-1,1:2)')/norm(path(i,1:2)'-path(i-1,1:2)')*...
                                  (path(i,1:2)'-path(i-1,1:2)')/norm(path(i,1:2)'-path(i-1,1:2)');
       
        V_abs_i(i) = norm((path(i,1:2)'-path(i-1,1:2)')./path(i-1,3));
        V_rel_i(i) = norm((path(i,1:2)'-path(i-1,1:2)')./path(i-1,3)-ocean_init(:,i-1));
        E_c(i) =   (c_d*V_rel_i(i)^3*path(i-1,3));
end
energy_init  = sum(E_c);
disp('--------------------------------------------')
disp('--------------------------------------------')
temp6 = ['Energy incurred in STOMP path = ', num2str(energy_init/1000), ' Kilo Joules'];
disp(temp6)

 
STOMP_path  = path;
STOMP_energy = sum(E_c);

% figure
% plot(V_rel_i)
% title('relative velocity profile') 
% figure
% plot(V_abs_i)
% title('absolute velocity profile')

end

























