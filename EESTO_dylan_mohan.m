function [EESTO_path, EESTO_energy,cost_EESTO,V_rel_i,V_abs_i] = EESTO_dylan_mohan(v_max,K,num_its,decay_fact,...
                                          N,threshold,mag_step,q_x_m,q_y_m,u,v,mag,path,c_d,...
                                          x_min,x_max,y_min,y_max,u_true,v_true)

                                      
%% parameters   
h = 10;
weight = 0.001; 
decay_it = decay_fact;
change_factor = .0005;

num_start = 10;

%% Initial calculations for EESTO algorithm

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

%% Running EESTO

% Creating variables for plotting puposes
tot_cost = zeros(num_start,1);
m = 1;
flag = true;
cost_EESTO = [];

% run the following loop untill convergence
while (flag && m <= num_its)
    
 
    % creating the covariance matrix
    cov_mat = Dinv * Tinv * Tinv * DinvTran;

    % this is needed to handle numerical errors
    cov_mat = (cov_mat + cov_mat.') / 2;
    
    % Creating the scaling matrix for position update
    scale_M = max(cov_mat) * N;
    M = zeros(N,N);
    for i = 1:N
        M(:,i) = cov_mat(:,i) / scale_M(i);
    end
    
    % Creating the scaling matrix for time update
    scale_T = max(Tinv)  * N;  %* .1
    M_t = zeros(N,N);
    for i = 1:N
        M_t(:,i) = Tinv(:,i) / scale_T(i);
    end
    
    % Array for creating the pertubations of the path
    cov_array(:,:,1) = cov_mat;
    cov_array(:,:,2) = cov_mat;
    cov_array(:,:,3) = Tinv;
    
    % Variable to hold all noisy paths generated
    noisy_paths = zeros(K,3,N);
    eps_mat = zeros(K,3,N);
    
    % Generating the pertubations to the initial path
    for i = 1:K-1
        means = zeros(3,N);
        temp_eps = mvnrnd(means,cov_array) * decay_it * mag_step;
        temp_eps(3,:) = temp_eps(3,:) ; % * .1
        temp_eps = temp_eps.';
        % Ensuring that start and end goals do not move
        temp_eps(1,1) = 0;
        temp_eps(1,2) = 0;
        temp_eps(N,1) = 0;
        temp_eps(N,2) = 0;
        temp_eps(N,3) = 0;
        temp_noisy_path = path(:,:) + temp_eps;
        eps_mat(i,:,:) = temp_eps.';
        noisy_paths(i,:,:) = temp_noisy_path.';
        % Looping to ensure that no paths have a negative travel time
        for j = 1:N-1
            if (noisy_paths(i,3,j) <= 0)
                x_dist = noisy_paths(i,1,j) - noisy_paths(i,1,j+1);
                y_dist = noisy_paths(i,2,j) - noisy_paths(i,2,j+1);
                noisy_paths(i,3,j) = 2 * sqrt(x_dist^2 + y_dist^2) / v_max;               
            end
        end
    end
    temp_eps = zeros(3,N).';
    temp_noisy_path = path(:,:) + temp_eps;
    eps_mat(K,:,:) = temp_eps.';
    noisy_paths(K,:,:) = temp_noisy_path.';

    % visualizing all the noisy trajectories around the current best trajectory
    figure(5)
    clf
    set(gca,'Color',[0.8 0.8 0.8]);
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
    plot(path(:,1),path(:,2),'w-x')
    hold off
    
    % Loop to calculate all the costs for the noisy paths - Second cost
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
    update_vector(:,3) = M_t * update_vector(:,3);
    update_vector(1,:) = [0,0,update_vector(1,3)];
    update_vector(N,:) = [0,0,0];
    
    % update the trajectory
    new_path = path(:,:) + update_vector;
    path(:,:) = new_path;
    
    % Handling the two end cases
    T(1,1) = 1 / (path(1,3) ^ 2);
    Tinv(1,1) = path(1,3) ^ 2;
    T(N,N) = 1 / (path(N-1,3) ^ 2);
    Tinv(N,N) = path(N-1,3) ^ 2;

    % Taking the square of the average time 
    for i = 2:N-1
        T(i,i) = 1 / (((path(i-1,3) + path(i,3)) / 2) ^ 2);
        Tinv(i,i) = (((path(i-1,3) + path(i,3)) / 2) ^ 2);
    end
    
    
    % Creating R matrix
    R = D.'* T * T * D;
    
    % calculating the cost of the new optimal trajectory
    tot_cost(m) = weight * (.5 * path(:,1).'*R*path(:,1) + .5 * path(:,2).'*R*path(:,2));
    for j = 2:N-1
        waypoint1 = zeros(1,3);
        waypoint2 = zeros(1,3);
        waypoint3 = zeros(1,3);

        waypoint1(1) = path(j-1,1);
        waypoint1(2) = path(j-1,2);
        waypoint1(3) = path(j-1,3);

        waypoint2(1) = path(j,1);
        waypoint2(2) = path(j,2);
        waypoint2(3) = path(j,3);

        waypoint3(1) = path(j+1,1);
        waypoint3(2) = path(j+1,2);
        waypoint3(3) = path(j+1,3);

        tot_cost(m) = tot_cost(m) + cost_with_currents_modified(waypoint1,waypoint2,waypoint3,...
                                                         u,v,v_max, [max(max(q_y_m)),max(max(q_x_m))],...
                                                         q_x_m(1,:),q_y_m(:,1)');  
    end   
    
    % Updating the decay factor
    decay_it = decay_it * decay_fact;
    
    cost_EESTO = [cost_EESTO tot_cost(m)];
        
    if m > 1
        if (~((tot_cost(m-1) - tot_cost(m)) > (change_factor * tot_cost(m))) && (tot_cost(m) < threshold))
            flag = false;
            display('stopping for cost reasons')
        end
    end
    m = m + 1;
    pause(.5);
end


%% finding the energy incurred in the STOMP trajectory
[ocean_x, ocean_y] = vel_ocean_opt(path(:,1), path(:,2), u, v, q_x_m(1,:),q_y_m(:,1)');
ocean_init = [ocean_x' ; ocean_y'];
V_rel_i =  zeros(1,N);
V_abs_i =  zeros(1,N);
E_c = zeros(N,1);
for i = 2:N
        V_abs_i(i) = norm((path(i,1:2)'-path(i-1,1:2)')./path(i-1,3));
        V_rel_i(i) = norm((path(i,1:2)'-path(i-1,1:2)')./path(i-1,3)-ocean_init(:,i-1));
        E_c(i) =   (c_d*V_rel_i(i)^3*path(i-1,3));
end
energy_init  = sum(E_c);
disp('--------------------------------------------')
disp('--------------------------------------------')
temp6 = ['Energy incurred in EESTO path = ', num2str(energy_init/1000), ' Kilo Joules'];
disp(temp6)

 
EESTO_path  = path;
EESTO_energy = sum(E_c);

% figure
% plot(V_rel_i)
% title('relative velocity profile') 
% figure
% plot(V_abs_i)
% title('absolute velocity profile')


end











