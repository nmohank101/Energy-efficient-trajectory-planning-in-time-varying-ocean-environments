% Function to do Cost Function with Currents


function [cost] = cost_with_currents_modified(w1, w2, w3, u, v, v_max, siz,X_loc,Y_loc) %, gp_u, gp_v)
                                   
%function [cost] = cost_with_currents_expectation_test(w1, w2, w3, v_max, size)

    % Constants for calculations
    cost = 0;
    L = 10;
    factor = 1;
    c_d = 3;
    
    l = size(u);
    
    % Constant for object checking
    O = 1000;
    
    % Making the points
    x_1 = w1(1);
    y_1 = w1(2);
    t_1 = w1(3);
    
    x_2 = w2(1);
    y_2 = w2(2);
    t_2 = w2(3);
    
    x_3 = w3(1);
    y_3 = w3(2);
    
    % Calculating distances between points
    dist_12 = sqrt((x_1 - x_2)^2 + (y_1 - y_2)^2);
    dist_23 = sqrt((x_2 - x_3)^2 + (y_2 - y_3)^2);
    
    t_13 = t_1 + dist_23 / dist_12 * t_1;
   
    % Translating the waypoint location to our discritized current field
    % Need to add bound checking here
    x_index_1 = abs(int8(l(1) / siz(2) * x_1)) + 1;
    y_index_1 = abs(int8(l(2) / siz(1) * y_1))+ 1;
    
    x_index_2 = abs(int8(l(1) / siz(2) * x_2)) + 1;
    y_index_2 = abs(int8(l(2) / siz(1) * y_2)) + 1;
    
    x_index_3 = abs(int8(l(1) / siz(2) * x_3)) + 1;
    y_index_3 = abs(int8(l(2) / siz(1) * y_3)) + 1;
    try
    if (isnan(u(x_index_1,y_index_1)) || isnan(u(x_index_2, y_index_2)) ||...
        isnan(u(x_index_3, y_index_3)) || isnan(v(x_index_1, y_index_1)) ||...
        isnan(v(x_index_2, y_index_2)) || isnan(v(x_index_3, y_index_3)))
        
    % Should make this a function of how deep it is into the object
        cost = O^2;
    else
    
        % Calculating the average current between the two way points
        vel_ocean_temp1 = find_ocean_vel(x_1,y_1,u,v,X_loc,Y_loc);
        vel_ocean_temp2 = find_ocean_vel(x_2,y_2,u,v,X_loc,Y_loc);
        vel_ocean_temp3 = find_ocean_vel(x_3,y_3,u,v,X_loc,Y_loc);
        
        u_avg_12 = (vel_ocean_temp1(1) + vel_ocean_temp2(1))/2 ;
        v_avg_12 = (vel_ocean_temp1(2) + vel_ocean_temp2(2))/2 ;
        
        u_avg_23 = (vel_ocean_temp2(1) + vel_ocean_temp3(1))/2 ;
        v_avg_23 = (vel_ocean_temp2(2) + vel_ocean_temp3(2))/2 ;
         
        u_avg_13 = (vel_ocean_temp1(1) + vel_ocean_temp3(1))/2 ;
        v_avg_13 = (vel_ocean_temp1(2) + vel_ocean_temp3(2))/2 ;
        
        
        
        
        
        % Calculating the required velocity to travel between the points
        u_req_12 = (x_2 - x_1) / t_1 - u_avg_12;
        v_req_12 = (y_2 - y_1) / t_1 - v_avg_12;
        
        u_req_23 = (x_3 - x_2) / t_2 - u_avg_23;
        v_req_23 = (y_3 - y_2) / t_2 - v_avg_23;
          
        u_req_13 = (x_3 - x_1) / t_13 - u_avg_13;
        v_req_13 = (y_3 - y_1) / t_13 - v_avg_13;
        
        %abs_vel_12 = sqrt(((x_2 - x_1) / t_1)^2 + ((y_2 - y_1) / t_1)^2);
        abs_vel_23 = sqrt(((x_3 - x_2) / t_2)^2 + ((y_3 - y_2) / t_2)^2);
        %abs_vel_13 = sqrt(((x_3 - x_1) / t_13)^2 + ((y_3 - y_1) / t_13)^2);
        
        vel_req_12 = sqrt(u_req_12^2 + v_req_12^2);
        vel_req_23 = sqrt(u_req_23^2 + v_req_23^2);
        vel_req_13 = sqrt(u_req_13^2 + v_req_13^2);

        
        
        
        
        
        
        if (vel_req_23 <= v_max && abs_vel_23 >= v_max)
            %display('option 1')
            cost = 0;
        elseif (vel_req_23 <= v_max)
            %display('option 2')
            cost = cost + factor * exp(-(vel_req_23 - v_max));
        else
            %display('option 3')
            cost = cost + L + (vel_req_23 - v_max) * L^2;
        end

        % Adding in an energy cost
        % cost = cost + c_d * vel_req ^ 3 * t_1;
        
        % Energy cost with waypoint
        cost_with = c_d * vel_req_12 ^ 3 * t_1 + c_d * vel_req_23 ^ 3 * t_2;

        % Energy cost without waypoint
        cost_without = c_d * vel_req_13 ^ 3 * t_13;
    
        % Scaling the difference
        diff = cost_with - cost_without;
        
        if diff ~= 0
            diff = diff / (10 ^ floor(log10(abs(diff))));
        end
        
        % Adding in energy cost
        cost = cost + exp(diff);
    end
    catch
        keyboard
    end
end