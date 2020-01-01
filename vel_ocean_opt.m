function [vel_u, vel_v] = vel_ocean_opt(theta_x,theta_y,U,V,X_loc,Y_loc)
%
N = size(theta_x,1);
K = size(theta_y,2);

vel_u = zeros(N,K);
vel_v = zeros(N,K);

for i = 1:N
    for j = 1:K
        way_point = [theta_x(i,j);theta_y(i,j)];      
        x_limit = [find(X_loc <= way_point(1), 1, 'last' ), find(X_loc >= way_point(1), 1 )];
        y_limit = [find(Y_loc <= way_point(2), 1, 'last' ), find(Y_loc >= way_point(2), 1 )];
        if (length(x_limit)==2 && length(y_limit)==2) 
            vel_u(i,j) = 0.25*(U(y_limit(1),x_limit(1)) + U(y_limit(2),x_limit(1)) + U(y_limit(1),x_limit(2)) + U(y_limit(2),x_limit(2))) ;
            vel_v(i,j) = 0.25*(V(y_limit(1),x_limit(1)) + V(y_limit(2),x_limit(1)) + V(y_limit(1),x_limit(2)) + V(y_limit(2),x_limit(2))) ;                   
        elseif (length(x_limit)==1 && length(y_limit)==2)
            vel_u(i,j) = 0.5*(U(y_limit(1),x_limit(1)) + U(y_limit(2),x_limit(1))) ;
            vel_v(i,j) = 0.5*(V(y_limit(1),x_limit(1)) + V(y_limit(2),x_limit(1))) ;
        elseif (length(x_limit)==2 && length(y_limit)==1)           
            vel_u(i,j) = 0.5*(U(y_limit(1),x_limit(1)) + U(y_limit(1),x_limit(2))) ;
            vel_v(i,j) = 0.5*(V(y_limit(1),x_limit(1)) + V(y_limit(1),x_limit(2))) ;
        elseif (way_point(1) < min(X_loc)) || (way_point(1) > max(X_loc)) || (way_point(2) < min(Y_loc)) || (way_point(2) > max(Y_loc)) 
            vel_u(i,j) = inf;
            vel_v(i,j) = inf;
        end
    end
end