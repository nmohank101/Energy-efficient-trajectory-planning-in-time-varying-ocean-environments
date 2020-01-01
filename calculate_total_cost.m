function cost = calculate_total_cost(R,path,u,v,v_max,siz,X_loc,Y_loc,weight)

 cost = weight*(.5*path(:,1).'*R*path(:,1) + .5 * path(:,2).'*R*path(:,2));                                                   
  
 for i = 2:length(path)-1
     
     temp = cost_with_currents_modified(path(i-1,:),path(i,:),path(i+1,:),u,v,v_max, siz,...
                                                         X_loc,Y_loc);
     
     cost = cost + temp;
     
 end
                                                     
                                                     
                                                     
                                                     
                                                     
end
                                                     
                                                     