% Code for generating a finite difference matrix using forward and
% backwards differentiation rules


function [A] = finitediff1(n)
    
    A = zeros(n,n);
    
    A(1,1) = -2;
    A(n,n) = -2;
    
    A(1,2) = 1;
    A(n,n-1) = 1;
    
    for j = 2:(n / 2)
        A(j,j-1) = 1;
        A(j,j) = -2;
        A(j,j+1) = 1;
        
        A(n-j+1,n-j+2) = 1;
        A(n-j+1,n-j+1) = -2;
        A(n-j+1,n-j) = 1;
    end
    
    if mod(n, 2) == 1
        index = int8(n / 2);
        A(index,index-1) = 1;
        A(index,index) = -2;
        A(index,index+1) = 1;
    end
    
end