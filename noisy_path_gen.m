function [noisy_paths,eps_mat] = noisy_path_gen(N,K,cov_array,decay_it,mag_step,path)

noisy_paths = zeros(K,3,N);
eps_mat = zeros(K,3,N); 

for i = 1:K-1
        means = zeros(2,N);
        temp_eps = mvnrnd(means,cov_array) * decay_it * mag_step;
%         temp_eps(3,:) = temp_eps(3,:) * .1;
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
end
%   keeping the previous optimal trajectory in the current list of samples 
    temp_eps = zeros(3,N).';
    temp_noisy_path = path(:,:) + temp_eps;
    eps_mat(K,:,:) = temp_eps.';
    noisy_paths(K,:,:) = temp_noisy_path.';
end




