clear;
close all
clc;
rng(43232);

% This script reproduces Figure 3 from the manuscript "Toward a
% mathematical theory of the crystallographic phase retrieval problem" By
% Tamir Bendory and Dan Edidin

%% parameters

K_vec = 2:100;
num_trial = 1000;
collision = zeros(length(K_vec),num_trial);

%% main loop
for kk = 1:length(K_vec)
    K = K_vec(kk);
    N = K*100;
    fprintf('K = %.3g\n', K);   
    for iter = 1:num_trial
        ind_true = randperm(N);
        ind_true = ind_true(1:K);
        x_true = zeros(N,1);
        x_true(ind_true) = ones(K,1);
        y = abs(fft(x_true)).^2; % data
        a = real(ifft(y));
        a = a(2:(floor(N/2)+1));
        collision(kk,iter) = sum(a>1.5);
    end
end

%% plotting and saving

collision_mean = mean(collision, 2);
collision_count = (collision==0);
collision_count = sum(collision_count, 2);

figure;
plot(K_vec, collision_mean);
xlabel('K');
%ylabel('average number of collisions');
saveas(gcf,'XP_collisions_mean_fixed_ratio.png')
pdf_print_code(gcf, 'XP_collisions_mean_fixed_ratio', 11);

figure;
stem(K_vec, collision_count);
xlabel('K');
%ylabel('# of collision-free events');
saveas(gcf,'XP_collisions_count_fixed_ratio.png')
pdf_print_code(gcf, 'XP_collisions_count_fixed_ratio', 11);

