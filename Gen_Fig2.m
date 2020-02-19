clear;
close all
clc;
rng(43232); % seed

% This script reproduces Figure 2 from the manuscript "Toward a
% mathematical theory of the crystallographic phase retrieval problem" By
% Tamir Bendory and Dan Edidin

%% parameters

N = 1000; % signal's length
K_vec = 2:100; % sparsity
num_trial = 1000;
collision = zeros(length(K_vec), num_trial);

%% main loop

for kk = 1:length(K_vec)
    K = K_vec(kk);
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
collision_mean = mean(collision, 2); % number of collisions
collision_count = (collision==0); % collision-free events
collision_count = sum(collision_count, 2);

figure;
plot(K_vec, collision_mean);
xlabel('K');
%ylabel('average number of collisions');
saveas(gcf,'XP_collisions_mean.png')
pdf_print_code(gcf, 'XP_collisions_mean', 11);

figure;
stem(K_vec, collision_count);
xlabel('K');
%ylabel('# of collision-free events');
saveas(gcf,'XP_collisions_count.png')
pdf_print_code(gcf, 'XP_collisions_count', 11);

