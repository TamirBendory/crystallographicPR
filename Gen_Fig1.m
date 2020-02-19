clear;
close all
clc;
rng(40932); %seed

% This script reproduces Figure 1 from the manuscript "Toward a
% mathematical theory of the crystallographic phase retrieval problem" By
% Tamir Bendory and Dan Edidin

%% parameters

N = 500; %signal's length
K_vec = [5, 10, 20, 40]; % sparsity level
num_trial = 10^4; 
support_size = zeros(num_trial, length(K_vec));

%% main loop
for kk = 1:length(K_vec)
    K = K_vec(kk);
    for iter = 1:num_trial
        ind_true = randperm(N);
        ind_true = ind_true(1:K); % random support
        x_true = zeros(N,1);
        x_true(ind_true) = ones(K,1);
        y = abs(fft(x_true)).^2; % data
        a = real(ifft(y));
        a = a(1:(floor(N/2)+1)); %auto-correlation's support
        support_size(iter, kk) = sum(a>0.5);
    end
end

%% plotting and saving

figure;
histogram(support_size(:,1));
xlabel('|S-S|')
saveas(gcf,'XP_SS_K5.png')
pdf_print_code(gcf, 'XP_SS_K5', 11);

figure; 
histogram(support_size(:,2));
xlabel('|S-S|')
saveas(gcf,'XP_SS_K10.png')
pdf_print_code(gcf, 'XP_SS_K10', 11);

figure; 
histogram(support_size(:,3));
xlabel('|S-S|')
saveas(gcf,'XP_SS_K20.png')
pdf_print_code(gcf, 'XP_SS_K20', 11);

figure;
histogram(support_size(:,4));
xlabel('|S-S|')
saveas(gcf,'XP_SS_K40.png')
pdf_print_code(gcf, 'XP_SS_K40', 11);


