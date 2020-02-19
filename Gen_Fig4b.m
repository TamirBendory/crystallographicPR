clear;
close all;
clc;
rng(123);

% This script reproduces Figure 4(b) from the manuscript "Toward a
% mathematical theory of the crystallographic phase retrieval problem" By
% Tamir Bendory and Dan Edidin

%% parameters

N = 100; % signal's length
K_vec = 8:24; % sparsity
num_rep = 500;
last_iter_XP2 = zeros(length(K_vec), num_rep);
SS_XP2 = zeros(length(K_vec), num_rep);

% RRR parameters
max_iter = 1e5;
beta = .5;
stop_criterion = 'error';
th = 1e-8;
verbosity = 0;

%% main loop

for kk = 1:length(K_vec)
    K = K_vec(kk);
    for iter = 1: num_rep
        
        fprintf('N = %g, K = %g, iter = %g\n', N, K, iter);
        
        %generating the true signal
        ind_true = randperm(N);
        ind_true = ind_true(1:K);
        x_nn = zeros(N,1);
        x_nn(ind_true) = rand(K,1);
        y_nn = abs(fft(x_nn)); % data
        a = ifft(y_nn.^2);
        a = a(1:(floor(N/2)+1));
        if sum(abs(a)>1e-10)> K %cardinality of S-S
            SS_XP2(kk, iter) = 1;
        end
        
        %% RRR
        x_init = rand(N, 1); %random initialization
        [x_est, error, eta, last_iter] = RRR(y_nn, x_init, beta, max_iter, K, stop_criterion, th, x_nn, verbosity);
        last_iter_XP2(kk, iter) = last_iter;
    end
    save('SS_XP2','SS_XP2');    
    save('last_iter_XP2','last_iter_XP2');
    
end

%% plotting and saving

med_iter = median(last_iter_XP2,2);
ln = 1.2; 
figure; 
hold on;
plot(K_vec, last_iter_XP2(:,1:10), 'b*', 'markersize', 4);
plot(K_vec, med_iter, 'linewidth', ln);
set(gca, 'YScale', 'log')
xlabel('K');
ylabel('# iterations')
ylim([10^1,10^5])
saveas(gcf,'XP2.png')
pdf_print_code(gcf, 'XP2', 11);


