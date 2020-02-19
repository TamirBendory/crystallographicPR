clear;
close all;
clc;
rng(123);

% This script reproduces Figure 5 from the manuscript "Toward a
% mathematical theory of the crystallographic phase retrieval problem" By
% Tamir Bendory and Dan Edidin

%% parameters

N = 8;
K_vec = 3:4;
num_rep = 1000;
last_iter_XP1 = zeros(length(K_vec), num_rep);

% RRR parameters
max_iter = 1e6;
beta = .5;
stop_criterion = 'error';
th = 1e-8;
verbosity = 0;

%% main loop
for kk = 1:length(K_vec)
    K = K_vec(kk);
    iter = 0;
    while iter<num_rep
        
        %generating the true signal
        ind_true = randperm(N);
        ind_true = ind_true(1:K);
        x_nn = zeros(N,1);
        x_nn(ind_true) = rand(K,1);
        y_nn = abs(fft(x_nn)); % data
        a = ifft(y_nn.^2);
        a = a(1:(floor(N/2)+1));
        if sum(abs(a)>1e-10)> K % consider only cases of |S-S|>K
            iter = iter+1;
            fprintf('N = %g, K = %g, iter = %g\n', N, K, iter);
        else
             fprintf('skip\n');
            continue;
        end
        
        %% RRR
        x_init = rand(N, 1); %random initialization
        [x_est, error, eta, last_iter] = RRR(y_nn, x_init, beta, max_iter, K, stop_criterion, th, x_nn, verbosity);
        last_iter_XP1(kk, iter) = last_iter;
    end
    
    save('last_iter_XP1','last_iter_XP1');
    
end

%% plotting and saving

figure; 
K3 = last_iter_XP1(1,:);
ind = (K3<300);
sum(ind)/length(ind)
histogram(K3(ind), 30);
xlabel('# iterations');
set(gca, 'YScale', 'log')
saveas(gcf,'XP1K3.png')
pdf_print_code(gcf, 'XP1K3', 11);

figure;
K4 = last_iter_XP1(2,:);
ind = (K4<550);
sum(ind)/length(ind)
histogram(K4(ind), 30);
xlabel('# iterations');
set(gca, 'YScale', 'log')
saveas(gcf,'XP1K4.png')
pdf_print_code(gcf, 'XP1K4', 11);


