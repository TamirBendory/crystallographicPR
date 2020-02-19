function [x_est, error, eta, last_iter] = RRR(y, x_init, beta, max_iter, K, stop_criterion, th, x_true, verbosity)
% 
% The RRR algorithm
% 
% Inputs:
% 
% y: the observed data
% x_init: initial guess of the signal 
% beta: RRR parameter (step size)
% max_iter: maxmimum number of iterations for RRR
% K: expected sparsity level
% stop_criterion: the algorithm halts either when eta is sufficiently
% large (this measure can be computed in practice) or when the error
% compared to the ground truth is small (this measure cannot be computed in
% practice)
% th: threshold to halt RRR iterations. If the iterations do not
% attain this value, then the algorithm halts after max_itet iterations.
% x_true: the true signal; used *only* to compute the error*if the stopping
% criterion is the error. Not used otherwise. 
% verbosity: if 1, reports on the iterations progress
%
%
% Output: 
%
% x_est - signal's estimate
% error - the error compared with the ground truth (for all iterations)
% eta - the index eta (for all iterations)
% last_iter - number of RRR iterations 
%
% Tamir Bendory
% Last update: Feb 19 2020
% 

L = length(y);

% initialization
if ~exist('x_init', 'var') || isempty(x_init)
    x_est = randn(L,1);
else
    x_est = x_init; 
end
% default values
if ~exist('max_iter', 'var') || isempty(max_iter)
    max_iter = 1e4;
end
if ~exist('th', 'var') || isempty(th)
    th = 1e-8;
end

eta = zeros(max_iter,1);
error = zeros(max_iter,1);
last_iter = max_iter;

% main loop
for iter = 1:max_iter
    
    % one RRR iteration
    x1 =  P1(x_est, K);
    x2 = P2(2*x1 - x_est, y);
    x_est = x_est + beta*(x2-x1);
    
    % measuring objectives
    x_proj = P1(x_est, K);
    eta(iter) = norm(x_proj)^2./norm(x_est)^2;
    error(iter) = compute_error(x_proj, x_true);
    % reporting progress    
    if mod(iter, max_iter/1000) == 0 &&  verbosity == 1
        fprintf('iter = %g, eta = %.4g, error = %.4g\n', iter, eta(iter), error(iter));
    end
    
    if strcmp(stop_criterion, 'error')
        if error(iter)<th
            last_iter = iter;
            eta = eta(1:last_iter);
            error = error(1:last_iter);
            break;
        end
    else % stopping by eta
        if eta(iter)>th
            last_iter = iter;
            eta = eta(1:last_iter);
            error = error(1:last_iter);
            break;
        end
    end
end

% output 
x_est = P1(x_est, K);

end

%% Auxiliary functions

function x1 = P1(x, K)
% Projecting x onto the significant K entries
x1 = zeros(size(x));
[val, ind] = maxk(x, K, 'ComparisonMethod','abs');
x1(ind) = val;
end

function x2 = P2(x, y)
% projecting x onto the Fourier magnitude of y
x2 = ifft(y.*sign(fft(x)));
end

function err = compute_error(x_est, x_true)
% comparing the error between the two signals, while taking all symmetries
% into account
X_true = fft(x_true);
X_est = fft(x_est);
a1 = abs(ifft(X_true.*X_est)); % the absolute values take the sign symmetry into account
a2 = abs(ifft(X_true.*conj(X_est))); % the reflected signal
max_correlation = max([a1; a2]);
err = norm(x_est).^2 + norm(x_true).^2 - 2*max_correlation;
err = err/norm(x_true).^2;
end
