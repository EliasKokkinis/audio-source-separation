function [W, H, J] = nmf(V, K, nIterations, epsilon, beta, lambda, W0, H0)
% This function calculates the NMF of V with rank K based on the
% generalized beta divergence with sparsity parameter lambda.
%
% INPUTS:
%  V             - Input matrix (size FxN)
%  K             - Number of components/Rank of factorization
%  nIterations   - The number of iterations
%  epsilon       - The minimum value of change between iterations
%  beta          - Choose divergence (0 = Itakura-Saito, 1 = Kullback-Leibler, 2 = Euclidean)
%  lambda        - Sparsity parameter for H
%
% OUTPUTS:
%  W             - Spectral profile matrix (size FxK)
%  H             - Activation function matrix (size KxN)
%  J             - Cost function per iteration
%
% Author: Elias Kokkinis
% e-mail: elias@accusonus.com
%
% If you find any bugs, please let me know!

% Extract data dimensions
[F, N] = size(V);

% Random initialization of matrices
if isempty(W0)
    W = abs(randn(F, K));
else
    W = W0;
end

if isempty(H0)
    H = abs(randn(K, N));
else
    H = H0;
end

% Estimate spectrogram
Vh = W*H;
% Initialize vector to hold values of cost function
J = zeros(nIterations, 1);
J(1) = beta_distance(V, Vh, beta);

for i = 2 : nIterations
    % 'Negative' gradient
    Hn = W'*(V.*(Vh.^(beta - 2)));
    % 'Positive' gradient
    Hp = W'*(Vh.^(beta - 1)) + lambda;
    % Multiplicative update rule for H
    H = H.*(Hn./(Hp + 1e-12));
    
    % Estimate spectrogram
    Vh = W*H;
    
    % 'Negative' gradient
    Wn = (V.*(Vh.^(beta - 2)))*H';
    % 'Positive' gradient
    Wp = (Vh.^(beta - 1))*H';
    % Multiplicative update rule for W
    W = W.*(Wn./(Wp + 1e-12));
    
    % Estimate spectrogram
    Vh = W*H;
    
    % Estimate distance and check for convergence
    J(i) = beta_distance(V, Vh, beta);
    if abs(J(i) - J(i - 1)) < epsilon
        break;
    end
end
