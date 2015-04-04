function W = fastICA(X, measure, maxIter, epsilon)
% This function implementes the FastICA method
%
% INPUTS
%   X       - The preprocessed input data (N x samples)
%   measure - The measure to maximize ('kurtosis' or 'negentropy')
%   maxIter - The maximum number of iterations
%   epsilon - The minimum change between iterations
%
% OUTPUTS
%   W       - The unmixing matrix. It needs post-processing
%
% Author: Elias Kokkinis
% e-mail: elias@accusonus.com
%
% If you find any bugs, please let me know!

% Extract dimension
[N, samples] = size(X);
% Random initialization
W = randn(N);
% Previous estimation
W_old = zeros(size(W));

for i = 1 : maxIter
    % Symmetric orthogonalization.
    W = W * real(inv(W' * W)^(1/2));
    
    minAbsCos = min(abs(diag(W' * W_old)));
    if (1 - minAbsCos < epsilon)
        fprintf('FastICA converged after %d iterations!\n', i);
        break;
    end
    
    W_old = W;
    
    switch measure
        case 'kurtosis'
            W = X*(X'*W).^3 - 3*(ones(N, 1)*mean((X'*W).^2))*W;
        case 'negentropy'
            hypTan = tanh(X'*W);
            W = X*hypTan/samples - ones(N, 1)*sum(1 - hypTan.^2).*W/samples;
        otherwise
            error('Unsupported non-linearity');
    end
end

if i == maxIter
    fprintf('FastICA reached maximum number of iterations!\n');
end