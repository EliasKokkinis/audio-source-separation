function W = fastICA(X, nonLin, maxIter, epsilon)
% This function implementes the FastICA method
%
% INPUTS
%   X       - The preprocessed input data (N x samples)
%   nonLin  - The non-linearity/the measure to maximize ('kurtosis' or
%   'negentropy'
%   maxIter - The maximum number of iterations
%   epsilon - The minimum change between iterations
%
% OUTPUTS
%   W       - The unmixing matrix. It needs post-processing
%
% Author: Elias Kokkinis
% e-mail: elias@accusonus.com
%
% ...

% Extract dimension
[N, samples] = size(X);
% Random initialization
W = orth(randn(N));
% Previous estimation
W_old = zeros(size(W));

a1 = 1;

for i = 1 : maxIter
    % Symmetric orthogonalization.
    W = W * real(inv(W' * W)^(1/2));
    
    minAbsCos = min(abs(diag(W' * W_old)));
    if (1 - minAbsCos < epsilon)
        fprintf('FastICA converged after %d iterations!\n', i);
        break;
    end
    
    W_old = W;
    
    switch nonLin
        case 'kurtosis'
            W = X*(X'*W).^3 - 3*W;
        case 'negentropy'
            hypTan = tanh(a1 * X' * W);
            W = X*hypTan/samples - ones(N, 1)*sum(1 - hypTan.^2).*W/samples*a1;
        otherwise
            error('Unsupported non-linearity');
    end
end

if i == maxIter
    fprintf('FastICA reached maximum number of iterations!\n');
end