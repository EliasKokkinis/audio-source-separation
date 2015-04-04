function d = beta_distance(X, Y, beta)
% This function calculates the generalized beta divergence between two
% matrices
%
% INPUTS:
%  X    - A matrix of size MxN
%  Y    - A matrix of size MxN
%  beta - Choose divergence (0 = Itakura-Saito, 1 = Kullback-Leibler, 2 = Euclidean)
%
% OUTPUTS:
%  d    - The scalar value of the distance
%
% Author: Elias Kokkinis
% e-mail: elias@accusonus.com
%
% If you find any bugs, please let me know!
if nargin < 3
    error('Three input arguments are required!');
end

if (size(X, 1)~=size(Y, 1))||(size(X, 2)~=size(Y, 2))
    error('Matrix dimensions do not match!');
end

switch beta
    case 0
        d = sum(sum(X./Y-log(X./Y)-1));
    case 1
        d = sum(sum(X.*log(X./Y) + (Y - X)));
    otherwise
        d = sum(sum((1/(beta*(beta - 1)))*(X.^beta + (beta - 1).*(Y.^beta) - beta.*X.*(Y.^(beta - 1)))));
end