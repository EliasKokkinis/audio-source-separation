%% Clean-up
clear; close all; clc;
load ICA_dataset_3.mat

%% Pre-processing
% Extract dimensions
[N, samples] = size(X);

% Remove mean
M = repmat(mean(X, 2), [1 samples]);
Xn = X - M;

% Covariance matrix
C = cov(Xn');
% EVD
[E, D] = eig(C);
% Calculate whitening matrix
sqrtD = diag(sqrt(diag(D)));
Tw = inv(sqrtD)*E';
Td = E*sqrtD;

% Whiten the data
Z = Tw*Xn;

%% FastICA
B = fastICA(Z, 'negentropy', 100, 1e-6);

%% Post-processing
% Unmixing matrix
W = B'*Tw;
% Separated signals (don't forget to add back the mean!)
Y = W*X + (W*M);