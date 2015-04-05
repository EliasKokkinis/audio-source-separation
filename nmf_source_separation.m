%% NMF based source separation
clear; close all; clc;
delete('./*.wav');

%% Parameters
inputFile = '../audio/DL7.wav';
% STFT parameters
windowLength = 1024;
hopSize = 256;
analysisWindow  = hamming(windowLength, 'periodic');
synthesisWindow = hanning(windowLength, 'periodic')./hamming(windowLength, 'periodic');

%% Load data
% Read audio file
[x, sampleRate] = audioread(inputFile);
% Analyze into frames
inputFrames = owa(x, windowLength, hopSize, analysisWindow);
% Perform FFT
X = fft(inputFrames);
% Power spectrogram
V = abs(X(1:end/2 + 1, :)).^2;
% Ensure non-negative values
V(V<=0) = 1e-12;

%% NMF
K = 3;
beta = 0;
lambda = 1;
% Perform NMF
[W, H, J] = nmf(V, K, 100, 0, beta, lambda, [], []);

%% Filter and reconstruct
% Estimated spectrogram
Vh = W*H;

for k = 1 : K
   % Pseudo-Wiener mask
   mask = (W(:, k)*H(k, :))./Vh;
   % Apply mask and IFFT
   outputFrames = real(ifft([mask;mask(end-1:-1:2, :)].*X));
   % Synthesize signal
   component(:, k) = ows(outputFrames, hopSize, synthesisWindow);
   % Write audio
   audiowrite(sprintf('cmp%02d.wav', k), component(:, k), sampleRate);
end