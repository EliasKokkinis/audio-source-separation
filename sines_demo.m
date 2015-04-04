%% 
clear; close all; clc;

%% Prepare data
x = audioread('../audio/sine_demo.wav');
% Analyze into frames
inputFrames = owa(x, 256, 128, hamming(256, 'periodic'));
% Perform FFT
X = fft(inputFrames);
% Power spectrogram
V = abs(X(1:end/2 + 1, :));
% Ensure non-negative values
V(V<=0) = 1e-12;

%% NMF
% Perform NMF
[W, H, J] = nmf(V, 2, 50, 0, 2, 0, [], []);

%% Plot
figure;
imagesc(V); xlabel('Frames index'); ylabel('Frequency bins'); title('Spectrogram');
figure;
plot(W); xlabel('Frequency bins'); title('Spectral profiles');
figure;
plot(H'); xlabel('Frames index'); title('Activation functions');