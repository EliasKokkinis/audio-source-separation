%% NMF based source separation
clear; close all; clc;
delete('./*.wav');

%% Parameters
inputFile = '../audio/p4.wav';
% STFT parameters
windowLength = 1024;
hopSize = 512;
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
K = 4;
beta = 0;
lambda = 0;
% Perform NMF
[W, H, J] = nmf(V, K, 1, 0, beta, lambda, [], []);

close all;
figure(1);
subplot(2, 1, 1);
imagesc(10*log10(V)); title('Original spectrogram');
xlabel('Frame index'); ylabel('Frequency bin index'); set(gca, 'YDir', 'normal');
colorRange = caxis;

subplot(2, 1, 2);
imagesc(10*log10(W*H)); title('Estimated spectrogram');
xlabel('Frame index'); ylabel('Frequency bin index'); set(gca, 'YDir', 'normal');
caxis(colorRange);

% Figure 2 plots the factorization results (W and H)
figure(2);
subplot(2, 1, 1);
fvec = generateFrequencyVector(windowLength, sampleRate);
semilogx(fvec, W);
title('Spectral profiles'); xlabel('Frequency (Hz)'); xlim([fvec(1) fvec(end)]); ylabel('Magnitude');

subplot(2, 1, 2);
plot(H');
title('Activation functions'); xlabel('Frame index'); ylabel('Gain/Weight');
pause;

for i = 2 : 100
    [W, H, J] = nmf(V, K, 2, 0, beta, lambda, W, H);

    if mod(i, 1) == 0
        
        figure(1);
        subplot(2, 1, 1);
        imagesc(10*log10(V)); title('Original spectrogram');
        xlabel('Frame index'); ylabel('Frequency bin index'); set(gca, 'YDir', 'normal');
        colorRange = caxis;
        
        subplot(2, 1, 2);
        imagesc(10*log10(W*H)); title('Estimated spectrogram');
        xlabel('Frame index'); ylabel('Frequency bin index'); set(gca, 'YDir', 'normal');
        caxis(colorRange);
        
        %% Plot W and H
        % Figure 2 plots the factorization results (W and H)
        figure(2);
        subplot(2, 1, 1);
        fvec = generateFrequencyVector(windowLength, sampleRate);
        semilogx(fvec, W);
        title('Spectral profiles'); xlabel('Frequency (Hz)'); xlim([fvec(1) fvec(end)]); ylabel('Magnitude');
        
        subplot(2, 1, 2);
        plot(H');
        title('Activation functions'); xlabel('Frame index'); ylabel('Gain/Weight');
        pause;
    end
    
        
    disp(sprintf('Iteration %02d\tCost J = %f\n', i, J));
    
end