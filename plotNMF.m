%%
close all;

%% Plot spectrograms
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

%% Plot masks
% Figure 3 plots the component Wiener masks
figure (3);
for k = 1 : K
    subplot(K, 1, k);
    imagesc(10*log10(W(:, k)*H(k,:)));
    title(sprintf('Component mask %02d', k)); xlabel('Frame index'); ylabel('Frequency bin index'); set(gca, 'YDir', 'normal');
    caxis(colorRange);
end
