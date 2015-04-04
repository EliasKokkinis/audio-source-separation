%% Central Limit Theorem demonstration
% This is a simple script to demonstrate the CLT

%% Clean-up
clear; close all; clc;

%% Setup
iidSamples = 1000;
nVariables = 1000;
histBins = 30;
plotFor = [1 5 10 50 100 500 1000];

%% Go!
% Initialize variable X
X = zeros(1, iidSamples);

for N = 1 : nVariables
    
    % Generate new random variable
    % An exponential rv is used with mu = 1 (which implies lambda = 1 and
    % variance = 1 also)
    S = exprnd(1, 1, iidSamples);
    
    % Add to X
    X = X + S;
    % Take the average
    Z = X/N;
    % Estimate pdf
    [counts, Xaxis] = hist(Z, histBins);
    pdfX = counts/(iidSamples*(Xaxis(2) - Xaxis(1)));
    
    % Generate corresponding Gaussian pdf 
    % (remember to change values if you change the pdf of S)
    pdfGaussian = pdf('norm', Xaxis, 1, sqrt(1/N));
    
    if ismember(N, plotFor)
        plot(Xaxis, pdfX); hold on;
        plot(Xaxis, pdfGaussian, 'r'); hold off;
        % Decorations/information
        title(sprintf('N = %03d', N));
        xlabel('Random variable Z');
        ylabel('Probability density function');
        legend('Sum of random variables', 'Gaussian');
        grid on;
        drawnow;
    end
end