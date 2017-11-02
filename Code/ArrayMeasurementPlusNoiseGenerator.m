function [x, s, n] = ArrayMeasurementPlusNoiseGenerator(SNR,psi,L)
%% ArrayMeasurementPlusNoiseGenerator
% SNR -> [dB] Signal to noise ratio
% psi -> interelement spacing and source direction
% L   -> Number of elements in array

snr = 10^(SNR/10); % Convert SNR from dB to Absolute

l = (1:L)'-1;

s = exp(-1j*psi*l); % Generate signal values for each element

n = (randn(L,1) + 1j*randn(L,1))./sqrt(2*snr); % Generate complex normal noise values for each element

x = s + n; % construct signal + noise model for each element.