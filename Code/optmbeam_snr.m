format compact; clear; close all; clc;
for L = 16 % Number of elements
    for N = [10 100 10000] % Number of samples
        for SNR_dB = [3 10] % [dB] Signal to Noise Ratio
            for th_s = 30 % [deg] Source direction from normal of the array
                %%
                SNR = 10^(SNR_dB/10); % Absolute Signal to Noise Ratio
                psi = pi*sind(th_s); % Phase difference between elements
                
                wc = (1/L)*exp(-1j*psi*((1:L)-1)'); % convention beamformer in the look direction
                
                for N = 1:N
                    [x, S, n] = ArrayMeasurementPlusNoiseGenerator(SNR_dB,psi,L); % Generate array measurements
                    R(:,:,N) = x*x'; % Signal + Noise Correlation matrix for sample N
                    Rs(:,:,N) = S*S'; % Signal Correlation matrix for sample N
                    Rn(:,:,N) = n*n'; % Noise Correlation matrix for sample N
                    y(N) = wc'*x; % Array output
                end
                
                R = mean(R,3); % Signal + Noise Correlation matrix
                Rs = mean(Rs,3); % Signal Correlation matrix
                Rn = mean(Rn,3); % Noise Correlation matrix
                
                PS = real(wc'*Rs*wc); % Signal Power
                Pn = real(wc'*Rn*wc); % Noise Power
                
                P = mean(abs(y).^2); % Power through array output
                P = real(wc'*R*wc); % Power through correlation matrix
                
                SINRout = PS/Pn; % Absolute
                SINRout_dB = 10*log10(SINRout); % [dB] Decibel
                
                Gcalc = L;
                Gestm = (1/SNR)/Pn;
                
                tit1 = {'Signal + Noise Correlation matrix';['N=' num2str(N) ', SNR=' num2str(SNR_dB) ' | Calc. Gain= ' num2str(Gcalc) ' | Est. Gain= ' num2str(Gestm)]};
                tit2 = {'Noise Correlation matrix';['N=' num2str(N) ', SNR=' num2str(SNR_dB) ' | Calc. Gain= ' num2str(Gcalc) ' | Est. Gain= ' num2str(Gestm)]};
                
                figure
                imagesc(abs(R))
                title(tit1)
                colorbar
                figure
                imagesc(abs(Rn))
                colorbar
                title(tit2)
                %%
            end
        end
    end
end