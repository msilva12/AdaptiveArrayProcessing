format compact; clear; close all; clc
for L = 16 % Number of elements
    for N = [10 100 10000] % Number of samples
        for SNR_dB = [3 10] % [dB] Signal to Noise Ratio
            for th_s = 30 % [deg] Source direction from normal of the array
                for th_i = 60 % [deg] Interference direction from normal of the array
                    for pI = [0.1 1] % Interference Amplitude
                        %%
                        SNR = 10^(SNR_dB/10); % Absolute Signal to Noise Ratio
                        psi_s = pi*sind(th_s); % Phase difference between elements
                        psi_i = pi*sind(th_i);
                        SI = exp(-1j*psi_i*((1:L)-1)');
                        S0 = exp(-1j*psi_s*((1:L)-1)');
                        w = ((eye(L)/SNR)^-1)*S0/(S0'*((eye(L)/SNR)^-1)*S0); % convention beamformer in the look direction

                        for N = 1:N
                            [~, S, n] = ArrayMeasurementPlusNoiseGenerator(SNR_dB,psi_s,L); % Generate array measurements
                            [~, I, ~] = ArrayMeasurementPlusNoiseGenerator(SNR_dB,psi_i,L); % Generate array measurements
                            I = pI*exp(-1j*2*pi*rand).*I;
                            x = S + I + n;
                            R(:,:,N) = x*x'; % Signal + Noise Correlation matrix for sample N
                            Rs(:,:,N) = S*S'; % Signal Correlation matrix for sample N
                            RI(:,:,N) = I*I'; % Interference Correlation matrix for sample N
                            Rn(:,:,N) = n*n'; % Noise Correlation matrix for sample N
                            y(N) = w'*x; % Array output
                        end
                        
                        SI = mean(SI,3);
                        
                        R = mean(R,3); % Signal + Noise Correlation matrix
                        Rs = mean(Rs,3); % Signal Correlation matrix
                        RI = mean(RI,3); % Interference Correlation matrix
                        Rn = mean(Rn,3); % Noise Correlation matrix
                        
                        rho = 1-((S0'*SI)*(SI'*S0))/(L^2);
                        
                        PS_cal = real(w'*S0*S0'*w);
                        PS = real(w'*Rs*w); % Signal Power
                        
                        PI_cal = real((0.1)^2*w'*SI*SI'*w);
                        PI = real(w'*RI*w); % Interference Power
                        
                        Pn_cal = real((1/SNR)*w'*w);
                        Pn = real(w'*Rn*w); % Noise Power
                        
                        P_cal = PS_cal + Pn_cal + PI_cal;
                        P = mean(abs(y).^2); % Power through array output
                        P = real(w'*R*w); % Power through correlation matrix
                        
                        PN = PI+Pn;
                        
                        SINRin_calc  = 1/((pI^2)*(1-rho)+(1/SNR)); % Signal to Interference + Noise Ratio
                        SINRout_calc = 1/((pI^2)*(1-rho)+(1/SNR)/L);
                        
                        SINRout = PS/PN; % Absolute
                        SINRout_dB = 10*log10(SINRout); % [dB] Decibel
                        
                        Gcalc = ((pI^2)*L/(1/SNR))*(1+(1/SNR)/pI^2)*(rho+(1/SNR)/(pI^2*L))/(1+(1/SNR)/(pI^2*L));
                        Gestm = (pI^2+1/SNR)/PN;
                        
                        tit1 = {'Signal + Interference + Noise Correlation Matrix';['N=' num2str(N) ', SNR=' num2str(SNR_dB) ' | pI= ' num2str(pI) ' | Calc. Gain= ' num2str(Gcalc) ' | Est. Gain= ' num2str(Gestm)]};
                        tit2 = {'Noise Correlation matrix';['N=' num2str(N) ', SNR=' num2str(SNR_dB) ' | pI= ' num2str(pI) ' | Calc. Gain= ' num2str(Gcalc) ' | Est. Gain= ' num2str(Gestm)]};
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
    end
end