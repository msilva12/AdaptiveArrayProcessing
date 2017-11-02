%% Conventional Beamformer
%%
% <html><h3>Array and Source Parameters</h3></html>
%
% * Number of array elements: $$L=16$
% * Element spacing: $$d=\frac{\lambda}{2}m$
% * Propagation speed: $$c=3*10^8\frac{m}{s}$
% * Signal to Noise Ratio: $$SNR=[3, 10]dB$
% * Signal power: $$ps=1$
% * Noise power: $$\sigma_n^2=[0.5, 0.1]$
% * Number of samples: $$N=[10,100]$
% * Source direction: $$\theta_k=30^o$
% * Steering Vector in the Look direction: $$\textbf{S}_0=[exp({j2{\pi}f_0{\tau}_l}),\ldots,exp({j2{\pi}f_0{\tau}_L})]^T$
%
% Where $$\tau_l(\theta_k)=\frac{d}{c}(l-1)sin(\theta_k)$
%
% <html><h3>Conventional Beamformer</h3></html>
%
% Weight vector: $$\textbf{w}_c=\frac{1}{L}\textbf{S}_0$
%
% <html><h3>Correlation Matrice</h3></html>
%
% $$R = \langle\textbf{x}(t)\textbf{x}(t)^H\rangle$
%
% <html><h3>Array output</h3></html>
%
% $$y(t)=\textbf{w}_c^H\textbf{x}(t)$
%
% <html><h3>Simulation Algorithm</h3></html>
%
% # Generate signal for L element array
% # Use Coventional Beamformer weight vector
% # Calculate the array output y
% # Calculate output SNR
% # Calculate array Gain
%%
format compact;
%% Signal + Noise Only
%
% <html><h3>Array input</h3></html>
%
% $$\textbf{x}(t)=exp({j2{\pi}f_0t})\textbf{S}_0+\textbf{n}(t)$$
%
% <html><h3>Correlation Matrix</h3></html>
%
% $$R=ps*\textbf{S}_0\textbf{S}_0^H+{\sigma}_n^2{I}$$
%
% <html><h3>Signal to Noise Ratio</h3></html>
%
% $$SNR=\frac{ps}{\sigma_n^2}L$$
%
% <html><h3>Array Gain</h3></html>
%
% $$G=L$$
%
% <html><h3>Simulations</h3></html>
%
clear; close all; clc;
for L = 16 % Number of elements
    for N = [10 100] % Number of samples
        for SNR_dB = [3 10] % [dB] Signal to Noise Ratio
            for th_s = 30 % [deg] Source direction from normal of the array
                close all;
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
                
                SNRout = PS/Pn; % Absolute
                SNRout_dB = 10*log10(SNRout); % [dB] Decibel
                
                Gcalc = L;
                Gestm = SNRout/SNR;
                
                tit = {'Noise Correlation matrix';['N=' num2str(N) ', SNR=' num2str(SNR_dB) ' | Calc. Gain= ' num2str(Gcalc) ' | Est. Gain= ' num2str(Gestm)]};
                %                 fprintf('Simulation for:\nN=%-4i | SNR=%-4i | Calc. Gain=%-4.2f | Est. Gain=%-4.2f\n',N,SNR_dB,Gcalc,Gestm)
                
                figure
                imagesc(abs(R))
                title(tit)
                colorbar
                figure
                imagesc(abs(Rn))
                colorbar
                title(tit)
                %%
            end
        end
    end
end
%% Signal + Interference + Noise
%
% <html><h3>Array input</h3></html>
%
% $$\textbf{x}(t)=exp({j2{\pi}f_0t})\textbf{S}_0+pI*exp({j2{\pi}f_0t})\textbf{S}_I+\textbf{n}(t)$$
%
% <html><h3>Correlation Matrix</h3></html>
%
% $$R=ps*\textbf{S}_0\textbf{S}_0^H+pI*\textbf{S}_I\textbf{S}_I^H+{\sigma}_n^2{I}$$
%
% <html><h3>Signal to Noise Ratio</h3></html>
%
% $$SNR=ps/(pI(1-\rho)+\frac{\sigma_n^2}{L})$$
%
% <html><h3>Array Gain</h3></html>
%
% $$G=\frac{SNR_{out}}{SNR_{in}}$$
%
% <html><h3>Simulations</h3></html>
%
clear; close all; clc
for L = 16 % Number of elements
    for N = [10 100] % Number of samples
        for SNR_dB = [3 10] % [dB] Signal to Noise Ratio
            for th_s = 30 % [deg] Source direction from normal of the array
                for th_i = 60 % [deg] Interference direction from normal of the array
                    for pI = [0.1 1] % Interference Amplitude
                        close all;
                        %%
                        SNR = 10^(SNR_dB/10); % Absolute Signal to Noise Ratio
                        psi_s = pi*sind(th_s); % Phase difference between elements
                        psi_i = pi*sind(th_i);
                        SI = exp(-1j*psi_i*((1:L)-1)');
                        S0 = exp(-1j*psi_s*((1:L)-1)');
                        wc = (1/L)*S0; % convention beamformer in the look direction
                        
                        for N = 1:N
                            [~, S, n] = ArrayMeasurementPlusNoiseGenerator(SNR_dB,psi_s,L); % Generate array measurements
                            [~, I, ~] = ArrayMeasurementPlusNoiseGenerator(SNR_dB,psi_i,L); % Generate array measurements
                            I = exp(-1j*2*pi*rand).*I;
                            x = S + pI*I + n;
                            R(:,:,N) = x*x'; % Signal + Noise Correlation matrix for sample N
                            Rs(:,:,N) = S*S'; % Signal Correlation matrix for sample N
                            RI(:,:,N) = I*I'; % Interference Correlation matrix for sample N
                            Rn(:,:,N) = n*n' + pI*SI*SI'; % Noise Correlation matrix for sample N
                            y(N) = wc'*x; % Array output
                        end
                        
                        SI = mean(SI,3);
                        
                        R = mean(R,3); % Signal + Noise Correlation matrix
                        Rs = mean(Rs,3); % Signal Correlation matrix
                        RI = pI*SI*SI'; % Interference Correlation matrix
                        Rn = mean(Rn,3); % Noise Correlation matrix
                        
                        rho = 1-((S0'*SI)*(SI'*S0))/(L^2);
                        
                        PS = real(wc'*Rs*wc); % Signal Power
                        PI = real(wc'*RI*wc); % Interference Power
                        Pn = real(wc'*Rn*wc); % Noise Power
                        
                        P = mean(abs(y).^2); % Power through array output
                        P = real(wc'*R*wc); % Power through correlation matrix
                        
                        SNRout = PS/(PI+Pn); % Absolute
                        SNRout_dB = 10*log10(SNRout); % [dB] Decibel
                        
                        Gcalc = (PS/(pI*(1-rho)+(1/SNR)/L))/SNR;
                        Gestm = SNRout/SNR;
                        
                        tit = {'Noise Correlation matrix';['N=' num2str(N) ', SNR=' num2str(SNR_dB) ' | pI= ' num2str(pI) ' | Calc. Gain= ' num2str(Gcalc) ' | Est. Gain= ' num2str(Gestm)]};
                        
                        figure
                        imagesc(abs(R))
                        title(tit)
                        colorbar
                        figure
                        imagesc(abs(Rn))
                        colorbar
                        title(tit)
                        %%
                    end
                end
            end
        end
    end
end

%% References
% [1] L. Godara, Smart antennas. Boca Raton: CRC Press, 2004.