%% Conventional Beamformer
%% Background
% <html><h3>Conventional Beamformer</h3></html>
%
% Weight vector: $$\textbf{w}_c=\frac{1}{L}\textbf{S}_0$
%
% <html><h3>Correlation Matrix</h3></html>
%
% $$R = \langle\textbf{x}(t)\textbf{x}(t)^H\rangle$
%
% <html><h3>Array output</h3></html>
%
% $$y(t)=\textbf{w}_c^H\textbf{x}(t)$
%
%
% <html><h3>Simulation Algorithm</h3></html>
%
% # Generate signal for L element array
% # Use Coventional Beamformer weight vector
% # Calculate the array output y
% # Calculate output SNR
% # Calculate array Gain
%
%% Simulation Parameters
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
%% Simulation
%
convbeam_snr
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
%% Simulation
%
convbeam_sinr
%% References
% [1] L. Godara, Smart antennas. Boca Raton: CRC Press, 2004.